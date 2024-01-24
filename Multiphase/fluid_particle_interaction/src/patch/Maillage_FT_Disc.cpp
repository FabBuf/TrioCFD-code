/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#include <Maillage_FT_Disc.h>
#include <TRUST_Deriv.h>
#include <TRUSTVect.h>
#include <Domaine.h>
#include <Domaine_VF.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Motcle.h>
#include <Statistiques.h>
#include <EcritureLectureSpecial.h>
#include <Champ_base.h>
#include <Paroi_FT_disc.h>
#include <Sauvegarde_Reprise_Maillage_FT.h>
#include <communications.h>
#include <Comm_Group.h>
#include <SFichier.h>
#include <LecFicDistribue.h>
#include <Probleme_FT_Disc_gen.h>
#include <Dirichlet_entree_fluide_leaves.h>
#include <Dirichlet_homogene.h>
#include <Debog.h>
#include <Array_tools.h>
#include <Param.h>
#include <stat_counters.h>

#if TCL_MODEL
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Schema_Temps_base.h>
#endif
//#define PATCH_HYSTERESIS_V2
//#define PATCH_HYSTERESIS_V3
#include <ArrOfBit.h>
#include <Domaine_VDF.h>
#include <Domaine.h>
// //#define DEBUG_HYSTERESIS_V2
/*
 * define permettant de post-traiter les triangles et leurs elements miroirs autour d'un point sommet s0.
 * Pour les tracer :
\rm xsplit_*
sed    -e '/Postraitement_ft_lata/d' -e '/^First postprocessing/,/First postprocessing/{//!b};d' err > bb
sed -i -e '/^Cas de 3 sommets/,/of exiting...$/d' \
       -e '/^Misma/,/of exiting$/d' -e '/^$/d' -e '/^3[ \t]*$/d' -e 's/TAG/# TAG/' bb
split -d -l 5 bb xsplit_
gnuplot -p << EOF
list=system("ls xsplit_*"); splot for [file in list] file w lp t file
EOF
 *
 */

Implemente_instanciable_sans_constructeur(Maillage_FT_Disc_Data_Cache,"Maillage_FT_Disc_Data_Cache",Objet_U);
Entree& Maillage_FT_Disc_Data_Cache::readOn(Entree& is)
{
  assert(0);
  exit();
  return is;
}
Sortie& Maillage_FT_Disc_Data_Cache::printOn(Sortie& os) const
{
  assert(0);
  exit();
  return os;
}
Maillage_FT_Disc_Data_Cache::Maillage_FT_Disc_Data_Cache()
{
  clear();
}
void Maillage_FT_Disc_Data_Cache::clear()
{
  tag_surface_ = -1;
  tag_normale_ = -1;
  tag_courbure_ = -1;
  surface_facettes_.resize_array(0);
  normale_facettes_.resize(0,0);
  courbure_sommets_.resize_array(0);
}

/////////////////////////////////////////////////////////////////////////

Implemente_instanciable_sans_constructeur(Maillage_FT_Disc,"Maillage_FT_Disc",Ensemble_Lagrange_base);


/*! @brief Pour chaque sommet du maillage, s'il est sur un bord, on calcule costheta min et max (hysteresis) correspondant a la condition aux limites ou
 *
 *   se trouve le sommet.
 *   L'angle est constant par face de bord... possibilite de faire mieux
 *   pour un champ xyz
 *
 */
void Maillage_FT_Disc::calculer_costheta_minmax(DoubleTab& costheta) const
{
  const Equation_base& eq = equation_transport();
  const Domaine_Cl_dis_base& domaine_cl = eq.domaine_Cl_dis().valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF,eq.domaine_dis().valeur());

  const long nb_som = nb_sommets();
  costheta.resize(nb_som, 2);
  costheta=0;
  const double deg_to_rad = M_PI / 180.;

  DoubleVect tmp;

  for (long i = 0; i < nb_som; i++)
    {
      // Indice de la face de bord dans le domaine_VF
      const long num_face = sommet_face_bord_[i];
      if ((num_face < 0) || (num_face>=domaine_vf.nb_faces_bord()))
        continue;
      const Cond_lim& cl =
        domaine_cl.la_cl_de_la_face(num_face);
      const Cond_lim_base& cl_base = cl.valeur();

      double theta1 = 0., theta2 = 0.;
      if (sub_type(Paroi_FT_disc, cl_base))
        {
          const Paroi_FT_disc& cl_ft =
            ref_cast(Paroi_FT_disc, cl_base);
          // Indice de la face sur la frontiere
          const Frontiere& front         = cl_ft.frontiere_dis().frontiere();
          const long num_premiere_face  = front.num_premiere_face();
          const long num_face_frontiere = num_face - num_premiere_face;
          // Valeur d'angle imposee sur la face
          const Paroi_FT_disc::Type_modele type_cl = cl_ft.get_type_modele();
          switch(type_cl)
            {
            case Paroi_FT_disc::CONSTANT:
              {
                const Champ_front_base& champ_frontiere = cl_ft.champ_front().valeur();
                champ_frontiere.valeurs_face(num_face_frontiere, tmp);
                theta1 = tmp[0];
                theta2 = theta1;
                break;
              }
            case Paroi_FT_disc::HYSTERESIS:
              {
                const Champ_front_base& champ_frontiere = cl_ft.champ_front().valeur();
                champ_frontiere.valeurs_face(num_face_frontiere, tmp);
                theta1 = tmp[0];
                theta2 = tmp[1];
                break;
              }
            case Paroi_FT_disc::SYMETRIE:
              {
                theta1 = 90.;
                theta2 = 90.;
                break;
              }
            default:
              Cerr << "Erreur dans Maillage_FT_Disc::calculer_costheta_minmax: condition aux limites non implementee"
                   << finl;
              exit();
            }
        }
      else
        {
          Cerr << "Error in Maillage_FT_Disc::calculer_costheta_minmax:\n "
               << cl_base.que_suis_je() << " is not a valid boundary condition type."
               << finl;
          exit();
        }
      costheta(i, 0) = cos(theta1 * deg_to_rad);
      costheta(i, 1) = cos(theta2 * deg_to_rad);
    }
  desc_sommets_.echange_espace_virtuel(costheta);
}

/*! @brief renvoie l'angle solide qui sert a calculer les surfaces et les volumes en bidim_axi
 *
 */
double Maillage_FT_Disc::angle_bidim_axi()
{
  return M_PI * 2.;
}

/*! @brief constructeur par defaut.
 *
 */
Maillage_FT_Disc::Maillage_FT_Disc() :
  statut_(RESET),
  mesh_state_tag_(0),
  temps_physique_(0.),
  is_solid_particle_(0),
  niveau_plot_(-1),
  correction_contact_courbure_coeff_(2.),
  calcul_courbure_iterations_(2),
  niter_pre_lissage_(0),
  methode_calcul_courbure_contact_line_(STANDARD),
  weight_CL_(0.5)
{
  mesh_data_cache_.typer("Maillage_FT_Disc_Data_Cache");
}

/*! @brief renvoie une ref non const au cache de valeurs calculees sur le maillage (courbure, surface, normale, .
 *
 * ..)
 *   On fait un cast de l'objet en non const (voir commentaire sur
 *   mesh_data_cache_ dans Maillage_FT_Disc.h)
 *
 */
Maillage_FT_Disc_Data_Cache& Maillage_FT_Disc::mesh_data_cache() const
{
  const Maillage_FT_Disc_Data_Cache *ptr = &mesh_data_cache_.valeur();
  // Cast en non const ici !!!
  return *(Maillage_FT_Disc_Data_Cache *) ptr;
}

/*! @brief Cette methode change le statut du maillage et invalide le cache de valeurs calculees (surface, courbure, .
 *
 * ..)
 *   Il faut l'appeler chaque fois que le maillage est modifie
 *   et avant le prochain appel a get_update_xxx.
 *   Cette methode est appelee par les methode publiques non constantes
 *   de Maillage_FT_Disc et par les methodes de Remaillage, ...
 *
 */
void Maillage_FT_Disc::maillage_modifie(Statut_Maillage nouveau_statut)
{
  //Process::Journal()<<"maillage_modifie de "<<statut_<<" a "<<nouveau_statut<<finl;
  if (nouveau_statut < PARCOURU && statut_ >= PARCOURU)
    {
      intersections_elem_facettes_.reset();
      intersections_face_facettes_x_.reset();
      intersections_face_facettes_y_.reset();
      intersections_face_facettes_z_.reset();
    }
  statut_ = nouveau_statut;
  mesh_state_tag_++;
}

/*! @brief La construction par copie est interdite !
 *
 */
Maillage_FT_Disc::Maillage_FT_Disc(const Maillage_FT_Disc&): Ensemble_Lagrange_base()
{
  assert(0);
  exit();
}

/*! @brief L'operateur = est interdit !
 *
 */
const Maillage_FT_Disc& Maillage_FT_Disc::operator=(const Maillage_FT_Disc&)
{
  assert(0);
  exit();
  return *this;
}

//Lecture des informations necessaires a l initialisation des coordonnees
//des points qui constituent un ensemble Lagrangien
// - soit lecture dans un fichier (remplissage direct de sommets_)
// - soit lecture des informations concernant la distribution dans les sous domaines :
//     nombre de sous domaines a considerer
//  Pour chaque sous domaine :
//  lecture de son nom
//   si repartition aleatoire des points
//        lecture du nombre de marqueurs
//   si repartition uniforme des points
//      lecture du nombre de marqueurs dans chaque direction

Entree& Maillage_FT_Disc::readOn(Entree& is)
{
  Motcles les_mots(3);
  {
    les_mots[0] = "fichier";
    les_mots[1] = "sous_domaines";
    les_mots[2] = "sous_zones";
  }

  Motcle motlu, accolade_fermee="}", accolade_ouverte="{";
  is >> motlu;
  if(motlu!=accolade_ouverte)
    {
      Cerr << "On attendait une { a la lecture d'une " << que_suis_je() << finl;
      Cerr << "et non : " << motlu << finl;
      exit();
    }
  is >> motlu;

  while (motlu != accolade_fermee)
    {
      long rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 :
          {
            Nom nomfic;
            is >> nomfic;

            EFichier fic;
            Cerr << "Lecture de l ensemble des sommets du Maillage_FT_Disc dans le fichier " << nomfic << finl;
            if (!fic.ouvrir(nomfic))
              {
                Cerr << " Erreur a l'ouverture du fichier." << finl;
                exit();
              }
            fic >> sommets_lu_;

            break;
          }
        case 1 :
        case 2 :
          {
            long nb_sz;
            is>>nb_sz;
            nom_sz.dimensionner(nb_sz);
            nb_marqs_sz.resize(nb_sz);
            long dim = Objet_U::dimension;
            nb_marqs_par_dir.resize(nb_sz,dim);

            for (long i=0; i<nb_sz; i++)
              {
                is>>nom_sz[i];
                is>>motlu;
                if (motlu=="aleatoire")
                  is>>nb_marqs_sz(i);
                else if (motlu=="uniforme")
                  {
                    nb_marqs_sz(i) = 1;
                    for (long k=0; k<dim; k++)
                      {
                        is >> nb_marqs_par_dir(i,k);
                        if (nb_marqs_par_dir(i,k)<=0)
                          {
                            if (je_suis_maitre())
                              Cerr<<"Le nombre de marqueurs specifies dans chacune des "<<dim<<" directions doit etre strictement positif"<<finl;
                            exit();
                          }
                        nb_marqs_sz(i) *=  nb_marqs_par_dir(i,k);
                      }
                  }
                else
                  {
                    Cerr<<"Les seules options disponibles pour la distribution de marqueurs sont aleatoire et uniforme"<<finl;
                    Cerr<<motlu<<" n est pas reconnu"<<finl;
                    exit();
                  }

              }

            break;
          }
        default :
          {
            if (je_suis_maitre())
              {
                Cerr << "On ne comprend pas le mot : " << motlu << " dans " << que_suis_je() << finl;
                Cerr << "Les mots compris sont : " << finl;
                for (long i=0; i<les_mots.size(); i++)
                  Cerr<<les_mots[i]<<finl;
              }
            exit();
          }
        }
      is >> motlu;
    }

  return is;
}

Sortie& Maillage_FT_Disc::printOn(Sortie& os) const
{
  Cerr << "Erreur : ::printOn n'est pas code." << finl;
  assert(0);
  return os;
}


Sortie& Maillage_FT_Disc::printFa7(long fa7,long affsom, Sortie& os) const
{
  const ArrOfDouble& surface_facette_ = get_update_surface_facettes();
  const DoubleTab& normale_facette_ = get_update_normale_facettes();
  long isom,som,nbsom = facettes_.dimension(1);
  os<<"#fa7="<<fa7<<"  soms= ";
  for (isom=0 ; isom<nbsom ; isom++)
    {
      som = facettes_(fa7,isom);
      os<<"  "<<som<<"("<<sommet_PE_owner_[som]<<"-"<<sommet_num_owner_[som]<<")";
    }
  if (surface_facette_.size_array()==facettes_.dimension(1))
    {
      os<<"    surf= "<<surface_facette_[fa7]<<"  normale= "<<normale_facette_(fa7,0)<<" "<<normale_facette_(fa7,1);
      if (dimension==3)
        {
          os<<" "<<normale_facette_(fa7,2);
        }
    }
  os<<finl;
  if (affsom==1 && nbsom>2)
    {
      for (isom=0 ; isom<nbsom ; isom++)
        {
          som = facettes_(fa7,isom);
          printSom(som,os);
        }
      som = facettes_(fa7,0);
      printSom(som,os);
      os<<finl;
    }

  return os;
}
Sortie& Maillage_FT_Disc::printSom(long som,Sortie& os) const
{
  os<<"sommet "<<som;
  if (som>=0)
    {
      os<<"  pe-Som="<<sommet_PE_owner_[som]<<"-"<<sommet_num_owner_[som]<<"  coords= "<<sommets_(som,0)<<" "<<sommets_(som,1);
      if (dimension==3)
        {
          os<<" "<<sommets_(som,2);
        }
      os<<" elem="<<sommet_elem_[som]<<" face_x="<<sommet_face_(som,0)<<" face_y="<<sommet_face_(som,1)<<" face_z="<<sommet_face_(som,2)<<" faceB="<<sommet_face_bord_[som]; // EB : rajout sommet_face_[som]
      os<<finl;
    }

  return os;
}

void Maillage_FT_Disc::ecrire_plot(const Nom& nom,double un_temps, long niveau_requete) const
{
  if (niveau_requete>niveau_plot_)
    return;

  static long compteur_plot = 0;
  Nom nom_fic=Objet_U::nom_du_cas();
  nom_fic += "_";
  char str[14];
#ifndef INT_is_64_
  snprintf(str,14,"%03d",compteur_plot++);
#else
  snprintf(str,14,"%03ld",compteur_plot++);
#endif
  nom_fic += Nom(str);
  nom_fic += "_";
  if (nom!="")
    {
      nom_fic += nom;
      nom_fic += "_";
    }
  if (Process::nproc()>1)
    {
#ifndef INT_is_64_
      if (Process::nproc()<=1000)
        snprintf(str,14,"%03d_",me());
      else if (Process::nproc()<=10000)
        snprintf(str,14,"%04d_",me());
      else if (Process::nproc()<=100000)
        snprintf(str,14,"%05d_",me());
#else
      if (Process::nproc()<=1000)
        snprintf(str,14,"%03ld_",me());
      else if (Process::nproc()<=10000)
        snprintf(str,14,"%04ld_",me());
      else if (Process::nproc()<=100000)
        snprintf(str,14,"%05ld_",me());
#endif
      else
        {
          Cerr << "Error in Maillage_FT_Disc::ecrire_plot." << finl;
          Cerr << "Contact TRUST support." << finl;
          exit();
        }
      nom_fic += Nom(str);
    }
  nom_fic += (Nom(un_temps));
  nom_fic += ".plot";
  SFichier fic(nom_fic);
  fic.setf(std::ios::scientific);
  fic.precision(5);
  Process::Journal() << "ecriture de " << nom_fic << " au temps " << un_temps << finl;

  //balyage des facettes
  long fa7,isom,som, k;
  const long nbfacettes = facettes_.dimension(0);
  const long nb_som_par_facette = facettes_.dimension(1);
  FTd_vecteur3 cdg;
  const double coeff = 0.1;
  for (fa7=0 ; fa7<nbfacettes ; fa7++)
    {
      if (! facette_virtuelle(fa7))
        {
          fic << "# ";
          printFa7(fa7,0,fic);
          //calcul du cdg de la fa7
          for (k=0 ; k<dimension ; k++)
            {
              cdg[k] = 0.;
            }
          for (isom=0 ; isom<nb_som_par_facette ; isom++)
            {
              som = facettes_(fa7,isom);
              for (k=0 ; k<dimension ; k++)
                {
                  cdg[k] += sommets_(som,k);
                }
            }
          for (k=0 ; k<dimension ; k++)
            {
              cdg[k] /= nb_som_par_facette;
            }
          //impression des sommets de la facette, avec leger decalage vers le cdg
          for (isom=0 ; isom<nb_som_par_facette ; isom++)
            {
              som = facettes_(fa7,isom);
              for (k=0 ; k<dimension ; k++)
                {
                  fic << "   "<< sommets_(som,k);
                }
              for (k=0 ; k<dimension ; k++)
                {
                  fic << "   "<< sommets_(som,k) - coeff*(sommets_(som,k) - cdg[k]);
                }
              fic<<finl;
            }
          if (nb_som_par_facette>2)
            {
              isom = 0;
              som = facettes_(fa7,isom);
              for (k=0 ; k<dimension ; k++)
                {
                  fic << "   "<< sommets_(som,k);
                }
              for (k=0 ; k<dimension ; k++)
                {
                  fic << "   "<< sommets_(som,k) - coeff*(sommets_(som,k) - cdg[k]);
                }
              fic<<finl;
            }
          fic<<finl<<finl;
          //impression du cdg
          for (k=0 ; k<dimension ; k++)
            {
              fic << "   "<< cdg[k];
            }
          fic<<finl<<finl;
        }
    }

  fic.close();
}

/*! @brief Cette fonction permet de lire les parametres pour le maillage des interfaces
 *
 * @param (is) flot d'entree
 * @return (Entree) le flot d'entree
 */
Entree& Maillage_FT_Disc::lire_param_maillage(Entree& is)
{
  Cerr<<"Lecture des parametres de remaillage (Remaillage_FT::lire_param_remaillage)"<<finl;

  Param param("Maillage_FT_Disc::lire_param_maillage");
  param.ajouter("niveau_plot",&niveau_plot_);
  param.ajouter("correction_contact_courbure_coeff",&correction_contact_courbure_coeff_);
  param.ajouter("niter_pre_lissage",&niter_pre_lissage_);
  param.ajouter("calcul_courbure_iterations",&calcul_courbure_iterations_);
  param.ajouter("methode_calcul_courbure_contact_line", &methode_calcul_courbure_contact_line_);
  param.dictionnaire("standard", (long)STANDARD);
  param.dictionnaire("mirror", (long)MIRROR);
  param.dictionnaire("improved", (long)IMPROVED);
  param.dictionnaire("none", (long)NONE);
  param.dictionnaire("weighted", (long)WEIGHTED);
  param.dictionnaire("hysteresis", (long)HYSTERESIS);
  param.ajouter("weight_CL",&weight_CL_);
  param.lire_avec_accolades(is);

  return is;
}

/*! @brief on remplit refequation_transport_, schema_comm_domaine_ desc_sommets_.comm_group_ et desc_facettes_.comm_group_
 *
 * Precondition: le domaine_dis de l'equation doit etre complete (joints)
 */
void Maillage_FT_Disc::associer_equation_transport(const Equation_base& equation)
{
  const Transport_Interfaces_FT_Disc& eq = ref_cast(Transport_Interfaces_FT_Disc,equation);
  refequation_transport_ = eq;

  const Domaine_dis& domaine_dis = eq.domaine_dis();
  const Parcours_interface& parcours_interface = eq.parcours_interface();
  associer_domaine_dis_parcours(domaine_dis, parcours_interface);
}

void Maillage_FT_Disc::associer_domaine_dis_parcours(const Domaine_dis& domaine_dis, const Parcours_interface& parcours)
{
  refdomaine_dis_ = domaine_dis;
  refparcours_interface_ = parcours;

  // On recupere la liste des PE voisins
  ArrOfIntFT pe_list;
  for (const auto& itr : domaine_dis.domaine().faces_joint())
    {
      const Joint& joint = itr;
      const long pe_voisin = joint.PEvoisin();
      pe_list.append_array(pe_voisin);
    }
  // La liste des processeurs avec qui on communique dans ce schema sont tous
  // les voisins du maillage eulerien. Communications symetriques (on envoie
  // et on recoit a tous les processeurs voisins).
  schema_comm_domaine_.set_send_recv_pe_list(pe_list, pe_list);
}

/*! @brief vide toutes les structures du maillage, le statut passe a RESET.
 *
 */
void Maillage_FT_Disc::reset()
{
  sommets_.resize(0,dimension);
  facettes_.resize(0,dimension);
  voisins_.resize(0,dimension);
  sommet_elem_.resize_array(0);
  sommet_face_.resize(0,dimension); // EB /!\ On a 3 indicatrices aux faces. Un sommet appartient a 3 faces
  sommet_face_bord_.resize_array(0);
  sommet_PE_owner_.resize_array(0);
  sommet_num_owner_.resize_array(0);
  desc_sommets_.reset();
  desc_facettes_.reset();
  drapeaux_sommets_.resize_array(0);
  intersections_elem_facettes_.reset();
  intersections_face_facettes_x_.reset();
  intersections_face_facettes_y_.reset();
  intersections_face_facettes_z_.reset();
  mesh_data_cache().clear();

  maillage_modifie(RESET);
}

/*! @brief Recopie une partie du maillage source dans *this.
 *
 * Si niveau_copie == MINIMAL, seuls les membres de l'etat minimal sont copies.
 *
 */
void Maillage_FT_Disc::recopie(const Maillage_FT_Disc& source, Statut_Maillage niveau_copie)
{
  reset();
  if (niveau_copie == RESET)
    return; // C'est sans interet, mais en toute rigueur on traite le cas...

  // Copie des membres qui definissent l'etat minimal:
  refequation_transport_ = source.refequation_transport_;
  refdomaine_dis_ = source.refdomaine_dis_;
  refparcours_interface_ = source.refparcours_interface_;
  schema_comm_domaine_ = source.schema_comm_domaine_;
  temps_physique_ = source.temps_physique_;
  sommets_ = source.sommets_;
  facettes_ = source.facettes_;
  voisins_ = source.voisins_;
  sommet_elem_ = source.sommet_elem_;
  sommet_face_ = source.sommet_face_; // EB
  sommet_face_bord_ = source.sommet_face_bord_;
  sommet_PE_owner_ = source.sommet_PE_owner_;
  sommet_num_owner_ = source.sommet_num_owner_;
  facette_num_owner_ = source.facette_num_owner_;
  desc_sommets_ = source.desc_sommets_;
  desc_facettes_ = source.desc_facettes_;
  drapeaux_sommets_ = source.drapeaux_sommets_;

  maillage_modifie(MINIMAL);

  if (niveau_copie > MINIMAL)
    {
      Cerr << "Erreur dans Maillage_FT_Disc::recopie\n";
      Cerr << " le niveau de copie > MINIMAL n'est pas implemente" << finl;
      assert(0);
      exit();
    }
}

//Cette methode ajoute le maillage de l'interface passe en parametre
//Amene l'etat du maillage a MINIMAL
//skip_facettes = 1 dans le cas d une injection de particules
void Maillage_FT_Disc::ajouter_maillage(const Maillage_FT_Disc& maillage_tmp,long skip_facettes)
{
  assert(maillage_tmp.statut_ >= MINIMAL);
  assert(&maillage_tmp != this);
  maillage_tmp.check_mesh(1,0,skip_facettes);

  const long nb_sommets_tmp = maillage_tmp.nb_sommets();
  const long nb_facettes_tmp = maillage_tmp.nb_facettes();
  const DoubleTab& sommets_tmp = maillage_tmp.sommets();
  const IntTab& facettes_tmp = maillage_tmp.facettes();
  const ArrOfInt& sommet_elem_tmp = maillage_tmp.sommet_elem_;
  const IntTab& sommet_face_tmp = maillage_tmp.sommet_face_; // EB
  const ArrOfInt& sommet_face_bord_tmp = maillage_tmp.sommet_face_bord_;
  const ArrOfInt& sommet_PE_owner_tmp = maillage_tmp.sommet_PE_owner();
  //const ArrOfInt & sommet_num_owner_tmp = maillage_tmp.sommet_num_owner();
  const ArrOfInt& drapeaux_sommets_tmp = maillage_tmp.drapeaux_sommets_;

  const long nb_sommets_ = nb_sommets();
  const long nb_facettes_ = nb_facettes();
  const long nb_sommets_tot = nb_sommets_ + nb_sommets_tmp;
  const long nb_facettes_tot = nb_facettes_ + nb_facettes_tmp;
  const long nb_som_par_facette = facettes_tmp.dimension(1);

  //on redimensionne les tableaux
  sommets_.resize(nb_sommets_tot,dimension);
  facettes_.resize(nb_facettes_tot,nb_som_par_facette);
  sommet_elem_.resize_array(nb_sommets_tot);
  sommet_face_.resize(nb_sommets_tot,dimension); // EB /!\ On a 3 indicatrices aux faces. Un sommet appartient a 3 faces
  sommet_face_bord_.resize_array(nb_sommets_tot);
  sommet_PE_owner_.resize_array(nb_sommets_tot);
  sommet_num_owner_.resize_array(nb_sommets_tot);
  drapeaux_sommets_.resize_array(nb_sommets_tot);
  facette_num_owner_.resize_array(nb_facettes_tot);

  //on ajoute les facettes au maillage
  long fa7, fa7_tmp, isom,som,som_tmp, k;
  for (fa7_tmp=0 ; fa7_tmp<nb_facettes_tmp ; fa7_tmp++)
    {
      fa7 = fa7_tmp + nb_facettes_;
      for (isom=0 ; isom<nb_som_par_facette ; isom++)
        {
          facettes_(fa7,isom) = facettes_tmp(fa7_tmp,isom) + nb_sommets_;
        }
      facette_num_owner_[fa7] = fa7; // Puis echange esp.virt. a la fin
    }

  //on ajoute les sommets au maillage + leurs parametres
  for (som_tmp=0 ; som_tmp<nb_sommets_tmp ; som_tmp++)
    {
      som = som_tmp + nb_sommets_;
      for (k=0 ; k<dimension ; k++)
        {
          sommets_(som,k) = sommets_tmp(som_tmp,k);
        }
      sommet_elem_[som]      = sommet_elem_tmp[som_tmp];
      for (long dim=0; dim<dimension; dim++) sommet_face_(som,dim)      = sommet_face_tmp(som_tmp,dim); // EB
      sommet_face_bord_[som] = sommet_face_bord_tmp[som_tmp];
      sommet_PE_owner_[som]  = sommet_PE_owner_tmp[som_tmp];
      sommet_num_owner_[som] = som; // Puis echange esp.virt. a la fin
      drapeaux_sommets_[som] = drapeaux_sommets_tmp[som_tmp];
    }

  // Il faut maintenant mettre a jour les descripteurs
  // avec renumerotation des elements
  ArrOfInt elements_tmp;
  elements_tmp.set_smart_resize(1);
  {
    //descripteurs des sommets : espace distant
    const Desc_Structure_FT& desc_sommets_tmp =  maillage_tmp.desc_sommets();
    const Descripteur_FT& espace_tmp = desc_sommets_tmp.espace_distant();
    Descripteur_FT& espace_ = desc_sommets_.espace_distant();
    const ArrOfInt& pe_voisins_tmp = espace_tmp.pe_voisins();
    long ipe_tmp, pe_tmp, nb_pe_tmp = pe_voisins_tmp.size_array();
    for (ipe_tmp=0 ; ipe_tmp<nb_pe_tmp ; ipe_tmp++)
      {
        pe_tmp = pe_voisins_tmp[ipe_tmp];
        elements_tmp = espace_tmp.elements(pe_tmp);
        elements_tmp += nb_sommets_;
        espace_.ajoute_elements(pe_tmp,elements_tmp);
      }
    espace_.calcul_liste_pe_voisins();
  }
  {
    //descripteurs des sommets : espace virtuel
    const Desc_Structure_FT& desc_sommets_tmp =  maillage_tmp.desc_sommets();
    const Descripteur_FT& espace_tmp = desc_sommets_tmp.espace_virtuel();
    Descripteur_FT& espace_ = desc_sommets_.espace_virtuel();
    const ArrOfInt& pe_voisins_tmp = espace_tmp.pe_voisins();
    long ipe_tmp,pe_tmp, nb_pe_tmp = pe_voisins_tmp.size_array();
    for (ipe_tmp=0 ; ipe_tmp<nb_pe_tmp ; ipe_tmp++)
      {
        pe_tmp = pe_voisins_tmp[ipe_tmp];
        elements_tmp = espace_tmp.elements(pe_tmp);
        elements_tmp += nb_sommets_;
        espace_.ajoute_elements(pe_tmp,elements_tmp);
      }
    espace_.calcul_liste_pe_voisins();
  }
  {
    //descripteurs des facettes : espace distant
    const Desc_Structure_FT& desc_facettes_tmp =  maillage_tmp.desc_facettes();
    const Descripteur_FT& espace_tmp = desc_facettes_tmp.espace_distant();
    Descripteur_FT& espace_ = desc_facettes_.espace_distant();
    const ArrOfInt& pe_voisins_tmp = espace_tmp.pe_voisins();
    long ipe_tmp,pe_tmp, nb_pe_tmp = pe_voisins_tmp.size_array();
    for (ipe_tmp=0 ; ipe_tmp<nb_pe_tmp ; ipe_tmp++)
      {
        pe_tmp = pe_voisins_tmp[ipe_tmp];
        elements_tmp = espace_tmp.elements(pe_tmp);
        elements_tmp += nb_facettes_;
        espace_.ajoute_elements(pe_tmp,elements_tmp);
      }
    espace_.calcul_liste_pe_voisins();
  }
  {
    //descripteurs des facettes : espace virtuel
    const Desc_Structure_FT& desc_facettes_tmp =  maillage_tmp.desc_facettes();
    const Descripteur_FT& espace_tmp = desc_facettes_tmp.espace_virtuel();
    Descripteur_FT& espace_ = desc_facettes_.espace_virtuel();
    const ArrOfInt& pe_voisins_tmp = espace_tmp.pe_voisins();
    long ipe_tmp,pe_tmp, nb_pe_tmp = pe_voisins_tmp.size_array();
    for (ipe_tmp=0 ; ipe_tmp<nb_pe_tmp ; ipe_tmp++)
      {
        pe_tmp = pe_voisins_tmp[ipe_tmp];
        elements_tmp = espace_tmp.elements(pe_tmp);
        elements_tmp += nb_facettes_;
        espace_.ajoute_elements(pe_tmp,elements_tmp);
      }
    espace_.calcul_liste_pe_voisins();
  }
  //puis enfin mettre a jour leur schema de comm
  desc_sommets_.calcul_schema_comm(nb_sommets_tot);
  desc_facettes_.calcul_schema_comm(nb_facettes_tot);

  desc_sommets_.echange_espace_virtuel(sommet_num_owner_);
  desc_facettes_.echange_espace_virtuel(facette_num_owner_);

  maillage_modifie(MINIMAL);

  check_mesh();
}
/*! @brief Remplit la structure intersections_elem_facettes_.
 *
 * Le statut passe a PARCOURU.
 *
 * Precondition: statut >= MINIMAL
 */
void Maillage_FT_Disc::parcourir_maillage()
{
  if (statut_==RESET)
    {
      Cerr << "Error! Maillage_FT_Disc is empty. Contact TRUST support." << finl;
      Process::exit();
    }
  if (statut_ >= PARCOURU)
    return;
  // Remplit la structure intersections_elem_facettes et
  // ajoute des facettes sur les processeurs "pauvres"
  const Parcours_interface& p = refparcours_interface_.valeur();

  static const Stat_Counter_Id counter = statistiques().new_counter(3, "Parcours de l'interface", "FrontTracking");
  statistiques().begin_count(counter);
  p.parcourir(*this);
  statistiques().end_count(counter);
  maillage_modifie(PARCOURU);

}

// debut EB
void Maillage_FT_Disc::remplir_equation_plan_faces_aretes_internes(Domaine_dis& domaine_dis)
{
  Parcours_interface& p = refparcours_interface_.valeur();
  p.remplir_equation_plan_faces_aretes_internes(domaine_dis);
}

// fin EB

/*! @brief Complete les structures de donnees du maillage.
 *
 * Le statut passe a COMPLET.
 *
 * Precondition: statut >= MINIMAL
 */
void Maillage_FT_Disc::completer_maillage()
{
  assert(statut_ != RESET);
  parcourir_maillage();
  if (statut_ >= COMPLET)
    return;
  // Force le calcul de la normale et la surface des elements
  get_update_surface_facettes();
  get_update_normale_facettes();
  get_update_courbure_sommets();
  calculer_voisins();
  statut_ = COMPLET;
}


/*! @brief Calcul de la fonction indicatrice (on suppose que "indicatrice" a la structure d'un tableau de valeurs aux elements, on ne remplit
 *
 *  que les elements reels).
 *  La fraction volumique de la phase 1 dans les elements traverses par
 *  une interface est determinee a partir des donnees du parcours dans
 *   "intersections_elem_facettes_".
 *  Les autres elements sont remplis par une methode heuristique utilisant
 *  l'indicatrice_precedente.
 *
 * Precondition: statut >= PARCOURU
 *  Attention, l'algorithme est concu de sorte que l'on puisse utiliser le
 * meme tableau "indicatrice" et "indicatrice_precedente".
 */
void Maillage_FT_Disc::calcul_indicatrice(DoubleVect& indicatrice,
                                          const DoubleVect& indicatrice_precedente)
{
  assert(statut_ >= PARCOURU);

  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Calculer_Indicatrice", "FrontTracking");
  statistiques().begin_count(stat_counter);

  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine& ladomaine = domaine_dis.domaine();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());
  const long nb_elem = ladomaine.nb_elem();
  const long nb_elem_tot = ladomaine.nb_elem_tot();
  const IntTab& elem_faces = domaine_vf.elem_faces();
  const IntTab& face_voisins = domaine_vf.face_voisins();

  static ArrOfBit elements_calcules;
  elements_calcules.resize_array(nb_elem_tot);
  // On ne recalcule pas l'indicatrice sur la majorite du domaine,
  // uniquement les elements qui ne sont pas traverses et qui ont
  // une indicatrice qui n'est pas egale a 0 ou 1.
  elements_calcules = 1;

  indicatrice = indicatrice_precedente;
  // Mettre a zero les elements traverses, les elements voisins et les elements dont
  // l'indicatrice n'est ni a zero ni a un.
  {
    const long nb_elem_voisins = elem_faces.dimension(1);

    // Boucle sur les elements
    const ArrOfInt& index_elem = intersections_elem_facettes_.index_elem();
    assert(indicatrice.size() == nb_elem);
    long i;
    DoubleVect check(indicatrice);
    for (i = 0; i < nb_elem_tot; i++)
      {
        const double x = indicatrice_precedente[i];
        long check_voisins = ((x != 0.) && (x != 1.));
        if (i < nb_elem)
          {
            long index = index_elem[i];
            check_voisins |= (index >= 0);
            check(i) = check_voisins;
          }
      }
    check.echange_espace_virtuel();
    Debog::verifier("Maillage_FT_Disc::calcul_indicatrice check=",check);
    for (i = 0; i < nb_elem_tot; i++)
      {
        if (check(i))
          {
            elements_calcules.clearbit(i);
            // Boucle sur les voisins
            long j;
            for (j = 0; j < nb_elem_voisins; j++)
              {
                const long face = elem_faces(i, j);
                const long elem = face_voisins(face, 0) + face_voisins(face, 1) - i;
                if (elem >= 0 && elem < nb_elem_tot)
                  elements_calcules.clearbit(elem); // Voisin d'une interf => considere comme non calcule.
              }
          }
      }
  }

  // Ajout des contributions de volume
  // Les elements traverses par l'interface deviennent des elements_calcules
  {
    const ArrOfInt& index_elem =
      intersections_elem_facettes_.index_elem();
    assert(indicatrice.size() == nb_elem);

    // Boucle sur les elements
    for (long i = 0; i < nb_elem; i++)
      {

        long index = index_elem[i];
        double somme_contrib = 0.;
        // Boucle sur les facettes qui traversent cet element
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data = intersections_elem_facettes_.data_intersection(index);
            somme_contrib += data.contrib_volume_phase1_;
            index = data.index_facette_suivante_;
          };
        while (somme_contrib > 1.)
          somme_contrib -= 1.;
        while (somme_contrib < 0.)
          somme_contrib += 1.;
        if (somme_contrib > 0.) // Pour ne pas faire les elements pures
          {
            indicatrice[i] = somme_contrib;
            elements_calcules.setbit(i);
          }
      }
  }

  // Calcul de l'indicatrice au voisinage de l'interface a l'aide
  // de la fonction distance.
  // Il reste dans elements_calcules[i] == 0 les voisins de l'interface
  /*
  {
    const DoubleTab& distance = equation_transport().get_update_distance_interface().valeurs();
    long i;
    long error_count = 0;

    for (i = 0; i < nb_elem; i++)
      {

        if (elements_calcules[i] == 0)
          {
            double x = distance(i);
            // La distance a-t-elle ete calculee pour cet element ?
            if (x > -1e10)
              {
                double v = (x > 0.) ? 1. : 0.;
                indicatrice[i] = v;
              }
            else
              {
                // Probleme : un element a une indicatrice suspecte et
                // on ne peut pas l'evaluer avec la fonction distance
                // (augmenter le nombre d'iterations du calcul de distance ?)
                error_count++;
              }
          }
      }
    if (error_count)
      {
        Cerr << "[" << me() << "] calcul_indicatrice : error_count = " << error_count << finl;
      }
  }
  */
  indicatrice.echange_espace_virtuel();

  // Certains elements ont une indicatrice erronee (error_count).
  // Deuxieme correction pour tuer les elements isoles qui seraient faux.
  // Pour chaque element monophasique (non traverse par une interface)
  //  calculer la moyenne de l'indicatrice sur les elements monophasiques voisins,
  //  si moyenne >0.5, mettre a 1 sinon mettre a 0
  {
    // Pour que l'algo soit parallele, on met a jour a la fin et non au fur et a mesure
    // sinon le resultat depend de l'ordre de parcours des elements
    // Liste des elements a mettre a changer (colonne 0) et valeur a mettre (colonne 1)
    IntTab elems_to_change(0,2);
    elems_to_change.set_smart_resize(1);
    const long nb_faces_elem = elem_faces.line_size();
    const ArrOfInt& index_elem = intersections_elem_facettes_.index_elem();
    for (long elem = 0; elem < nb_elem; elem++)
      {
        // element non traverse par une interface ?
        if (index_elem[elem] >= 0)
          continue;

        double somme = 0.; //somme des indicatrices des elements monophasiques voisins
        long count = 0; //nombre d'elements monophasiques voisins
        for (long ivoisin = 0; ivoisin < nb_faces_elem; ivoisin++)
          {
            const long face = elem_faces(elem, ivoisin);
            const long elem_voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
            // element au bord du domaine ?
            if (elem_voisin < 0)
              continue;
            // element voisin non monophasique ?
            if (indicatrice[elem_voisin] != 0. && indicatrice[elem_voisin] != 1.)
              continue;
            // L'element voisin est monophasique
            somme += indicatrice[elem_voisin];
            count++;
          }
        // Bug fix Salim Hamidi 2019/25/02
        // Si count==1 l'algo etait considere non pertinent... pas toujours correction du bug
        // Correction testee en VDF mais deux cas test VEF ont fait des ecarts
        if(count == 1)
          {
            // TODO: A investiguer
            // Elem est une maille diphasique.
            // Pourquoi lui met-on la valeur de somme qui vaut ici indicatrice[elem_voisin]
            // qui est pure (indicatrice[elem_voisin] vaut 0 ou 1)
            elems_to_change.append_line(elem, (long)std::lrint(somme));
          }
        if (count > 1)
          {
            long indic;
            if (indicatrice[elem] == 1.)
              indic = 1;
            else if (indicatrice[elem] == 0.)
              indic = 0;
            else
              indic = -1;
            long new_indic = ((somme * 2) > count) ? 1 : 0;
            // on compare somme/count avec l'indicatrice
            if (indic != new_indic)
              {
                // L'indicatrice de cet element doit etre corrigee
                elems_to_change.append_line(elem, new_indic);
              }
          }
      }
    // Correction du tableau
    const long n = elems_to_change.dimension(0);
    for (long i = 0; i < n; i++)
      {
        const long elem = elems_to_change(i, 0);
        indicatrice[elem] = elems_to_change(i, 1);
      }
    if (n > 0)
      Journal() << "Calcul indicatrice: correction par voisinage de " << n << " elements" << finl;
  }
  indicatrice.echange_espace_virtuel();

  Debog::verifier("Maillage_FT_Disc::calcul_indicatrice indicatrice=",indicatrice);
  elements_calcules.resize_array(0);

  calcul_cg_fa7(); // EB // EB : a verifier si on doit vraiment mettre ca ici, semble inapproprie

  statistiques().end_count(stat_counter);
}
// debut EB
/*! @brief Calcul de la fonction indicatrice aux faces du maillage eulerien (on suppose que "indicatrice_face" a la structure d'un tableau de valeurs aux faces, on ne remplit
 *
 *  que les faces reelles). Pour les faces de joint appartenant a 2 processeur, chaque processeur calcul le taux de controle dans le demi volume de controle lui appartenant.
 *  On somme ensuite les contributions pour avoir le taux de presence global au volume de controle de la face.
 *  La fraction volumique de la phase 1 dans les elements traverses par
 *  une interface est determinee a partir des donnees du parcours dans
 *   "intersections_face_facettes_".
 *  Les autres faces sont remplies par une methode heuristique utilisant
 *  l'indicatrice_precedente.
 *
 * Precondition: statut >= PARCOURU
 *  Attention, l'algorithme est concu de sorte que l'on puisse utiliser le
 * meme tableau "indicatrice" et "indicatrice_precedente".
 */void Maillage_FT_Disc::calcul_indicatrice_face(const DoubleVect& indicatrice, DoubleVect& indicatrice_face,
                                                  const DoubleVect& indicatrice_face_precedente)
{
  Cerr << "Maillage_FT_Disc::calcul_indicatrice_face"<< finl;

  assert(statut_ >= PARCOURU);
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Calculer_Indicatrice_Face", "FrontTracking");
  statistiques().begin_count(stat_counter);
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());
  const long nb_face = domaine_vf.nb_faces();
  const long nb_face_tot = domaine_vf.nb_faces_tot();
  const IntTab& elem_faces = domaine_vf.elem_faces();
  static ArrOfBit faces_calculees;
  faces_calculees.resize_array(nb_face_tot);
  // On ne recalcule pas l'indicatrice sur la majorite du domaine,
  // uniquement les faces qui ne sont pas traversees et qui ont
  // une indicatrice qui n'est pas egale a 0 ou 1.
  // on recalcule egalement les faces_doubles
  const DoubleVect& volumes_entrelaces=domaine_vf.volumes_entrelaces();
  faces_calculees = 1;
  const ArrOfInt& faces_doubles = domaine_vf.faces_doubles();

  indicatrice_face = indicatrice_face_precedente;

  // Mettre a zero les faces traverses, les faces voisines et les faces dont
  // l'indicatrice n'est ni a zero ni a un.
  {
    const long nb_faces_voisines = elem_faces.dimension(1); // idem que pour les elements
    // Boucle sur les faces
    // Si l'indicatrice precedente est differente de 0 ou 1 ou que la face est traversee --> alors on recalcule
    const ArrOfInt& index_face_x = intersections_face_facettes_x_.index_face();
    const ArrOfInt& index_face_y = intersections_face_facettes_y_.index_face();
    const ArrOfInt& index_face_z = intersections_face_facettes_z_.index_face();
    assert(indicatrice_face.size() == nb_face);
    long i;
    DoubleVect check(indicatrice_face);
    for (i = 0; i < nb_face_tot; i++)
      {
        const double x = indicatrice_face_precedente[i];
        long check_voisins = ( ((x != 0.) && (x != 1.)));
        if (i < nb_face)
          {
            long orientation=domaine_vf.orientation(i);
            long index;
            if (orientation==0) index=index_face_x[i];
            else if (orientation==1) index=index_face_y[i];
            else if (orientation==2) index=index_face_z[i];
            else
              exit();
            check_voisins |= (index >= 0);
            check(i) = check_voisins;

          }
        // if (faces_doubles(i)) check(i)=1; // on recalcule automatiquement pour les faces doubles
      }

    //MD_Vector_tools::echange_espace_virtuel(check,MD_Vector_tools::EV_MAX);
    check.echange_espace_virtuel();
    Debog::verifier("Maillage_FT_Disc::calcul_indicatrice_face check=",check);

    // On identifie les faces voisines pour lesquelles on recalculera l'indicatrice
    for (i = 0; i < nb_face_tot; i++)
      {
        if (check(i))
          {
            faces_calculees.clearbit(i);
            const long elem0=domaine_vf.face_voisins(i,0);
            const long elem1=domaine_vf.face_voisins(i,1);
            const long ori=domaine_vf.orientation(i);
            // Boucle sur les voisins
            long j;
            for (j = 0; j < nb_faces_voisines; j++)
              {
                // On cherche l'element voisin
                const long direction = (j>=dimension);
                long num_elem;
                long face_voisine;
                if (j%dimension==ori)
                  {
                    num_elem = (direction==0) ? elem0 : elem1;
                    face_voisine = (num_elem>=0) ? domaine_vf.elem_faces(num_elem, ori+direction*dimension) : -1;
                  }
                else
                  {
                    const long elem_2 = (elem0>0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem0,j),direction) : -1;
                    const long elem_3 = (elem1>0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem1,j),direction) : -1;
                    num_elem = std::max(elem_2,elem_3);
                    // S'il n'y a pas de voisin, on est au bord du domaine
                    if (num_elem >= 0)
                      face_voisine =  (elem_2>0) ? domaine_vf.elem_faces(elem_2, ori+dimension) : domaine_vf.elem_faces(elem_3, ori);
                    else
                      face_voisine = -1;
                  }

                if (face_voisine >= 0)
                  faces_calculees.clearbit(face_voisine);

              }
          }
      }
  }

  // Ajout des contributions de volume
  {
    const ArrOfInt& index_face_x =
      intersections_face_facettes_x_.index_face();
    const ArrOfInt& index_face_y =
      intersections_face_facettes_y_.index_face();
    const ArrOfInt& index_face_z =
      intersections_face_facettes_z_.index_face();

    assert(indicatrice_face.size() == nb_face);
    // Boucle sur les faces
    for (long i = 0; i < nb_face; i++)
      {
        long orientation=domaine_vf.orientation(i);
        long index;
        if (orientation==0)
          {
            // Faces de normale x
            index = index_face_x[i];
            double somme_contrib = 0.;
            // Boucle sur les facettes qui traversent cette face
            while (index >= 0)
              {
                const Intersections_Face_Facettes_Data& data = intersections_face_facettes_x_.data_intersection(index);
                somme_contrib += data.contrib_volume_phase1_;
                index = data.index_facette_suivante_;
              };
            while (somme_contrib > 1.)
              somme_contrib -= 1.;
            while (somme_contrib < 0.)
              somme_contrib += 1.;
            if (somme_contrib > 0.)
              {
                indicatrice_face[i] = somme_contrib;
                faces_calculees.setbit(i);
              }
          }
        else if (orientation==1)
          {
            // Faces de normale y
            index = index_face_y[i];
            double somme_contrib = 0.;
            // Boucle sur les facettes qui traversent cette face
            while (index >= 0)
              {
                const Intersections_Face_Facettes_Data& data = intersections_face_facettes_y_.data_intersection(index);
                somme_contrib += data.contrib_volume_phase1_;
                index = data.index_facette_suivante_;
              };
            while (somme_contrib > 1.)
              somme_contrib -= 1.;
            while (somme_contrib < 0.)
              somme_contrib += 1.;
            if (somme_contrib > 0.)
              {
                indicatrice_face[i] = somme_contrib;
                faces_calculees.setbit(i);
              }
          }
        else if (orientation==2)
          {
            // Faces de normale z
            index = index_face_z[i];
            double somme_contrib = 0.;
            // Boucle sur les facettes qui traversent cette face
            while (index >= 0)
              {
                const Intersections_Face_Facettes_Data& data = intersections_face_facettes_z_.data_intersection(index);
                somme_contrib += data.contrib_volume_phase1_;
                index = data.index_facette_suivante_;
              };
            while (somme_contrib > 1.)
              somme_contrib -= 1.;
            while (somme_contrib < 0.)
              somme_contrib += 1.;
            if (somme_contrib > 0.)
              {
                indicatrice_face[i] = somme_contrib;
                faces_calculees.setbit(i);
              }
          }
        else
          exit(); // On a rien a faire la
      }
  }


  // Calcul de l'indicatrice au voisinage de l'interface a l'aide
  // de la fonction distance.

  {
    const DoubleTab& distance = equation_transport().get_update_distance_interface_faces().valeurs();
    long i;
    long error_count = 0;

    for (i = 0; i < nb_face; i++)
      {

        if (faces_calculees[i] == 0)
          {
            double x = distance(i);
            // La distance a-t-elle ete calculee pour cette face ?
            if (x > -1e10)
              {
                double v = (x > 0.) ? 1. : 0.;
                indicatrice_face[i] = v;
              }
            else
              {
                // Probleme : une face a une indicatrice suspecte et
                // on ne peut pas l'evaluer avec la fonction distance
                // (augmenter le nombre d'iterations du calcul de distance ?)
                error_count++;
              }
          }
        // if (faces_doubles(i) && fabs(domaine_vf.xv(i,0))<0.3e-3 && domaine_vf.xv(i,1)<7.8e-3 && domaine_vf.xv(i,1)>7.2e-3 &&  domaine_vf.xv(i,2)<-1.2e-3 ) Cerr << "Face double - indic " << indicatrice_face(i) << "\t" << domaine_vf.xv(i,0) << " " << domaine_vf.xv(i,1) << " " << domaine_vf.xv(i,2) <<  finl;
      }
    if (error_count)
      {
        Cerr << "[" << me() << "] calcul_indicatrice_face : error_count = " << error_count << finl;
      }
  }


  //indicatrice_face.echange_espace_virtuel();
  long n_proc=Process::nproc();
  if (n_proc>1)
    {
      for (long face=0; face<domaine_vf.nb_faces_tot(); face++)
        {
          //double face_monophasique = (indicatrice_face(face) !=0 || indicatrice_face(face)!=1) ? 0.5:1;

          double coeff=0.5;
          if (faces_doubles(face) && face<nb_face)
            {
              indicatrice_face(face) *=volumes_entrelaces(face)*coeff;
            }
          if (face>=nb_face) indicatrice_face(face)=0;
        }
      MD_Vector_tools::echange_espace_virtuel(indicatrice_face,MD_Vector_tools::EV_SOMME_ECHANGE);

      for (long face=0; face<nb_face; face++)
        {
          if (faces_doubles(face))
            {
              indicatrice_face(face) /=volumes_entrelaces(face);
            }
        }
    }
  indicatrice_face.echange_espace_virtuel();

  // Certaines faces ont une indicatrice erronee (error_count).
  // Deuxieme correction pour tuer les faces isolees qui seraient fausses.
  // Pour chaque face monophasique (non traversee par une interface)
  //  calculer la moyenne de l'indicatrice sur les faces monophasiques voisins,
  //  si moyenne >0.5, mettre a 1 sinon mettre a 0
  {
    // Pour que l'algo soit parallele, on met a jour a la fin et non au fur et a mesure
    // sinon le resultat depend de l'ordre de parcours des faces
    // Liste des faces a mettre a changer (colonne 0) et valeur a mettre (colonne 1)
    IntTab faces_to_change(0,2);
    faces_to_change.set_smart_resize(1);
    const long nb_faces_face = elem_faces.line_size();
    const ArrOfInt& index_face_x = intersections_face_facettes_x_.index_face();
    const ArrOfInt& index_face_y = intersections_face_facettes_y_.index_face();
    const ArrOfInt& index_face_z = intersections_face_facettes_z_.index_face();
    for (long face = 0; face < nb_face; face++)
      {
        // face non traversee par une interface ?
        long orientation = domaine_vf.orientation(face);
        if (orientation==0)
          {
            if (index_face_x[face] >= 0)
              continue;
          }
        else if (orientation==1)
          {
            if (index_face_y[face] >= 0)
              continue;
          }
        else if (orientation==2)
          {
            if (index_face_z[face] >= 0)
              continue;
          }
        else
          exit();

        double somme = 0.; //somme des indicatrices des faces monophasiques voisins
        long count = 0; //nombre de faces monophasiques voisines
        const long elem0=domaine_vf.face_voisins(face,0);
        const long elem1=domaine_vf.face_voisins(face,1);
        const long ori=domaine_vf.orientation(face);

        for (long ivoisin = 0; ivoisin < nb_faces_face; ivoisin++)
          {
            // On cherche l'element voisin
            const long direction = (ivoisin>=dimension);
            long num_elem;
            long face_voisine;

            if (ivoisin%dimension==ori)
              {
                num_elem = (direction==0) ? elem0 : elem1;
                face_voisine = (num_elem>=0) ? domaine_vf.elem_faces(num_elem, ori+direction*dimension) : -1;
              }
            else
              {
                const long elem_2 = (elem0>0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem0,ivoisin),direction) : -1;
                const long elem_3 = (elem1>0) ? domaine_vf.face_voisins(domaine_vf.elem_faces(elem1,ivoisin),direction) : -1;
                num_elem = std::max(elem_2,elem_3);
                // S'il n'y a pas de voisin, on est au bord du domaine
                if (num_elem >= 0)
                  face_voisine =  (elem_2>0) ? domaine_vf.elem_faces(elem_2, ori+dimension) : domaine_vf.elem_faces(elem_3, ori);
                else
                  face_voisine = -1;
              }
            // face au bord du domaine ?
            if (face_voisine < 0)
              continue;
            // face voisine non monophasique ?
            if (indicatrice_face[face_voisine] != 0. && indicatrice_face[face_voisine] != 1.)
              continue;
            // La face voisine est monophasique
            somme += indicatrice_face[face_voisine];
            count++;
          }
        // Si count==1 l'algo etait considere non pertinent... pas toujours correction du bug
        // Correction testee en VDF mais deux cas test VEF ont fait des ecarts
        if(count == 1)
          {
            faces_to_change.append_line(face, (long)std::lrint(somme));
          }
        if (count > 1)
          {
            long indic;
            if (indicatrice_face[face] == 1.)
              indic = 1;
            else if (indicatrice_face[face] == 0.)
              indic = 0;
            else
              indic = -1;
            long new_indic = ((somme * 2) > count) ? 1 : 0;
            // on compare somme/count avec l'indicatrice
            if (indic != new_indic)
              {
                // L'indicatrice de cet element doit etre corrigee
                faces_to_change.append_line(face, new_indic);
              }
          }
      }

    // Correction du tableau
    const long n = faces_to_change.dimension(0);
    for (long i = 0; i < n; i++)
      {
        const long face = faces_to_change(i, 0);
        if (!faces_doubles(face)) indicatrice_face[face] = faces_to_change(i, 1);
      }
    if (n > 0)
      Journal() << "Calcul indicatrice face : correction par voisinage de " << n << " faces" << finl;
  }


  indicatrice_face.echange_espace_virtuel();

  Debog::verifier("Maillage_FT_Disc::calcul_indicatrice_face indicatrice_face=",indicatrice_face);
  faces_calculees.resize_array(0);

  statistiques().end_count(stat_counter);
}

// EB idem mais aux aretes
/*! @brief Calcul de la fonction indicatrice aux aretes du maillage eulerien (on suppose que "indicatrice_arete" a la structure d'un tableau de valeurs aux faces, on ne remplit
 *
 *  que les aretes reelles). Pour les aretes de joint appartenant a 4 (au max) processeurs, chaque processeur calcul le taux de controle dans le demi volume de controle lui appartenant.
 *  On somme ensuite les contributions pour avoir le taux de presence global au volume de controle de l'arete.
 *  ATTENTION : ce calcul n'est valable en parallele que pour les partitionnement "bien ordonnes". C'est a dire, que l'on peut definir le nombre de procs du domaine par Npx * Npy *Npz avec Npx, Npy, Npz, le nombre de procs suivant x,y,z.
 *  Avec Npx, constant suivant y et z, Npy constant suivant x et z, Npz constant suivant x et y.
 *  La fraction volumique de la phase 1 dans les elements traverses par
 *  une interface est determinee a partir des donnees du parcours dans
 *   "intersections_arete_facettes_".
 *  Les autres aretes sont remplies par une methode heuristique utilisant
 *  l'indicatrice_precedente.
 *
 * Precondition: statut >= PARCOURU
 *  Attention, l'algorithme est concu de sorte que l'on puisse utiliser le
 * meme tableau "indicatrice" et "indicatrice_precedente".
*/
void Maillage_FT_Disc::calcul_indicatrice_arete(const DoubleVect& indicatrice, DoubleVect& indicatrice_arete,
                                                const DoubleVect& indicatrice_arete_precedente)
{
  Cerr << "Maillage_FT_Disc::calcul_indicatrice_arete"<< finl;
  assert(statut_ >= PARCOURU);
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Calculer_Indicatrice_Arete", "FrontTracking");
  statistiques().begin_count(stat_counter);
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis.valeur());
  const long nb_aretes_reelles = domaine_vdf.nb_aretes_reelles();
  const IntVect& orientation_aretes=domaine_vdf.orientation_aretes();
  const IntTab& Elem_Aretes=domaine_vdf.domaine().elem_aretes();
  const IntTab& Qdm=domaine_vdf.Qdm();
  const IntVect& type_arete=domaine_vdf.type_arete();
  static ArrOfBit aretes_calculees;
  const DoubleVect& volumes_aretes=domaine_vdf.volumes_aretes();
  //const ArrOfInt& aretes_multiples = domaine_vdf.aretes_multiples();
  aretes_calculees.resize_array(nb_aretes_reelles);

  indicatrice_arete = indicatrice_arete_precedente;

  {
    const long nb_aretes_voisines = domaine_vdf.elem_faces().dimension(1);
    // Boucle sur les faces
    // Si l'indicatrice precedente est differente de 0 ou 1 ou que la face est traversee --> alors on recalcule
    const ArrOfInt& index_arete_x = intersections_arete_facettes_x_.index_arete();
    const ArrOfInt& index_arete_y = intersections_arete_facettes_y_.index_arete();
    const ArrOfInt& index_arete_z = intersections_arete_facettes_z_.index_arete();

    assert(indicatrice_arete.size_array() == nb_aretes_reelles);
    long i;
    DoubleVect check(indicatrice_arete);
    for (i = 0; i < nb_aretes_reelles; i++)
      {
        if (type_arete(i)!=2) continue; // on ne calcule que pour les aretes_internes
        const double x = indicatrice_arete_precedente[i];
        long check_voisins = ( ((x != 0.) && (x != 1.)));

        long orientation=(dimension-1)-orientation_aretes(i);
        long index;
        if (orientation==0) index=index_arete_x[i];
        else if (orientation==1) index=index_arete_y[i];
        else if (orientation==2) index=index_arete_z[i];
        else
          exit();

        check_voisins |= (index >= 0);
        check(i) = check_voisins;

      }

    // On identifie les aretes voisines pour lesquelles on recalculera l'indicatrice
    for (i = 0; i < nb_aretes_reelles; i++)
      {
        if (type_arete(i)!=2) continue; // on ne calcule que pour les aretes_internes
        if (check(i))
          {
            aretes_calculees.clearbit(i);
            const long ori_arete=(dimension-1)-orientation_aretes(i);

            long face1=Qdm(i,0);
            long face2=Qdm(i,1);

            long elem1=domaine_vdf.face_voisins(face1,0);
            long elem2=domaine_vdf.face_voisins(face2,0);
            long elem3=domaine_vdf.face_voisins(face1,1);
            //long elem4=domaine_vdf.face_voisins(face2,1);

            long face1_av=domaine_vdf.elem_faces(elem1, ori_arete+dimension);
            long face1_arr=domaine_vdf.elem_faces(elem1, ori_arete);

            long elem_av=domaine_vdf.face_voisins(face1_av,1);
            long elem_arr=domaine_vdf.face_voisins(face1_arr,0);

            //Cerr << elem1+elem2+elem3+elem_av+elem_arr<< finl;

            // Boucle sur les voisins
            long j;
            for (j = 0; j < nb_aretes_voisines; j++)
              {
                // On cherche l'element voisin
                long arete_voisine;

                if (j==ori_arete) arete_voisine=Elem_Aretes(elem_av,orientation_aretes(i)); // voir numerotation Aretes dans Hexaedre.cpp
                else if (j==ori_arete+dimension) arete_voisine=Elem_Aretes(elem_arr,orientation_aretes(i));

                else if (ori_arete==0) // arete YZ
                  {
                    if (j==1) arete_voisine=Elem_Aretes(elem1,11);
                    else if (j==1+dimension) arete_voisine=Elem_Aretes(elem2,2);
                    else if (j==2) arete_voisine=Elem_Aretes(elem1,8);
                    else arete_voisine=Elem_Aretes(elem3,2);
                  }
                else if (ori_arete) // arete XZ
                  {
                    if (j==0) arete_voisine=Elem_Aretes(elem1,10); // voir numerotation Aretes dans Hexaedre.cpp
                    else if (j==0+dimension) arete_voisine=Elem_Aretes(elem2,1);
                    else if (j==2) arete_voisine=Elem_Aretes(elem1,7);
                    else arete_voisine=Elem_Aretes(elem3,1);
                  }
                else // arete XY
                  {
                    if (j==0) arete_voisine=Elem_Aretes(elem1,9); // voir numerotation Aretes dans Hexaedre.cpp
                    else if (j==0+dimension) arete_voisine=Elem_Aretes(elem2,0);
                    else if (j==1) arete_voisine=Elem_Aretes(elem1,6);
                    else arete_voisine=Elem_Aretes(elem3,0);
                  }

                if (arete_voisine >= 0 && arete_voisine<nb_aretes_reelles) aretes_calculees.clearbit(arete_voisine);

              }

          }
      }


  }

  // Ajout des contributions de volume
  {
    const ArrOfInt& index_arete_x =
      intersections_arete_facettes_x_.index_arete();
    const ArrOfInt& index_arete_y =
      intersections_arete_facettes_y_.index_arete();
    const ArrOfInt& index_arete_z =
      intersections_arete_facettes_z_.index_arete();

    assert(indicatrice_arete.size() == nb_aretes_reelles);
    // Boucle sur les aretes

    for (long i = 0; i < nb_aretes_reelles; i++)
      {
        if (type_arete(i)!=2) continue; // on ne calcule que pour les aretes_internes
        const long ori_arete=(dimension-1)-orientation_aretes(i);
        long index;
        if (ori_arete==0)
          {

            // Faces de normale x
            index = index_arete_x[i];
            double somme_contrib = 0.;
            // Boucle sur les facettes qui traversent cette face
            while (index >= 0)
              {
                const Intersections_Arete_Facettes_Data& data = intersections_arete_facettes_x_.data_intersection(index);
                somme_contrib += data.contrib_volume_phase1_;
                index = data.index_facette_suivante_;
              };
            while (somme_contrib > 1.)
              somme_contrib -= 1.;
            while (somme_contrib < 0.)
              somme_contrib += 1.;
            if (somme_contrib > 0.)
              {
                indicatrice_arete[i] = somme_contrib;
                aretes_calculees.setbit(i);
              }

          }
        else if (ori_arete==1)
          {

            // Faces de normale y
            index = index_arete_y[i];
            double somme_contrib = 0.;
            // Boucle sur les facettes qui traversent cette face
            while (index >= 0)
              {
                const Intersections_Arete_Facettes_Data& data = intersections_arete_facettes_y_.data_intersection(index);
                somme_contrib += data.contrib_volume_phase1_;
                index = data.index_facette_suivante_;
              };
            while (somme_contrib > 1.)
              somme_contrib -= 1.;
            while (somme_contrib < 0.)
              somme_contrib += 1.;
            if (somme_contrib > 0.)
              {
                indicatrice_arete[i] = somme_contrib;
                aretes_calculees.setbit(i);
              }
          }
        else if (ori_arete==2)
          {
            // Faces de normale z
            index = index_arete_z[i];
            double somme_contrib = 0.;
            // Boucle sur les facettes qui traversent cette face
            while (index >= 0)
              {
                const Intersections_Arete_Facettes_Data& data = intersections_arete_facettes_z_.data_intersection(index);
                somme_contrib += data.contrib_volume_phase1_;
                index = data.index_facette_suivante_;
              };
            while (somme_contrib > 1.)
              somme_contrib -= 1.;
            while (somme_contrib < 0.)
              somme_contrib += 1.;
            if (somme_contrib > 0.)
              {
                indicatrice_arete[i] = somme_contrib;
                aretes_calculees.setbit(i);
              }
          }
        else
          exit(); // On a rien a faire la
      }
  }

  // Calcul de l'indicatrice au voisinage de l'interface a l'aide
  // de la fonction distance.
  {
    const DoubleTab& distance = equation_transport().get_update_distance_interface_aretes();
    long i;
    long error_count = 0;

    for (i = 0; i < nb_aretes_reelles; i++)
      {
        if (type_arete(i)!=2) continue; // on ne calcule que pour les aretes_internes
        if (aretes_calculees[i] == 0)
          {
            double x = distance(i);
            // La distance a-t-elle ete calculee pour cette face ?
            if (x > -1e10)
              {
                double v = (x > 0.) ? 1. : 0.;
                indicatrice_arete[i] = v;
              }
            else
              {
                // Probleme : une face a une indicatrice suspecte et
                // on ne peut pas l'evaluer avec la fonction distance
                // (augmenter le nombre d'iterations du calcul de distance ?)
                error_count++;
              }
          }
        // if (faces_doubles(i) && fabs(domaine_vf.xv(i,0))<0.3e-3 && domaine_vf.xv(i,1)<7.8e-3 && domaine_vf.xv(i,1)>7.2e-3 &&  domaine_vf.xv(i,2)<-1.2e-3 ) Cerr << "Face double - indic " << indicatrice_face(i) << "\t" << domaine_vf.xv(i,0) << " " << domaine_vf.xv(i,1) << " " << domaine_vf.xv(i,2) <<  finl;
      }
    if (error_count)
      {
        Cerr << "[" << me() << "] calcul_indicatrice_arete : error_count = " << error_count << finl;
      }
  }
  /*
    for (long arete=0; arete<domaine_vdf.nb_aretes_tot(); arete++)
      {
        double coeff=0.5;

  	  if (aretes_multiples(arete)==1) coeff=0.5;
  	  else if (aretes_multiples(arete)==2) coeff=1/3;
  	  else if (aretes_multiples(arete)==3) coeff=0.25;


        if ((aretes_multiples(arete)>0) && arete<nb_aretes_reelles)
          {
            //if (indicatrice_arete(arete)<1 && orientation_aretes(arete)==0) Cerr << "indic " << indicatrice_arete(arete) << "\tcg "
            //                                                                     << domaine_vdf.xa(arete,0) << " " << domaine_vdf.xa(arete,1) << " " << domaine_vdf.xa(arete,2) << finl;
            indicatrice_arete(arete) *=volumes_aretes(arete)*coeff;
          }
        if (arete>=nb_aretes_reelles) indicatrice_arete(arete)=0;
      }

    MD_Vector_tools::echange_espace_virtuel(indicatrice_arete,MD_Vector_tools::EV_SOMME_ECHANGE);

    for (long arete=0; arete<nb_aretes_reelles; arete++)
      {
        if (aretes_multiples(arete)>0)
          {
            indicatrice_arete(arete) /=volumes_aretes(arete);
          }
      }
  */
  // Certaines aretes ont une indicatrice erronee (error_count).
  // Deuxieme correction pour tuer les aretes isolees qui seraient fausses.
  // Pour chaque arete monophasique (non traversee par une interface)
  //  calculer la moyenne de l'indicatrice sur les aretes monophasiques voisins,
  //  si moyenne >0.5, mettre a 1 sinon mettre a 0
  {
    IntTab aretes_to_change(0,2);
    aretes_to_change.set_smart_resize(1);
    const long nb_faces_arete = 2*dimension;
    const ArrOfInt& index_arete_x = intersections_arete_facettes_x_.index_arete();
    const ArrOfInt& index_arete_y = intersections_arete_facettes_y_.index_arete();
    const ArrOfInt& index_arete_z = intersections_arete_facettes_z_.index_arete();

    for (long arete=0; arete<nb_aretes_reelles; arete++)
      {
        if (type_arete(arete)!=2) continue; // on ne calcule que pour les aretes_internes
        // face non traversee par une interface ?
        const long ori_arete=(dimension-1)-orientation_aretes(arete);
        if (ori_arete==0)
          {
            if (index_arete_x[arete] >= 0)
              continue;
          }
        else if (ori_arete==1)
          {
            if (index_arete_y[arete] >= 0)
              continue;
          }
        else if (ori_arete==2)
          {
            if (index_arete_z[arete] >= 0)
              continue;
          }
        else
          exit();

        double somme = 0.; //somme des indicatrices des aretes monophasiques voisines
        long count = 0; //nombre d'aretes monophasiques voisines


        long face1=Qdm(arete,0);
        long face2=Qdm(arete,1);
        //long face3=Qdm(arete,2);
        //long face4=Qdm(arete,3);

        long elem1=domaine_vdf.face_voisins(face1,0);
        long elem2=domaine_vdf.face_voisins(face2,0);
        long elem3=domaine_vdf.face_voisins(face1,1);
        //long elem4=domaine_vdf.face_voisins(face2,1);

        long face1_av=domaine_vdf.elem_faces(elem1, ori_arete+dimension);
        long face1_arr=domaine_vdf.elem_faces(elem1, ori_arete);

        long elem_av=domaine_vdf.face_voisins(face1_av,1);
        long elem_arr=domaine_vdf.face_voisins(face1_arr,0);


        for (long ivoisin=0; ivoisin<nb_faces_arete; ivoisin++)
          {
            long arete_voisine;
            if (ivoisin==ori_arete) arete_voisine= (elem_av>=0) ? Elem_Aretes(elem_av,orientation_aretes(arete)) : -1; // voir numerotation Aretes dans Hexaedre.cpp
            else if (ivoisin==ori_arete+dimension) arete_voisine=(elem_arr>=0) ? Elem_Aretes(elem_arr,orientation_aretes(arete)):-1;

            else if (ori_arete==0) // arete YZ
              {
                if (ivoisin==1) arete_voisine= (elem1>=0) ? Elem_Aretes(elem1,11) : -1;
                else if (ivoisin==1+dimension) arete_voisine=(elem2>=0) ? Elem_Aretes(elem2,2) : -1;
                else if (ivoisin==2) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,8) : -1;
                else arete_voisine=(elem3>=0) ? Elem_Aretes(elem3,2) : -1;
              }
            else if (ori_arete) // arete XZ
              {
                if (ivoisin==0) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,10) : -1; // voir numerotation Aretes dans Hexaedre.cpp
                else if (ivoisin==0+dimension) arete_voisine=(elem2>=0) ? Elem_Aretes(elem2,1) : -1;
                else if (ivoisin==2) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,7) : -1;
                else arete_voisine=(elem3>=0) ? Elem_Aretes(elem3,1) : -1;
              }
            else // arete XY
              {
                if (ivoisin==0) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,9) : -1; // voir numerotation Aretes dans Hexaedre.cpp
                else if (ivoisin==0+dimension) arete_voisine=(elem2>=0) ? Elem_Aretes(elem2,0) : -1;
                else if (ivoisin==1) arete_voisine=(elem1>=0) ? Elem_Aretes(elem1,6) : -1;
                else arete_voisine=(elem3>=0) ? Elem_Aretes(elem3,0) : -1;
              }

            // face au bord du domaine ?
            if (arete_voisine < 0 || arete_voisine>nb_aretes_reelles)
              continue;
            // face voisine non monophasique ?
            if (indicatrice_arete[arete_voisine] != 0. && indicatrice_arete[arete_voisine] != 1.)
              continue;


            somme += indicatrice_arete[arete_voisine];
            count++;
          }
        // Si count==1 l'algo etait considere non pertinent... pas toujours correction du bug
        // Correction testee en VDF mais deux cas test VEF ont fait des ecarts
        if(count == 1)
          {
            aretes_to_change.append_line(arete, (long)std::lrint(somme));
          }
        if (count > 1)
          {
            long indic;
            if (indicatrice_arete[arete] == 1.)
              indic = 1;
            else if (indicatrice_arete[arete] == 0.)
              indic = 0;
            else
              indic = -1;
            long new_indic = ((somme * 2) > count) ? 1 : 0;
            // on compare somme/count avec l'indicatrice
            if (indic != new_indic)
              {
                // L'indicatrice de cet element doit etre corrigee
                aretes_to_change.append_line(arete, new_indic);
              }
          }
      }

    // Correction du tableau
    const long n = aretes_to_change.dimension(0);
    for (long i = 0; i < n; i++)
      {
        const long arete = aretes_to_change(i, 0);
        if (1) indicatrice_arete[arete] = aretes_to_change(i, 1); // if !faces_doubles(face)
      }
    if (n > 0)
      Journal() << "Calcul indicatrice arete : correction par voisinage de " << n << " aretes" << finl;
  }

  Debog::verifier("Maillage_FT_Disc::calcul_indicatrice_arete indicatrice_arete=",indicatrice_arete);
  aretes_calculees.resize_array(0);
  statistiques().end_count(stat_counter);

}
// fin EB
/*! @brief Deplace les sommets de l'interface d'un vecteur "deplacement" fourni, Change eventuellement les sommets de processeur, cree eventuellement
 *
 *   des lignes de contact et detecte les collisions.
 *
 * @param (deplacement) un tableau de taille (nb_sommets(), Objet_U::dimension) contenant le vecteur deplacement de chaque sommet. Les deplacements des sommets virtuels sont ignores.
 */
void Maillage_FT_Disc::transporter(const DoubleTab& deplacement)
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Transporter_maillage", "FrontTracking");
  statistiques().begin_count(stat_counter);

  assert(deplacement.dimension(0) == sommets_.dimension(0));
  assert(deplacement.dimension(1) == sommets_.dimension(1));

  long i;
  const long dim = Objet_U::dimension;
  const long nb_som = nb_sommets();
  // La taille est surestimee
  ArrOfInt liste_sommets;
  liste_sommets.set_smart_resize(1);
  liste_sommets.resize_array(nb_som, Array_base::NOCOPY_NOINIT);
  // Le tableau de deplacement des sommets reels
  DoubleTab vecteur;
  vecteur.set_smart_resize(1);
  vecteur.resize(nb_som, dim, Array_base::NOCOPY_NOINIT);
  ArrOfIntFT liste_sommets_sortis;
  ArrOfIntFT numero_face_sortie;

  long n = 0;
  for (i = 0; i < nb_som; i++)
    {
      if (! sommet_virtuel(i))
        {
          liste_sommets[n] = i;
          for (long j = 0; j < dim; j++)
            vecteur(n, j) = deplacement(i, j);
          n++;
        }
    }
  // Taille finale
  liste_sommets.resize_array(n);
  vecteur.resize_dim0(n);

  // Cette methode s'occupe simplement des sommets :
  deplacer_sommets(liste_sommets,
                   vecteur,
                   liste_sommets_sortis,
                   numero_face_sortie);

  // Lorsque le premier sommet d'une facette est recu ou envoyer,
  // il faut changer son proprietaire.
  corriger_proprietaires_facettes();

  maillage_modifie(MINIMAL);
  statistiques().end_count(stat_counter);
}

//-Remplit un tableau temporaire sommets_tmp de positions
//(voir Ensemble_Lagrange_base::remplir_sommets_tmp)
//-Genere la structure complete a partir de ce tableau de positions
//(voir remplir_structure())
void Maillage_FT_Disc::generer_structure()
{
  DoubleTab sommets_tmp;
  remplir_sommets_tmp(sommets_tmp);
  remplir_structure(sommets_tmp);
}

//On remplit la structure de l ensemble Lagrangien
// -sommets_ -sommet_elem_ -sommet_face_bord_
// -sommet_PE_owner_ -sommet_num_owner_ -drapeaux_sommets_

// Demarche similaire a celle appliquee dans Marching_Cubes::construire_iso
// -constrution de def_noeud
// -on remplit la structure de l ensemble Lagrangien :
// sommets_ sommet_elem_ sommet_face_bord_
// sommet_PE_owner_ sommet_num_owner_ drapeaux_sommets_

//som_init_util_ est rempli (=1 si le sommet i appartient au processeur, 0 sinon)

void Maillage_FT_Disc::remplir_structure(const DoubleTab& soms)
{

  //  Au depart la signification de def_noeud est la suivante:
  //  def_noeud(i,0) = 0 //pas de sens ici
  //  def_noeud(i,1) = 0 //pas de sens ici
  //  def_noeud(i,2) = numero_PE (-1 sinon)
  //  def_noeud(i,3) = numero de sommet
  //  def_noeud(i,4) = numero de l element qui contient le sommet (-1 sinon)
  Cerr << "Maillage_FT_Disc::remplir_structure " << finl;
  IntTab def_noeud(0, 4);
  def_noeud.set_smart_resize(1);

  reset();
  construire_noeuds(def_noeud,soms); // ici on a pas encore rempli def_noeud(i,5), def_noeud(i,6) et def_noeud(i,7)

  Descripteur_FT& espace_distant = desc_sommets_.espace_distant();
  Descripteur_FT& espace_virtuel = desc_sommets_.espace_virtuel();
  //long nb_som_tot=def_noeud.dimension(0);
  //def_noeud.resize(nb_som_tot,8);
  //  def_noeud(i,5) = numero de la face de normale x qui contient le sommet (-1 sinon) // EB
  //  def_noeud(i,6) = numero de la face de normale y qui contient le sommet (-1 sinon) // EB
  //  def_noeud(i,7) = numero de la face de normale z qui contient le sommet (-1 sinon) // EB
  // debut EB
  // Maintenant que l'on connait le num des elements qui contiennent les sommets,
  // on peut construire sommet_face de la maniere suivante :
  // --------------------------------
  //	|	 |    !	x |   !   |
  //	|---e1---|---e2---|---e3--|
  //	|	 |    !	  |   !   |
  // ---------------------------------
  //	| 1  |	2|	  |	   |	y
  //	|--------|--------|--------|	^
  //	| 3  |	4|	  |	   |    |
  // ---------------------------------	---> x
  // x : sommet lagrangien
  // representation du volume de controle d'une face eulerienne
  // !	|  !
  // !	|  !


  // Faces de normale x :
  // Si x appartient a la zone 2 ou 4, alors le sommet appartient a la face 2 de l'element e2
  // Si x appartient a la zone 1 ou 3, alors, le sommmet appartient a la face 0 de l'element e2
  // Faces de normale y :
  // Si x appartient a la zone 1 ou 2, alors le sommet appartient a la face 4 de l'element e2
  // Si x appartient a la zone 3 ou 4, alors, le sommmet appartient a la face 1 de l'element e2
  // Idem avec Z


  // Fin EB
  //const Zone_dis& zone_dis = refzone_dis_.valeur(); // EB
  //const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis.valeur());  // EB
  //const long nb_faces_reelles=zone_vf.nb_faces(); // EB

  long nb_noeuds = def_noeud.dimension(0);
  for (long noeud=0; noeud<nb_noeuds; noeud++)
    {
      def_noeud(noeud,3) =  def_noeud(noeud,4);
      def_noeud(noeud,4) = -1;
      // debut EB
      /*
      def_noeud(noeud,5) = -1;
      def_noeud(noeud,6) = -1;
      def_noeud(noeud,7) = -1;

      long mon_elem=def_noeud(noeud, 3);
      // debut EB
      long elem_virt=-1;
      IntVect check_som_in_face(2*dimension);
      if (mon_elem>=0) check_som_in_face=1;
      else // Le sommet ne m'appartient pas. Je ne connais pas l'element qui le contient -> je le cherche
        {
          elem_virt=zone_vf.zone().chercher_elements(sommets_(noeud,0),sommets_(noeud,1),sommets_(noeud,2)); // element eulerien contenant le sommet virtuel
          assert(elem_virt>=0);
          // Si le sommet est dans le volume de controle d'une face double, on met check_som_in_face a 1.
          check_som_in_face=0;
          for (long dim=0; dim<2*dimension; dim++) // on parcourt les faces du volume de controle de l'element. Si une des faces est reelles, alors le sommet est dans le volume de controle de la premiere couche de joint
            {
              long la_face=zone_vf.elem_faces(elem_virt,dim);
              if (la_face<nb_faces_reelles)  // la face est reelle mais l'element est virtuel -> face_double
                {
                  check_som_in_face(dim)=1;
                  mon_elem=elem_virt;
                }
            }
        }
      if (mon_elem>=0)
        {
          IntVect faces_elem_eulerien(2*dimension);
          for (long dim=0; dim<2*dimension; dim++) faces_elem_eulerien(dim) = zone_vf.elem_faces(mon_elem,dim); // on recupere les faces de l'element eulerien
          double pos_sommet, coord_elem;
          for (long dim=0; dim<dimension; dim++)
            {
              pos_sommet=sommets_(noeud,dim); // on recupere la position du sommet suivant l'axe "dim"
              coord_elem=zone_vf.xp(mon_elem,dim);
              if (pos_sommet < coord_elem && check_som_in_face(dim)) def_noeud(noeud,5+dim)=faces_elem_eulerien(dim);
              else if (pos_sommet >= coord_elem && check_som_in_face(dim+dimension)) def_noeud(noeud,5+dim)=faces_elem_eulerien(dim+dimension);
              assert(def_noeud(noeud,5+dim)>=0);
            }
        }
      */
      // fin EB
    }
  espace_distant.calcul_liste_pe_voisins();
  espace_virtuel.calcul_liste_pe_voisins();
  desc_sommets_.calcul_schema_comm(nb_noeuds);

  // Ici def_noeud a change de definition:
  //  def_noeud(i,0) = 0 //pas de sens ici
  //  def_noeud(i,1) = 0 //pas de sens ici
  //  def_noeud(i,2) = numero_PE (-1 sinon)
  //  def_noeud(i,3) = numero de l'element (-1 sinon)
  //  def_noeud(i,4) = -1 //Pas de noeaud considere sur une face de bord au depart
  //  def_noeud(i,5) = numero de la face de normale x (-1 sinon) // EB
  //  def_noeud(i,6) = numero de la face de normale y (-1 sinon) // EB
  //  def_noeud(i,7) = numero de la face de normale z (-1 sinon) // EB

  //Construction de som_init_util pour la structure construite
  som_init_util_.resize_array(nb_noeuds);
  som_init_util_ = 0;
  for (long noeud=0; noeud<nb_noeuds; noeud++)
    if (def_noeud(noeud,3)!=-1)
      som_init_util_[noeud] = 1;

  //On remplit sommets_
  //renum donne la correspondance entre un numero de sommet retenu
  //et son indice dans def_noeud
  IntTab renum;
  calculer_coord_noeuds(def_noeud,soms,renum);

  {
    long nbsommets = sommets_.dimension(0);
    sommet_PE_owner_. resize_array(nbsommets);
    sommet_num_owner_.resize_array(nbsommets);
    sommet_elem_.     resize_array(nbsommets);
    sommet_face_.     resize(nbsommets,dimension); // EB
    sommet_face_bord_.resize_array(nbsommets);
    drapeaux_sommets_.resize_array(nbsommets);
    desc_sommets_.calcul_schema_comm(nbsommets);
    desc_sommets_.remplir_element_pe(sommet_PE_owner_);
    // Le descripteur des facettes est vide, mais on calcule quand meme
    // le schema de comm pour qu'il soit valide

    for (long i = 0; i < nbsommets; i++)
      {
        sommet_num_owner_[i] = i;
        const long elem = def_noeud(renum(i), 3);
        const long face = def_noeud(renum(i), 4);
        /*const long face_x = def_noeud(renum(i), 5); // EB
        const long face_y = def_noeud(renum(i), 6); // EB
        const long face_z = def_noeud(renum(i), 7); // EB */
        sommet_elem_[i] = elem;
        sommet_face_bord_[i] = face;
        /*sommet_face_(i,0)=face_x; // EB
        //Cerr << "face_x : " << face << " cg : " << domaine_vdf.xv(face_x,0) << " " << domaine_vdf.xv(face_x,1) << " " << domaine_vdf.xv(face_x,2)<< finl;
        sommet_face_(i,1)=face_y; // EB
        sommet_face_(i,2)=face_z; // EB */
        drapeaux_sommets_[i] = 0;
      }
    desc_sommets_.echange_espace_virtuel(sommet_num_owner_);
    desc_sommets_.echange_espace_virtuel(sommet_face_bord_);
  }


  statut_ = Maillage_FT_Disc::MINIMAL;

  Journal() << "Maillage_FT_Disc::construire_points" << finl;
  Journal() << sommets_.dimension(0) << " noeuds, ";

  maillage_modifie(Maillage_FT_Disc::MINIMAL);
}


//On applique une procedure pour determiner
// -le processeur a qui appartient un sommet
// -le numero d element qui contient ce sommet
// -le numera de la face qui contient ce sommet // EB
//On remplit ensuite def_noeud

void Maillage_FT_Disc::construire_noeuds(IntTab& def_noeud,const DoubleTab& soms)
{
  long som;
  const Domaine& madomaine = mon_dom_.valeur();
  const long nb_elements_reels = madomaine.nb_elem();
  const long nb_sommets_tot = soms.dimension(0);

  //Procedure pour determiner a quel processeur appartiennent les sommets (tmp3)
  //et quel element les contiennent (tmp2)
  ///////////////////////////////////////////////////////////////////////////

  long nb_som_tot;
  const long dim = Objet_U::dimension;
  const long moi = Process::me();
  const long nbproc = Process::nproc();
  const long chunk_size = 65536*nbproc ; // Lire les sommets par blocs

  // Lecture des indices des elements contenant le sommet
  // long erreur_sommets_exterieurs = 0;

  DoubleTab tmp;
  // tmp2 : pour chaque sommet, -1 si je ne le garde pas, sinon numero de l'element qui contient le sommet
  ArrOfInt tmp2;
  // tmp3 : -1 ou numero du processeur qui garde le sommet.
  ArrOfInt tmp3;
  // tmp4 : pour chaque sommet, -1 si je ne le garde pas, sinon numero de la face qui contient le sommet // EB
  ArrOfInt tmp4; // EB
  long i = 0;

  nb_som_tot = soms.dimension(0);
  if (nb_som_tot>chunk_size)
    {
      Cerr << "The maximum of particules that you can use in your datafile is equal to " << chunk_size/nbproc << " / processor" << finl;
      Cerr << "If you want to exceed this limitation, you must parallelize your calculation with an apropriate number of processor" << finl;
      exit();
    }

  while (i < nb_som_tot)
    {
      const long n_to_read = std::min(chunk_size, nb_som_tot - i);
      tmp.resize(n_to_read, dim);
      tmp2.resize_array(n_to_read);
      tmp3.resize_array(n_to_read);
      tmp4.resize_array(n_to_read); // EB
      for (long j = 0; j < n_to_read; j++)
        for (long k = 0; k < dim; k++)
          tmp(j, k) = soms(i+j, k);

      madomaine.chercher_elements(tmp, tmp2);
      // Ne pas tenir compte des sommets dans les elements virtuels
      for (long j = 0; j < n_to_read; j++)
        if (tmp2[j] >= nb_elements_reels)
          tmp2[j] = -1;
      // Verifier qu'un processeur et un seul prend la main.
      // Le premier processeur initialise tmp3 en mettant son numero
      // pour les sommets qu'il possede et -1 pour les autres.
      // Puis la liste est passee au processeur suivant qui la complete.
      if (moi == 0)
        {
          // Preparer la premiere liste
          for (long j = 0; j < n_to_read; j++)
            if (tmp2[j] >= 0)
              tmp3[j] = moi;
            else
              tmp3[j] = -1;
        }
      else
        {
          // Recevoir une liste d'attribution du processeur precedent
          // et la completer
          recevoir(tmp3, moi - 1, moi, 0 /* canal */);
          for (long j = 0; j < n_to_read; j++)
            if (tmp3[j] >= 0) // Le sommet a ete pris par un autre processeur
              tmp2[j] = -1;
            else if (tmp2[j] >= 0) // Le sommet n'a pas ete pris et est chez moi
              tmp3[j] = moi;
        }
      if (moi < nbproc - 1)
        {
          envoyer(tmp3, moi, moi + 1, 0 /* canal */);
        }
      else
        {
          // Le dernier processeur prend les sommets qui restent
          for (long j = 0; j < n_to_read; j++)
            if (tmp3[j] < 0)
              {
                tmp3[j] = moi;
                // Pour l'instant on ne sait pas traiter des maillages qui depassent du
                // domaine:
                //            erreur_sommets_exterieurs=1;
              }
        }

      i += n_to_read;
    }

  ////////////////////////////////////////////////////////////////////////////
  def_noeud.resize(nb_sommets_tot, 5);
  long nb_noeuds = 0;

  for (som = 0; som < nb_sommets_tot; som++)
    {
      // Remplissage de def_noeud pour chaque point traites
      nb_noeuds=som;
      long PE_element;
      long numero_element_a_stocker;
      PE_element = tmp3[som];
      numero_element_a_stocker = tmp2[som];
      def_noeud(nb_noeuds, 0) = 0;
      def_noeud(nb_noeuds, 1) = 0;
      def_noeud(nb_noeuds, 2) = PE_element;
      def_noeud(nb_noeuds, 3) = nb_noeuds;
      def_noeud(nb_noeuds, 4) = numero_element_a_stocker;

    }

}

//On remplit sommets_ (coordonnees des points a suivre)
//On cree renum pour connaitre la correspondance
//entre nouveau numero de sommet et ancien dans def_noeud

void Maillage_FT_Disc::calculer_coord_noeuds(const IntTab& def_noeud,const DoubleTab& soms,IntTab& renum)
{
  // Raccourci vers les coordonnees des sommets du maillage eulerien
  const long nb_noeuds = def_noeud.dimension(0);
  DoubleTab& coord_noeuds = sommets_;
  coord_noeuds.resize(nb_noeuds, dimension);
  long noeud;
  long dimension3 = (dimension == 3);

  long indice=0;
  renum.resize(0);
  for (noeud = 0; noeud < nb_noeuds; noeud++)
    {
      if (def_noeud(noeud,3) != -1)
        {
          coord_noeuds(indice, 0) = soms(noeud,0);
          coord_noeuds(indice, 1) = soms(noeud,1);
          if (dimension3)
            coord_noeuds(indice, 2) = soms(noeud,2);
          renum.resize(indice+1);
          renum(indice)=noeud;
          indice++;
        }
    }

  if (!dimension3)
    coord_noeuds.resize(indice,2);
  else
    coord_noeuds.resize(indice,3);

}

//Methode qui assure le deplacement des sommets (particules)
//sans considerations sur les facettes
//Une fois  le deplacement effectue on applique une procedure de suppression des
//sommets virtuels et des sommets se trouvant sur une frontiere ouverte
void Maillage_FT_Disc::transporter_simple(const DoubleTab& deplacement)
{
  assert(deplacement.dimension(0) == sommets_.dimension(0));
  assert(deplacement.dimension(1) == sommets_.dimension(1));

  long i;
  const long dim = Objet_U::dimension;
  const long nb_som = nb_sommets();
  // La taille est surestimee
  ArrOfIntFT liste_sommets(nb_som);
  // Le tableau de deplacement des sommets reels
  DoubleTabFT vecteur(nb_som, dim);
  ArrOfIntFT liste_sommets_sortis;
  ArrOfIntFT numero_face_sortie;

  long n = 0;
  for (i = 0; i < nb_som; i++)
    {
      if (! sommet_virtuel(i))
        {
          liste_sommets[n] = i;
          for (long j = 0; j < dim; j++)
            vecteur(n, j) = deplacement(i, j);
          n++;
        }
    }
  // Taille finale
  liste_sommets.resize_array(n);
  vecteur.resize(n, dim);

  long skip_facettes = 1;
  deplacer_sommets(liste_sommets,
                   vecteur,
                   liste_sommets_sortis,
                   numero_face_sortie,
                   skip_facettes);

  maillage_modifie(MINIMAL);
}

//Description :
// Cette methode calcule les equations des plans passant par les aretes de la facette fa7
//   et perpendiculaires a la facette.
// Les 3 premiers coefficients (a,b,c) designent la normale a ce plan, orientee vers l'interieur de la facette.
long Maillage_FT_Disc::calcule_equation_plan_areteFa7(
  long fa7, long isom,
  double& a, double& b, double& c, double& d) const
{
  long res = 1;
  const long nb_faces = facettes_.dimension(1);
  const DoubleTab& normale_facettes = get_update_normale_facettes();

  long som0,som1;
  double x0, y0, z0, x1, y1, z1;
  double inverse_norme;
  //recuperation des sommets de l'arete
  som0 = facettes_(fa7,isom);
  som1 = facettes_(fa7,(isom+1)%nb_faces);
  //recuperation des coordonnees de l'arete
  x0 = sommets_(som0,0);
  y0 = sommets_(som0,1);
  z0 = sommets_(som0,2);
  x1 = sommets_(som1,0);
  y1 = sommets_(som1,1);
  z1 = sommets_(som1,2);
  //calcul de la normale a l'arete (dans le plan de la facette)
  a = (y1-y0) * normale_facettes(fa7,2) - (z1-z0) * normale_facettes(fa7,1);
  b = (z1-z0) * normale_facettes(fa7,0) - (x1-x0) * normale_facettes(fa7,2);
  c = (x1-x0) * normale_facettes(fa7,1) - (y1-y0) * normale_facettes(fa7,0);
  inverse_norme = 1. / sqrt(a*a + b*b + c*c);
  a *= -inverse_norme;
  b *= -inverse_norme;
  c *= -inverse_norme;
  d = - (a * x0 + b * y0 + c * z0);

  return res;
}


//cette fonction calcule la normale unitaire a la facette, a partir des coordonnees de ses sommets.
//elle renvoie la surface de la facette.
double Maillage_FT_Disc::calcul_normale_3D(long num_facette, double norme[3]) const
{
  double l = -1.;
  long s0 = facettes_(num_facette,0);
  long s1 = facettes_(num_facette,1);
  long s2 = facettes_(num_facette,2);
  double x0 = sommets_(s0,0);
  double y0 = sommets_(s0,1);
  double z0 = sommets_(s0,2);
  double dx1 = sommets_(s1,0) - x0;
  double dy1 = sommets_(s1,1) - y0;
  double dz1 = sommets_(s1,2) - z0;
  double dx2 = sommets_(s2,0) - x0;
  double dy2 = sommets_(s2,1) - y0;
  double dz2 = sommets_(s2,2) - z0;

  //calcul normale : pdt vectoriel
  norme[0] = dy1 * dz2 - dy2 * dz1;
  norme[1] = dz1 * dx2 - dz2 * dx1;
  norme[2] = dx1 * dy2 - dx2 * dy1;
  l = sqrt(norme[0] * norme[0] + norme[1] * norme[1] + norme[2] * norme[2]);
#ifdef _AFFDEBUG
  {
    Process::Journal()<<"Maillage_FT_Disc::calcul_normale_3D fa7= "<<num_facette<<finl;
    Process::Journal()<<"   som 0 coords= "<<x0<<" "<<y0<<" "<<z0<<finl;
    Process::Journal()<<"   som 1 coords= "<<sommets_(s1,0)<<" "<<sommets_(s1,1)<<" "<<sommets_(s1,2)<<finl;
    Process::Journal()<<"   som 2 coords= "<<sommets_(s2,0)<<" "<<sommets_(s2,1)<<" "<<sommets_(s2,2)<<finl;
    Process::Journal()<<"   som 0 coords= "<<x0<<" "<<y0<<" "<<z0<<finl;
    Process::Journal()<<" ->normale= "<<norme[0]<<" "<<norme[1]<<" "<<norme[2]<<"  norme= "<<l<<finl;
  }
#endif
  if (l != 0.)
    {
      double inv_l = 1. / l;
      norme[0] *= inv_l;
      norme[1] *= inv_l;
      norme[2] *= inv_l;
    }

  return l*0.5;
}

/*! @brief Calcule la grandeur demandee, stocke le resultat dans un tableau interne a la classe et renvoie le resultat.
 *
 * Si le maillage
 *   n'a pas change depuis le dernier calcul (mesh_tag identique)
 *   alors on ne recalcule pas la valeur.
 *   Attention, cette methode doit etre appelee simultanement sur
 *   tous les processeurs.
 *
 */
const ArrOfDouble& Maillage_FT_Disc::get_update_surface_facettes() const
{
  Maillage_FT_Disc_Data_Cache& data_cache = mesh_data_cache();
  const long tag = get_mesh_tag();
  if (data_cache.tag_surface_ != tag)
    {
      data_cache.tag_surface_ = tag;
      data_cache.tag_normale_ = tag;
      calcul_surface_normale(data_cache.surface_facettes_,
                             data_cache.normale_facettes_);
    }
  return data_cache.surface_facettes_;
}
// debut EB
void Maillage_FT_Disc::calcul_cg_fa7()
{
  const long nb_fa7=nb_facettes();
  cg_fa7_.resize(nb_fa7,dimension);
  for (long fa7=0; fa7<nb_fa7; fa7++)
    {
      long s0 = facettes_(fa7,0);
      long s1 = facettes_(fa7,1);
      long s2=-1;
      double x2=1e15,y2=1e15;
      double z0=1e15,z1=1e15,z2=1e15;
      if (dimension==3) s2 = facettes_(fa7,2);
      double x0 = sommets_(s0,0);
      double y0 = sommets_(s0,1);
      if (dimension==3) z0 = sommets_(s0,2);
      double x1 = sommets_(s1,0);
      double y1 = sommets_(s1,1);
      if (dimension==3) z1 = sommets_(s1,2);
      if (dimension==3) x2 = sommets_(s2,0);
      if (dimension==3) y2 = sommets_(s2,1);
      if (dimension==3) z2 = sommets_(s2,2);
      cg_fa7_(fa7,0)=(x0+x1+x2)/dimension;
      cg_fa7_(fa7,1)=(y0+y1+y2)/dimension;
      if (dimension==3) cg_fa7_(fa7,2)=(z0+z1+z2)/dimension;
    }
}
// fin EB
/*! @brief Calcule la grandeur demandee, stocke le resultat dans un tableau interne a la classe et renvoie le resultat.
 *
 * Si le maillage
 *   n'a pas change depuis le dernier calcul (mesh_tag identique)
 *   alors on ne recalcule pas la valeur.
 *   Attention, cette methode doit etre appelee simultanement sur
 *   tous les processeurs.
 *
 */
const DoubleTab& Maillage_FT_Disc::get_update_normale_facettes() const
{
  Maillage_FT_Disc_Data_Cache& data_cache = mesh_data_cache();
  const long tag = get_mesh_tag();
  if (data_cache.tag_normale_ != tag)
    {
      data_cache.tag_surface_ = tag;
      data_cache.tag_normale_ = tag;
      calcul_surface_normale(data_cache.surface_facettes_,
                             data_cache.normale_facettes_);
    }
  return data_cache.normale_facettes_;
}

double Maillage_FT_Disc::compute_surface_and_normale_element(const long elem, const bool normalize, double surface, double normale[3]) const
{
  const ArrOfDouble& surface_facettes = get_update_surface_facettes();
  const DoubleTab& normale_facettes = get_update_normale_facettes();
  const Intersections_Elem_Facettes& intersections = intersections_elem_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();
  normale[0]=0.;
  normale[1]=0.;
  normale[2]=0.;// = {0.,0.,0.}; // Normale sortante de I=0 vers I=1
  surface = 0.;
  {
    long index = index_elem[elem];
    if (index<0)
      return surface;
    // Boucle sur les faces qui traversent l'element:
    while (index >= 0)
      {
        const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
        const double fraction_surf = data.fraction_surface_intersection_ * surface_facettes[data.numero_facette_];
        for (long dir= 0; dir<Objet_U::dimension; dir++)
          normale[dir] += fraction_surf * normale_facettes(data.numero_facette_, dir);
        surface += fraction_surf;
        index = data.index_facette_suivante_;
      }
  }

  if (normalize)
    {
      const double norm = sqrt(normale[0]*normale[0]+normale[1]*normale[1]+normale[2]*normale[2]);
      if (norm> Objet_U::precision_geom * Objet_U::precision_geom) // c'est a peu pres une surface pour le moment donc norm tres petit.
        {
          for (long dir= 0; dir<Objet_U::dimension; dir++)
            normale[dir] /= norm;
        }
    }
  return surface;
}

double Maillage_FT_Disc::compute_normale_element(const long elem, const bool normalize, ArrOfDouble& normale) const
{
  const ArrOfDouble& surface_facettes = get_update_surface_facettes();
  const DoubleTab& normale_facettes = get_update_normale_facettes();
  const Intersections_Elem_Facettes& intersections = intersections_elem_facettes();
  const ArrOfInt& index_elem = intersections.index_elem();
  normale=0.;// Normale sortante de I=0 vers I=1
  {
    long index = index_elem[elem];
    if (index<0)
      return 0.;
    // Boucle sur les faces qui traversent l'element:
    while (index >= 0)
      {
        const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
        const double fraction_surf = data.fraction_surface_intersection_ * surface_facettes[data.numero_facette_];
        for (long dir= 0; dir<Objet_U::dimension; dir++)
          normale[dir] += fraction_surf * normale_facettes(data.numero_facette_, dir);
        //surface_totale += fraction_surf;
        index = data.index_facette_suivante_;
      }
    //Cerr << "Surface dans l'elem : " <<  surface_totale << finl;
  }

  const double norm = norme_array(normale);
  if (normalize)
    {
      if (norm> Objet_U::precision_geom * Objet_U::precision_geom) // c'est a peu pres une surface pour le moment donc norm tres petit.
        {
          for (long dir= 0; dir<Objet_U::dimension; dir++)
            normale[dir] /= norm;
        }
    }
  return norm;
}

const ArrOfDouble& Maillage_FT_Disc::get_surface_facettes() const
{
  const Maillage_FT_Disc_Data_Cache& data_cache = mesh_data_cache();
  return data_cache.surface_facettes_;
}

const DoubleTab& Maillage_FT_Disc::get_normale_facettes() const
{
  const Maillage_FT_Disc_Data_Cache& data_cache = mesh_data_cache();
  return data_cache.normale_facettes_;
}

/*! @brief Calcule la grandeur demandee, stocke le resultat dans un tableau interne a la classe et renvoie le resultat.
 *
 * Si le maillage
 *   n'a pas change depuis le dernier calcul (mesh_tag identique)
 *   alors on ne recalcule pas la valeur.
 *   Attention, cette methode doit etre appelee simultanement sur
 *   tous les processeurs.
 *
 */
const ArrOfDouble& Maillage_FT_Disc::get_update_courbure_sommets() const
{
  Maillage_FT_Disc_Data_Cache& data_cache = mesh_data_cache();
  const long tag = get_mesh_tag();
  if (data_cache.tag_courbure_ != tag)
    {
      data_cache.tag_courbure_ = tag;
      calcul_courbure_sommets(data_cache.courbure_sommets_, 1 /*1st call*/);
#if (defined(PATCH_HYSTERESIS_V2) || defined(PATCH_HYSTERESIS_V3) )
      // pour evaluer correctement l'angle de la ligne de contact...
      long ncalls = calcul_courbure_iterations_; // 2;
      for (long i=2; i<=ncalls; i++)
        calcul_courbure_sommets(data_cache.courbure_sommets_, i /*ith call*/);
#endif
    }
  return data_cache.courbure_sommets_;
}

/*! @brief Calcul de la surface et de la normale aux facettes du maillage.
 *
 * Stocke le resultat dans les tableaux en parametres.
 *
 */
void Maillage_FT_Disc::calcul_surface_normale(ArrOfDouble& surface, DoubleTab& normale) const
{
  assert(statut_ >= MINIMAL);
  long i;
  long nfacettes_nulles = 0;
  const long nbfacettes = facettes_.dimension(0);
  surface.resize_array(nbfacettes);
  normale.resize(nbfacettes, dimension);
  if (dimension == 2)
    {
      for (i = 0; i < nbfacettes; i++)
        {
          long s0 = facettes_(i,0);
          long s1 = facettes_(i,1);
          double dx = sommets_(s1,0) - sommets_(s0,0);
          double dy = sommets_(s1,1) - sommets_(s0,1);
          double l = sqrt(dx * dx + dy * dy);
          double inv_l;
          if (l == 0.)
            inv_l = 1.;
          else
            inv_l = 1. / l;

          if (bidim_axi)
            l *= angle_bidim_axi()*(sommets_(s1,0) + sommets_(s0,0))*0.5;

          surface[i] = l;
          if (l == 0.)
            {
              l = 1.;
              dy = -1.;
              nfacettes_nulles++;
            }
          // GB. 2020/17/06. In order to create a unit "normale",
          // it is necessary to have inv_l as 1./sqrt(dx * dx + dy * dy)
          // and not as : 1. / l after the correction!
          // That is why 'inv_l' is calculated before 'if (bidim_axi)'
          // double inv_l = 1. / l;
          normale(i, 0) = - dy * inv_l;
          normale(i, 1) = dx * inv_l;
        }
    }
  else
    {
      double norme[3], s;
      for (i = 0; i < nbfacettes; i++)
        {
          s = calcul_normale_3D(i,norme);
          surface[i] = s;
          if (s == 0.)
            {
              s = 1.;
              nfacettes_nulles++;
#ifdef _AFFDEBUG
              Objet_U::Journal() <<"Fa7 de surf nulle ";
              printFa7(i,0,Objet_U::Journal());
#endif
            }
          normale(i, 0) = norme[0];
          normale(i, 1) = norme[1];
          normale(i, 2) = norme[2];
        }
    }
  if (nfacettes_nulles > 0)
    Objet_U::Journal() << "Facettes de surface nulle : " << nfacettes_nulles << finl;
}

void Maillage_FT_Disc::calculer_voisins()
{
  assert(statut_ >= MINIMAL);
  Cerr<<"ATTENTION : Maillage_FT_Disc::calculer_voisins non implementee !!!"<<finl;
  //assert(0);
}

/*! @brief renvoie le tableau des sommets (reels et virtuels) dimension(0) = nombre de sommets,
 *
 *  dimension(1) = dimension du probleme (2 ou 3)
 *  contenu = coordonnees des sommets
 *
 */
const DoubleTab& Maillage_FT_Disc::sommets() const
{
  assert(statut_ >= MINIMAL);
  return sommets_;
}
// EB
const DoubleTab& Maillage_FT_Disc::cg_fa7() const
{
  assert(statut_ >= MINIMAL);
  return cg_fa7_;
}
/*! @brief renvoie le nombre de sommets (reels et virtuels) (egal a sommets().
 *
 * dimension(0))
 *
 */
long Maillage_FT_Disc::nb_sommets() const
{
  if (statut_< MINIMAL)
    return 0;

  return sommets_.dimension(0);
}

/*! @brief renvoie le tableau des facettes (reelles et virtuelles) dimension(0) = nombre de facettes,
 *
 *  dimension(1) = nombre de sommets par facette (2 en 2D, 3 en 3D)
 *  contenu = numero des sommets de chaque facette dans le tableau des sommets
 *
 */
const IntTab& Maillage_FT_Disc::facettes() const
{
  assert(statut_ >= MINIMAL);
  return facettes_;
}

/*! @brief renvoie le nombre de facettes (reelles et virtuelles) (egal a facettes().
 *
 * dimension(0))
 *
 */
long Maillage_FT_Disc::nb_facettes() const
{
  if (statut_< MINIMAL)
    return 0;

  return facettes_.dimension(0);
}

/*! @brief Cette methode teste si les facettes sont voisines : Des facettes sont voisines si :
 *
 *     -elles ont 1 sommet commun en 2D
 *     -elles ont 2 sommets communs en 3D
 *   La methode renvoie -1 si les facettes ne sont pas voisines.
 *   Elle renvoie l'indice (local) de l'arete (dans fa70) par laquelle les facettes sont voisines
 *  A OPTIMISER eventuellement
 *
 */
long Maillage_FT_Disc::facettes_voisines(long fa70, long fa71,
                                         long& iarete0, long& iarete1) const
{
  long isom0,isom1, som0,som1;
  long res = 0;
  iarete1 = iarete0 = -1;
  const long nb_som_par_fa7 = facettes_.dimension(1);
  long compteur = 0, trouve;

  //on marque aussi les sommets communs
  long somsCommuns0[3];
  long somsCommuns1[3];
  somsCommuns0[0] = somsCommuns0[1] = somsCommuns0[2] = 0;
  somsCommuns1[0] = somsCommuns1[1] = somsCommuns1[2] = 0;

  for (isom0=0 ; isom0<nb_som_par_fa7 ; isom0++)
    {
      som0 = facettes_(fa70,isom0);
      trouve = 0;
      for (isom1=0 ; isom1<nb_som_par_fa7 && trouve==0 ; isom1++)
        {
          som1 = facettes_(fa71,isom1);
          if (som1==som0)
            {
              compteur++;
              trouve = 1;
              somsCommuns0[isom0]++;
              somsCommuns1[isom1]++;
            }
        }
    }

  //s'assure que le nb de sommets communs est correct
  if (compteur==(nb_som_par_fa7-1))
    {
      //les facettes sont donc bien voisines
      //et recupere les indices des aretes de voisinage dans les facettes
      res = 1;
      if (dimension==2)
        {
          iarete0 = (somsCommuns0[0]==1) ? (0) : (1);
          iarete1 = (somsCommuns1[0]==1) ? (0) : (1);
          assert(facettes_(fa70,iarete0) == facettes_(fa71,iarete1));
        }
      else
        {
          //on determine les aretes communes, en fonction des sommets communs
          if (somsCommuns0[0]==1)
            {
              if (somsCommuns0[1]==1)
                {
                  iarete0 = 0;
                }
              else
                {
                  iarete0 = 2;
                }
            }
          else
            {
              iarete0 = 1;
            }
          if (somsCommuns1[0]==1)
            {
              if (somsCommuns1[1]==1)
                {
                  iarete1 = 0;
                }
              else
                {
                  iarete1 = 2;
                }
            }
          else
            {
              iarete1 = 1;
            }
#ifdef _AFFDEBUG
          {
            Process::Journal()<<"test facettes_voisines compteur="<<compteur<<finl;
            Process::Journal()<<"  ";
            printFa7(fa70,0,Process::Journal());
            Process::Journal()<<"  ";
            printFa7(fa71,0,Process::Journal());
            Process::Journal()<<"   somsCommuns0="<<somsCommuns0[0]<<" "<<somsCommuns0[1]<<" "<<somsCommuns0[2]<<finl;
            Process::Journal()<<"   somsCommuns1="<<somsCommuns1[0]<<" "<<somsCommuns1[1]<<" "<<somsCommuns1[2]<<finl;
            Process::Journal()<<"   aretes  voisines iarete0="<<iarete0<<" iarete1="<<iarete1<<finl;
            Process::Journal()<<"     test = "<<(facettes_(fa70,iarete0)==facettes_(fa71,iarete1)?1:0)<<" res="<<res<<finl;
          }
#endif
          //attention : il faut aussi que le sens de parcours de l'arete soit opposse
          if (facettes_(fa70,iarete0)==facettes_(fa71,iarete1))
            {
              //les 2 facettes parcourenet l'arete dans le meme sens : elles ne sont pas voisines
              //(ou voisines mais pas compatibles au sens du maillage)
              res = 0;
            }
          else
            {
              if (! (facettes_(fa70,iarete0) == facettes_(fa71,(iarete1+1)%nb_som_par_fa7)) ||
                  ! (facettes_(fa71,iarete1) == facettes_(fa70,(iarete0+1)%nb_som_par_fa7)) )
                {
                  Process::Journal()<<"PB dans facettes_voisines : fa7 voisines iarete0="<<iarete0<<" iarete1="<<iarete1<<finl;
                  Process::Journal()<<"  ";
                  printFa7(fa70,0,Process::Journal());
                  Process::Journal()<<"  ";
                  printFa7(fa71,0,Process::Journal());
                  Process::Journal()<<"   somsCommuns0="<<somsCommuns0[0]<<" "<<somsCommuns0[1]<<" "<<somsCommuns0[2]<<finl;
                  Process::Journal()<<"   somsCommuns1="<<somsCommuns1[0]<<" "<<somsCommuns1[1]<<" "<<somsCommuns1[2]<<finl;
                  Process::Journal()<<"   aretes  voisines iarete0="<<iarete0<<" iarete1="<<iarete1<<finl;
                  Process::Journal()<<"     test = "<< (long) (facettes_(fa70,iarete0)==facettes_(fa71,iarete1)?1:0)<<" res="<<res<<finl;
                }
              assert( facettes_(fa70,iarete0) == facettes_(fa71,(iarete1+1)%nb_som_par_fa7) );
              assert( facettes_(fa71,iarete1) == facettes_(fa70,(iarete0+1)%nb_som_par_fa7) );
            }
        }
    }

  return res;
}

//Description :
// Cette calcule les voisinages des facettes, en parcourant les intersections avec les elements
// Le maillage doit etre en etat PARCOURU.
//A OPTIMISER
long Maillage_FT_Disc::calculer_voisinage_facettes(IntTab& fa7Voisines,
                                                   const Intersections_Elem_Facettes* ief) const
{
  Process::Journal()<<"Maillage_FT_Disc::calculer_voisinage_facettes"<<finl;
  long res = 1;
  if (ief==NULL)
    {
      ief = & intersections_elem_facettes_;
      assert(statut_ >= PARCOURU);
    }
  fa7Voisines.resize(nb_facettes(),facettes_.dimension(1));
  fa7Voisines = -1;

  ArrOfIntFT fa7s_elem;
  long i, nb_fa7, ifa70,ifa71,fa70,fa71, iarete0,iarete1;
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine& ladomaine = domaine_dis.domaine();
  const long nb_elem = ladomaine.nb_elem(); // Nombre d'elements reels
  for (i=0 ; i<nb_elem ; i++)
    {
      ief->get_liste_facettes_traversantes(i,fa7s_elem);
      nb_fa7 = fa7s_elem.size_array();
      //on va ensuite tester tous les couples de fa7 posssibles
      for (ifa70=0 ; ifa70<nb_fa7 ; ifa70++)
        {
          fa70 = fa7s_elem[ifa70];
          for (ifa71=ifa70+1 ; ifa71<nb_fa7 ; ifa71++)
            {
              fa71 = fa7s_elem[ifa71];
              //on teste d'abord si les facettes sont voisines
              if (facettes_voisines(fa70,fa71,iarete0,iarete1))
                {
                  //les facettes sont bien voisines : on memorise ca dans le tableau de voisinage
                  fa7Voisines(fa70,iarete0) = fa71;
                  fa7Voisines(fa71,iarete1) = fa70;
                }
            }
        }
    }//fin du balayage des elements

  return res;
}

/*! @brief Pour postraitement uniquement : la signification precise des drapeaux est de la cuisine interne a la classe.
 *
 * ..
 *
 */
const ArrOfInt& Maillage_FT_Disc::drapeaux_sommets() const
{
  assert(statut_ >= MINIMAL);
  return drapeaux_sommets_;
}

/*! @brief renvoie le descripteur des sommets (espace_distant/virtuel)
 *
 */
const Desc_Structure_FT& Maillage_FT_Disc::desc_sommets() const
{
  assert(statut_ >= MINIMAL);
  return desc_sommets_;
}
/*! @brief renvoie le descripteur des facettes (espace_distant/virtuel)
 *
 */
const Desc_Structure_FT& Maillage_FT_Disc::desc_facettes() const
{
  assert(statut_ >= MINIMAL);
  return desc_facettes_;
}

/*! @brief pour postraitement, renvoie le numero du PE proprietaire des sommets
 *
 */
const ArrOfInt& Maillage_FT_Disc::sommet_PE_owner() const
{
  assert(statut_ >= MINIMAL);
  return sommet_PE_owner_;
}
/*! @brief pour postraitement, renvoie le numero des sommets sur le PE proprietaire des sommets
 *
 */
const ArrOfInt& Maillage_FT_Disc::sommet_num_owner() const
{
  assert(statut_ >= MINIMAL);
  return sommet_num_owner_;
}

/*! @brief pour postraitement, remplit le tableau en parametre avec le numero du PE proprietaire de chaque facette.
 *
 */
void Maillage_FT_Disc::facette_PE_owner(ArrOfInt& facette_pe) const
{
  assert(statut_ >= MINIMAL);
  facette_pe.resize_array(facettes_.dimension(0));
  desc_facettes_.remplir_element_pe(facette_pe);
}

/*! @brief pour postraitement, renvoie sommet_elem_
 *
 */
const ArrOfInt& Maillage_FT_Disc::sommet_elem() const
{
  assert(statut_ >= MINIMAL);
  return sommet_elem_;
}
// debut EB
// Description: pour postraitement, renvoie sommet_face_
const IntTab& Maillage_FT_Disc::sommet_face() const
{
  assert(statut_ >= MINIMAL);
  return sommet_face_;
}
// fin EB
/*! @brief pour postraitement, renvoie sommet_face_bord_
 *
 */
const ArrOfInt& Maillage_FT_Disc::sommet_face_bord() const
{
  assert(statut_ >= MINIMAL);
  return sommet_face_bord_;
}

/*! @brief return temps_physique_ (temps_physique_ ne sert a rien en interne.
 *
 * ..)
 *
 */
double Maillage_FT_Disc::temps() const
{
  return temps_physique_;
}

/*! @brief return temps_physique_ = t
 *
 */
double Maillage_FT_Disc::changer_temps(double t)
{
  temps_physique_ = t;
  return t;
}

/*! @brief return mesh_state_tag_
 *
 */
long Maillage_FT_Disc::get_mesh_tag() const
{
  return mesh_state_tag_;
}

/*! @brief return som_init_util_
 *
 */
const ArrOfInt& Maillage_FT_Disc::som_init_util() const
{
  return som_init_util_;
}

long Maillage_FT_Disc::sauvegarder(Sortie& os) const
{
  long special, afaire;
  const long format_xyz = EcritureLectureSpecial::is_ecriture_special(special, afaire);
  if (format_xyz)
    {
      const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());
      Sauvegarde_Reprise_Maillage_FT::ecrire_xyz(*this, domaine_vf, os);
      return 0;
    }
  else
    {
      long bytes = 0;
      os << que_suis_je() << finl;
      os << mesh_state_tag_ << finl;
      os << temps_physique_ << finl;
      // On augmente la precision d'ecriture du maillage Front Tracking
      // au moment de l'ecriture et on restaure l'ancienne apres
      // Probleme rencontre sur la FA819
      long old_precision = (long)os.get_ostream().precision(14);
      os << sommets_;
      bytes += 8 * sommets_.size_array();
      os.precision(old_precision);
      os << facettes_;
      bytes += 4 * facettes_.size_array();
      os << voisins_;
      bytes += 4 * voisins_.size_array();
      os << sommet_elem_;
      bytes += 4 * sommet_elem_.size_array();
      os << sommet_face_; // EB
      bytes += 4 * sommet_face_.size_array(); // EB
      os << sommet_face_bord_;
      bytes += 4 * sommet_face_bord_.size_array();
      os << sommet_PE_owner_;
      bytes += 4 * sommet_PE_owner_.size_array();
      os << sommet_num_owner_;
      bytes += 4 * sommet_num_owner_.size_array();
      os << facette_num_owner_;
      bytes += 4 * facette_num_owner_.size_array();
      os << desc_sommets_;
      os << desc_facettes_;
      os << drapeaux_sommets_;
      bytes += 4 * drapeaux_sommets_.size_array();
      return bytes;
    }
}

long Maillage_FT_Disc::reprendre(Entree& is)
{
  reset();
  const long format_xyz =
    EcritureLectureSpecial::is_lecture_special();
  if (format_xyz)
    {
      if (refdomaine_dis_.non_nul())
        {
          const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
          const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());
          Sauvegarde_Reprise_Maillage_FT::lire_xyz(*this, &domaine_vf, &is, 0);
        }
      else
        {
          // On est dans "avancer", equation non associee
          Sauvegarde_Reprise_Maillage_FT::lire_xyz(*this, 0, &is, 0);
        }
    }
  else
    {
      Nom motlu;
      is >> motlu;
      if (motlu != que_suis_je())
        {
          Cerr << "Erreur dans Maillage_FT_Disc::reprendre\n";
          Cerr << " On attendait le motcle " << que_suis_je();
          Cerr << "\n On a trouve " << motlu << finl;
          Process::exit();
        }
      is >> mesh_state_tag_;
      is >> temps_physique_;
      is >> sommets_;
      is >> facettes_;
      is >> voisins_;
      is >> sommet_elem_;
      is >> sommet_face_; // EB
      is >> sommet_face_bord_;
      is >> sommet_PE_owner_;
      is >> sommet_num_owner_;
      is >> facette_num_owner_;
      is >> desc_sommets_;
      is >> desc_facettes_;
      is >> drapeaux_sommets_;
      maillage_modifie(MINIMAL);
      if (refdomaine_dis_.non_nul())
        {
          // L'equation de transport a ete associee : on reprend le vrai
          // maillage (c'est pas une lecture bidon pour passer au bloc suivant
          // dans avancer_fichier...)
          desc_sommets_.calcul_schema_comm(sommets_.dimension(0));
          desc_facettes_.calcul_schema_comm(facettes_.dimension(0));
          check_mesh();
        }
    }
  return 1;
}

/*! @brief fonction qui cree un nouveau sommet par copie d'une existant utilise dans Remailleur_Collision_FT_Collision_Seq
 *
 */
long Maillage_FT_Disc::copier_sommet(long som)
{
  maillage_modifie(MINIMAL);
  if (me()!=sommet_PE_owner_[som])
    {
      //si je ne suis pas le proprietaire du sommet de reference : on abandonne
      return -1;
    }

  long NVsom = nb_sommets();
  long nb_sommets_tot = NVsom+1;
  //on redimensionne les tableaux
  sommets_.resize(nb_sommets_tot,dimension);
  sommet_elem_.resize_array(nb_sommets_tot);
  sommet_face_.resize(nb_sommets_tot,dimension); // EB /!\ On a 3 indicatrices aux faces. Un sommet appartient a 3 faces
  sommet_face_bord_.resize_array(nb_sommets_tot);
  sommet_PE_owner_.resize_array(nb_sommets_tot);
  sommet_num_owner_.resize_array(nb_sommets_tot);
  drapeaux_sommets_.resize_array(nb_sommets_tot);

  //copie des donnes dui sommet som dans NVsom
  long k;
  for (k=0 ; k<dimension ; k++)
    {
      sommets_(NVsom,k)      = sommets_(som,k);
    }
  sommet_elem_[NVsom]      = sommet_elem_[som];
  for (long dim=0; dim<dimension; dim++) sommet_face_(NVsom,dim)      = sommet_face_(som,dim); // EB
  sommet_face_bord_[NVsom] = sommet_face_bord_[som];
  sommet_PE_owner_[NVsom]  = sommet_PE_owner_[som];
  sommet_num_owner_[NVsom] = NVsom;
  drapeaux_sommets_[NVsom] = drapeaux_sommets_[som];

  return NVsom;
}

//idem ci-dessus, mais rend le sommet interne (sommet_face_bord_=-1)
long Maillage_FT_Disc::copier_sommet_interne(long som)
{
  long NVsom = copier_sommet(som);

  if (NVsom>=0)
    {
      //rend le sommet interne
      sommet_face_bord_[som] = -1;
    }

  return NVsom;
}


/*! @brief Envoi des sommets dont le numero est donne dans liste_sommets au processeur dont le numero est donne dans liste_pe (on cree un
 *
 *  sommet virtuel sur ce processeur).  Parmi les sommets de la liste_sommets,
 *  seuls ceux qui ne sont pas encore dans l'espace virtuel sont crees.
 *
 * Precondition: Les membres suivants doivent etre valides:
 * * sommets_,
 * * desc_sommets_,
 * * drapeaux_sommets_,
 * * sommet_elem_,
 * * sommet_face_, // EB
 * * sommet_PE_owner_,
 * * sommet_num_owner_
 *
 * Postcondition: Les memes tableaux sont valides a la sortie de la fonction.
 *
 * @param (liste_sommets) une liste de numeros de sommets REELS qui peut contenir des doublons.
 * @param (liste_pe) une liste de meme taille que liste_sommets contenant le numero du processeur sur lequel il faut creer un noeud virtuel. Le processeur ne doit pas etre moi.
 * @param (comm) un schema de comm valide, dans lequel send_pe_list contient au moins les processeurs cites dans liste_pe.
 */
void Maillage_FT_Disc::creer_sommets_virtuels(const ArrOfInt& liste_sommets,
                                              const ArrOfInt& liste_pe,
                                              const Schema_Comm_FT& comm)
{
  if (Comm_Group::check_enabled()) check_sommets();

  const long dimension3 = (dimension == 3);
  long i;
  const long n = liste_sommets.size_array();
  Descripteur_FT& espace_distant = desc_sommets_.espace_distant();

  comm.begin_comm();

  for (i = 0; i < n; i++)
    {
      long pe_destination = liste_pe[i];
      assert(pe_destination != me());
      long numero_sommet = liste_sommets[i];
      // Si le sommet est deja connu par le destinataire, on ne l'envoie pas.
      if (espace_distant.contient_element(pe_destination, numero_sommet))
        continue;

      // Ajout du sommet a l'espace distant
      espace_distant.ajoute_element(pe_destination, numero_sommet);

      // Envoi a l'autre processeur
      Sortie& send_buffer = comm.send_buffer(pe_destination);
      // On verifie que le sommet est bien reel
      assert(! sommet_virtuel(numero_sommet));
      long drapeau = drapeaux_sommets_[numero_sommet];
      long face_bord = sommet_face_bord_[numero_sommet];
      double x = sommets_(numero_sommet,0);
      double y = sommets_(numero_sommet,1);
      double z = dimension3 ? sommets_(numero_sommet, 2) : 0.;
      send_buffer << numero_sommet << drapeau << face_bord << x << y << z;
    }

  comm.echange_taille_et_messages();

  Descripteur_FT& espace_virtuel = desc_sommets_.espace_virtuel();
  // Liste des processeurs de qui on a recu des messages
  const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
  const long recv_pe_size = recv_pe_list.size_array();
  long indice_pe;
  for (indice_pe = 0; indice_pe < recv_pe_size; indice_pe++)
    {
      const long pe_source = recv_pe_list[indice_pe];
      Entree& recv_buffer = comm.recv_buffer(pe_source);

      do
        {
          // Numero du sommet sur le PE proprietaire :
          long numero_sur_pe_source, drapeau, face_bord;
          double x, y, z;
          recv_buffer >> numero_sur_pe_source >> drapeau >> face_bord >> x >> y >> z;
          if (recv_buffer.eof())
            break;

          // Creation du sommet
          long nsom = sommets_.dimension(0);
          if (dimension3)
            sommets_.append_line(x, y, z);
          else
            sommets_.append_line(x, y);
          espace_virtuel.ajoute_element(pe_source, nsom);
          sommet_elem_.append_array(-1); // C'est un sommet virtuel => -1
          if (dimension==3 && sommet_face_.dimension(0)>0 && sommet_face_.dimension(1)>0) sommet_face_.append_line(-1,-1,-1); // EB : C'est un sommet virtuel => -1. /!\ Si faces doubles
          else if (dimension==2 && sommet_face_.dimension(0)>0 && sommet_face_.dimension(1)>0) sommet_face_.append_line(-1,-1);
          sommet_face_bord_.append_array(face_bord);
          sommet_PE_owner_.append_array(pe_source);
          sommet_num_owner_.append_array(numero_sur_pe_source);
          drapeaux_sommets_.append_array(drapeau);
        }
      while (1);
    }
  comm.end_comm();

  // Recalcul des listes de voisins et du schema de comm
  espace_distant.calcul_liste_pe_voisins();
  espace_virtuel.calcul_liste_pe_voisins();
  desc_sommets_.calcul_schema_comm(sommets_.dimension(0));

  if (Comm_Group::check_enabled()) check_sommets();
}

/*! @brief Cree chez moi les sommets virtuels specifies par le couple (pe,num) si le sommet n'existe pas encore.
 *
 *   Le sommet est suppose etre reel sur "pe" et le sommet virtuel est cree sur "me()".
 *   (c'est l'inverse de "creer_sommets_virtuels" qui cree des sommets virtuels chez
 *    chez le pe demande).
 *   Fait appel a creer_sommets_virtuels, donc meme preconditions.
 *
 * @param (request_sommets_pe) une liste de numeros de processeurs DIFFERENTS de me()
 * @param (request_sommets_num) une liste de numeros de sommets. Chaque couple pe[i],num[i] designe un sommet reel d'indice "num" sur le processeur "pe". Pour chaque sommet, on cree un sommet virtuel sur me() s'il n'existe pas deja.
 */
void Maillage_FT_Disc::creer_sommets_virtuels_numowner(const ArrOfInt& request_sommets_pe,
                                                       const ArrOfInt& request_sommets_num)
{
  const long nbproc = Process::nproc();
  // Liste des processeurs a qui on envoyer des requetes
  //   et de qui on va recevoir des noeuds:
  ArrOfIntFT send_pe_list;
  // Liste des processeurs de qui on va recevoir des requetes
  //   et a qui on va envoyer des noeuds:
  ArrOfIntFT recv_pe_list;
  // Pour chaque processeur, est-il dans send_pe_list ?
  ArrOfIntFT pe_markers(nbproc);
  pe_markers = 0;
  {
    long i;
    const long n = request_sommets_pe.size_array();
    assert(n == request_sommets_num.size_array());
    for (i = 0; i < n; i++)
      {
        const long pe = request_sommets_pe[i];
        assert(pe != Process::me());
        if (pe_markers[pe] == 0)
          {
            pe_markers[pe] = 1;
            send_pe_list.append_array(pe);
          }
      }
  }
  // Calcule la send_pe_list en fonction de la recv_pe_list.
  // C'est long (communication all_to_all) mais c'est le plus simple
  // car les processeurs de la recv_pe_list ne savent pas a priori
  // de quels processeurs ils vont recevoir des messages.
  {
    long sz = send_pe_list.size_array();
    ArrOfInt s(sz);
    ArrOfInt r;
    long i;
    for (i = 0; i < sz; i++) s[i] = send_pe_list[i];
    reverse_send_recv_pe_list(s, r);
    sz = r.size_array();
    recv_pe_list.resize_array(sz);
    for (i = 0; i < sz; i++) recv_pe_list[i] = r[i];
  }

  // Envoi des requetes de creation de noeuds aux proprietaires des noeuds.
  {
    Schema_Comm_FT schema;
    schema.set_send_recv_pe_list(send_pe_list, recv_pe_list);
    schema.begin_comm();
    {
      long i;
      const long n = request_sommets_pe.size_array();
      for (i = 0; i < n; i++)
        {
          const long pe = request_sommets_pe[i];
          const long som = request_sommets_num[i];
          schema.send_buffer(pe) << som;
        }
    }
    schema.echange_taille_et_messages();

    // Liste des requetes de creation de noeuds recues:
    ArrOfIntFT liste_sommets_pe;  // Sur quel proc faut-il creer un noeud virt ?
    ArrOfIntFT liste_sommets_num; // Numero local du noeud a envoyer
    {
      long i_pe;
      const long n_pe = recv_pe_list.size_array();
      for (i_pe = 0; i_pe < n_pe; i_pe++)
        {
          const long pe = recv_pe_list[i_pe];
          Entree& buffer = schema.recv_buffer(pe);
          long num_sommet;
          while (1)
            {
              buffer >> num_sommet;
              if (buffer.eof())
                break;
              assert(num_sommet >= 0 && num_sommet < nb_sommets());
              liste_sommets_pe.append_array(pe);
              liste_sommets_num.append_array(num_sommet);
            }
        }
    }
    schema.end_comm();
    // Envoi des sommets virtuels requis
    // On intervertit les listes "envoi" et "reception"
    schema.set_send_recv_pe_list(recv_pe_list, send_pe_list);
    creer_sommets_virtuels(liste_sommets_num,
                           liste_sommets_pe,
                           schema);
  }
}

/*! @brief Realise l'ensemble des echanges de noeuds specifies dans liste_sommets, liste_elem_virtuel_arrivee et deplacement.
 *
 *  Met a jour les sommets, facettes, espaces distants
 *  et virtuels. Le maillage retourne a l'etat minimal. Le buffer d'echange est
 *  reinitialise.
 *  On stocke dans liste_nouveaux_sommets la liste des numeros des sommets reels
 *  fraichement arrives sur le domaine, ainsi que le deplacement restant pour
 *  chacun de ces sommets (tableau a 3 colonnes (meme en dimension 2 !), de taille
 *  le nombre de nouveaux sommets).
 *
 * Precondition:
 * Les membres suivants doivent etre valides a l'entree de la fonction:
 *  (Meme liste que pour creer_sommets_virtuels(...) )
 * Postcondition:
 *  Les memes membres sont valides a la sortie.
 *  Attention: on ne met PAS a jour les facettes, la regle de propriete des
 *            facettes n'est donc pas verifiee au retour de la fonction.
 *            Il faut appeler corriger_proprietaire_facettes.
 *
 */
void Maillage_FT_Disc::echanger_sommets_PE(const ArrOfInt& liste_sommets,
                                           const ArrOfInt& liste_elem_virtuel_arrivee,
                                           const ArrOfInt& liste_face_virtuelle_arrivee,
                                           const DoubleTab& deplacement,
                                           ArrOfInt& liste_nouveaux_sommets,
                                           DoubleTab& deplacement_restant,
                                           long skip_facettes)
{
  //Cerr <<"Maillage_FT_Disc::echanger_sommets_PE" << finl;

  if (Comm_Group::check_enabled()) check_mesh(1 /* error is fatal */,
                                                1 /* do not test facette_owner */,
                                                skip_facettes);
  assert(statut_ >= MINIMAL);

  long i;
  // Nombre de sommets a echanger
  const long nechange = liste_sommets.size_array();

  // Numeros des processeurs a qui on envoie les noeuds
  ArrOfIntFT liste_pe(nechange);
  // Numeros des elements ou arrivent les noeuds (numero local sur ce proc)
  ArrOfIntFT liste_elem_arrivee(nechange);
  // Numeros des faces ou arrivent les noeuds "ligne de contact"
  ArrOfIntFT liste_face_arrivee(nechange);
  // EB des faces xyz qui ou arrivent les noeuds 'numero local sur ce proc)
  IntTabFT liste_faces_xyz(nechange,dimension);

  // Creation des noeuds virtuels sur le processeur d'arrivee s'ils n'existent
  // pas encore.
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine& ladomaine = domaine_dis.domaine();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());
  const IntTab& elem_virt_pe_num = ladomaine.elem_virt_pe_num();
  const IntTab& face_virt_pe_num = domaine_vf.face_virt_pe_num();
  const long nb_elements_reels = ladomaine.nb_elem();
  const long nb_faces_reelles = domaine_vf.nb_faces();

  for (i = 0; i < nechange; i++)
    {
      // Numero de l'element virtuel ou arrive le sommet
      const long numero_element = liste_elem_virtuel_arrivee[i];
      // Numero du pe qui possede l'element eulerien
      // (si l'index est negatif, c'est que l'element d'arrivee n'est pas virtuel !)
      const long index = numero_element - nb_elements_reels;
      const long pe_destination = elem_virt_pe_num(index, 0);
      // Numero de l'element ou arrive le sommet sur l'autre pe
      const long element_d_arrivee = elem_virt_pe_num(index, 1);

      // Numero de la face de bord d'arrivee sur l'autre pe
      // si le sommet est sur une ligne de contact.
      const long numero_face = liste_face_virtuelle_arrivee[i];
      // On verifie qu'on connait la face d'arrivee pour les sommets de bord
      assert( ((sommet_ligne_contact(liste_sommets[i]) != 0) && (numero_face >= 0))
              || ((sommet_ligne_contact(liste_sommets[i]) == 0) && (numero_face < 0)));
      long face_d_arrivee = -1;
      if (numero_face >= 0)
        {
          // Calcul du numero local de la face d'arrivee sur le processeur voisin :
          const long index_face = numero_face - nb_faces_reelles;
          face_d_arrivee = face_virt_pe_num(index_face, 1);
        }

      liste_pe[i]           = pe_destination;
      liste_elem_arrivee[i] = element_d_arrivee;
      liste_face_arrivee[i] = face_d_arrivee;
    }

  // Les echanges de sommets se font entre voisins par des joints du maillage
  // fixe. On prend donc le schema de communication de la zone.
  // Envoi des sommets aux destinataires (pour l'instant, creation de sommets
  // virtuels sur le processeur destination).
  creer_sommets_virtuels(liste_sommets, liste_pe, schema_comm_domaine_);

  // Mise a jour du nouveau PE proprietaire, de l'element contenant le sommet et de la
  // face de bord :
  //  * Remplissage du tableau sommet_PE_owner avec le nouveau proprietaire
  //    et -1 si le sommet ne change pas de main,
  //  * Remplissage de sommet_elem_ avec le numero de l'element d'arrivee,
  //  * Remplissage de sommet_face_bord_ avec le numero de la face,
  //  * Appel a la methode echanger_elements qui transforme
  //    les espaces distants et virtuels en fonction du nouveau proprietaire.
  sommet_PE_owner_ = -1;
  for (i = 0; i < nechange; i++)
    {
      const long sommet = liste_sommets[i];
      const long pe_destination = liste_pe[i];
      const long element_d_arrivee = liste_elem_arrivee[i];
      const long face_d_arrivee = liste_face_arrivee[i];
      // Verifie que le meme sommet ne figure pas deux fois dans la liste :
      assert(sommet_PE_owner_[sommet] < 0);
      assert(pe_destination != Process::me());
      sommet_PE_owner_[sommet] = pe_destination;
      // element_d_arrivee est un numero d'element eulerien reel sur le pe_destination
      sommet_elem_[sommet] = element_d_arrivee;
      sommet_face_bord_[sommet] = face_d_arrivee;
    }
  desc_sommets_.echange_espace_virtuel(sommet_PE_owner_);
  // Envoi du numero d'element contenant les sommets a tous les processeurs
  desc_sommets_.echange_espace_virtuel(sommet_elem_);
  desc_sommets_.echange_espace_virtuel(sommet_face_bord_);
  // Mise a jour des espaces distants et virtuels (c'est la qu'on change
  // reelement de proprietaire)
  desc_sommets_.echanger_elements(sommet_PE_owner_);
  // Recalcul de sommet_PE_owner a partir des descripteurs (on l'avait ecrase
  // avec -1 partout sauf pour les sommets qui changent de procs).
  desc_sommets_.remplir_element_pe(sommet_PE_owner_);
  // On remet a -1 le numero d'element si le sommet n'est pas a moi
  const long nbsommets = sommets_.dimension(0);
  const long moi = me();
  for (i = 0; i < nbsommets; i++)
    {
      if (sommet_PE_owner_[i] != moi)
        {
          sommet_elem_[i] = -1;
          //for (long dim=0; dim<dimension; dim++) sommet_face_(i,dim) = -1; // EB
        }
      else
        assert(sommet_elem_[i] >= 0 && sommet_elem_[i] < nb_elements_reels);
    }
  // Recalcul de sommet_num_owner
  for (i = 0; i < nbsommets; i++)
    sommet_num_owner_[i] = i;
  desc_sommets_.echange_espace_virtuel(sommet_num_owner_);

  // Envoi du numero du sommet et du deplacement restant
  {
    liste_nouveaux_sommets.resize_array(0);
    deplacement_restant.resize(0, Objet_U::dimension);
    // Choix d'un schema de communication:
    // Apres l'echange des sommets, le sommet a envoyer est virtuel sur l'expediteur
    // et reel sur le destinataire : on prend le schema
    // "les procs ayant des sommets virtuels parlent aux procs ayant des sommets reels"
    // (on pourrait reduire cet ensemble).
    const Schema_Comm_FT& comm = desc_sommets_.schema_comm_inverse();

    comm.begin_comm();
    // Remplissage des buffers
    for (i = 0; i < nechange; i++)
      {
        const long pe_destination = liste_pe[i];
        const long sommet = liste_sommets[i];
        const long num_owner = sommet_num_owner_[sommet];
        const double x = deplacement(i, 0);
        const double y = deplacement(i, 1);
        const double z = deplacement(i, 2);
        Sortie& send_buffer = comm.send_buffer(pe_destination);
        send_buffer << num_owner << x << y << z;
      }
    comm.echange_taille_et_messages();
    // Recuperation des donnees
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long recv_pe_size = recv_pe_list.size_array();
    long n_recus = 0;
    const long dim = Objet_U::dimension;
    for (long indice_pe = 0; indice_pe < recv_pe_size; indice_pe++)
      {
        const long pe_source = recv_pe_list[indice_pe];
        Entree& recv_buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long sommet;
            double x, y, z;
            recv_buffer >> sommet >> x >> y >> z;
            if (recv_buffer.eof())
              break;
            liste_nouveaux_sommets.append_array(sommet);
            deplacement_restant.resize(n_recus+1, dim);
            deplacement_restant(n_recus, 0) = x;
            deplacement_restant(n_recus, 1) = y;
            if (dim==3)
              deplacement_restant(n_recus, 2) = z;
            n_recus++;
          }
      }
    comm.end_comm();
  }

  statut_ = MINIMAL;

  if (Comm_Group::check_enabled()) check_mesh(1 /* error is fatal */,
                                                1 /* do not test facette_owner */,
                                                skip_facettes);
}
void Maillage_FT_Disc::update_sommet_face()
{
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF,domaine_dis.valeur());
  const long nb_som=sommets_.dimension(0);
  sommet_face_.resize(nb_som,dimension);

  long mon_elem;

  for (long som=0; som<nb_som; som++)
    {
      for (long dim=0; dim<dimension; dim++) sommet_face_(som,dim)=-1;
      mon_elem=sommet_elem_(som);

      if (mon_elem>=0)
        {
          IntVect faces_elem_eulerien(2*dimension);
          for (long dim=0; dim<2*dimension; dim++) faces_elem_eulerien(dim) = domaine_vf.elem_faces(mon_elem,dim); // on recupere les faces de l'element eulerien
          double pos_sommet, coord_elem;
          for (long dim=0; dim<dimension; dim++)
            {
              pos_sommet=sommets_(som,dim); // on recupere la position du sommet suivant l'axe "dim"
              coord_elem=domaine_vf.xp(mon_elem,dim);
              if (pos_sommet < coord_elem ) sommet_face_(som,dim)=faces_elem_eulerien(dim); //&& check_som_in_face(dim)
              else if (pos_sommet >= coord_elem) sommet_face_(som,dim)=faces_elem_eulerien(dim+dimension); // && check_som_in_face(dim+dimension)
              assert(sommet_face_(som,dim)>=0);
            }
        }

    }
}

void Maillage_FT_Disc::update_sommet_arete()
{
  assert(dimension==3);
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine& domaine = domaine_dis.domaine();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());
  const long nb_som=sommets_.dimension(0);
  sommet_arete_.resize(nb_som,dimension);
  long nb_arete_elem=12;
  long mon_elem;
  IntVect aretes_elem_eulerien(nb_arete_elem);
  const IntTab& elem_aretes = domaine.elem_aretes();
  long arete_bas_gauche,arete_bas_droite,arete_haut_gauche,arete_haut_droite;
  double cgx,cgy,cgz;
  double somx,somy,somz;
  for (long som=0; som<nb_som; som++)
    {
      for (long dim=0; dim<dimension; dim++) sommet_arete_(som,dim)=-1;
      mon_elem=sommet_elem_(som);
      if (mon_elem>=0)
        {
          // cg elem
          cgx=domaine_vf.xp(mon_elem,0);
          cgy=domaine_vf.xp(mon_elem,1);
          cgz=domaine_vf.xp(mon_elem,2);

          // position sommet
          somx=sommets_(som,0);
          somy=sommets_(som,1);
          somz=sommets_(som,2);

          // arete XY
          arete_bas_gauche=elem_aretes(mon_elem,3);
          arete_bas_droite=elem_aretes(mon_elem,6);
          arete_haut_gauche=elem_aretes(mon_elem,9);
          arete_haut_droite=elem_aretes(mon_elem,0);

          if (somx<cgx)
            {
              if (somy<cgy) sommet_arete_(som,2)=arete_bas_gauche;
              else sommet_arete_(som,2)=arete_haut_gauche;
            }
          else
            {
              if (somy<cgy) sommet_arete_(som,2)=arete_bas_droite;
              else sommet_arete_(som,2)=arete_haut_droite;
            }

          // arete XZ
          arete_bas_gauche=elem_aretes(mon_elem,4);
          arete_bas_droite=elem_aretes(mon_elem,7);
          arete_haut_gauche=elem_aretes(mon_elem,10);
          arete_haut_droite=elem_aretes(mon_elem,1);

          if (somx<cgx)
            {
              if (somz<cgz) sommet_arete_(som,1)=arete_bas_gauche;
              else sommet_arete_(som,1)=arete_haut_gauche;
            }
          else
            {
              if (somz<cgz) sommet_arete_(som,1)=arete_bas_droite;
              else sommet_arete_(som,1)=arete_haut_droite;
            }

          // arete YZ
          arete_bas_gauche=elem_aretes(mon_elem,5);
          arete_bas_droite=elem_aretes(mon_elem,8);
          arete_haut_gauche=elem_aretes(mon_elem,11);
          arete_haut_droite=elem_aretes(mon_elem,2);

          if (somy<cgy)
            {
              if (somz<cgz) sommet_arete_(som,0)=arete_bas_gauche;
              else sommet_arete_(som,0)=arete_haut_gauche;
            }
          else
            {
              if (somz<cgz) sommet_arete_(som,0)=arete_bas_droite;
              else sommet_arete_(som,0)=arete_haut_droite;
            }

        }

    }
}
static void ordonner_sommets_facettes(IntTab& facettes,
                                      const ArrOfInt& sommet_PE_owner,
                                      const ArrOfInt& sommet_num_owner)
{
  assert(facettes.dimension(1) == 3);
  const long nb_facettes = facettes.dimension(0);
  long i;
  for (i = 0; i < nb_facettes; i++)
    {
      long s0,s1,s2;
      s0 = facettes(i,0);
      s1 = facettes(i,1);
      s2 = facettes(i,2);
      if (s0 == s1) // Facette invalide a effacer
        s2 = s0;    // il faut qu'elle reste invalide apres permutation
      long pe0,pe1,pe2;
      pe0 = sommet_PE_owner[s0];
      pe1 = sommet_PE_owner[s1];
      pe2 = sommet_PE_owner[s2];
      long num0,num1,num2;
      num0 = sommet_num_owner[s0];
      num1 = sommet_num_owner[s1];
      num2 = sommet_num_owner[s2];
      while (pe1 < pe0 || (pe1 == pe0 && num1 < num0)
             || pe2 < pe0 || (pe2 == pe0 && num2 < num0))
        {
          long tmp;
          tmp = s0;
          s0 = s1;
          s1 = s2;
          s2 = tmp;
          tmp = pe0;
          pe0 = pe1;
          pe1 = pe2;
          pe2 = tmp;
          tmp = num0;
          num0 = num1;
          num1 = num2;
          num2 = tmp;
        }
      facettes(i,0) = s0;
      facettes(i,1) = s1;
      facettes(i,2) = s2;
    }
}

/*! @brief Sans changer les sommets existants ni la numerotation des facettes, on change le proprietaire des facettes de sorte que ce soit aussi le proprietaire
 *
 *  du premier sommet de la facette. Pour cela on doit eventuellement creer des
 *  sommets virtuels supplementaire et des facettes.
 *  Aucune facette n'est supprimee.
 *  Certaines facettes reelles sont creees,
 *  Certaines facettes reelles deviennent virtuelles.
 *  Le maillage retourne a l'etat minimal.
 *  Cette methode est utilisee lors de l'algorithme Marching-Cubes, du transport
 *  et du remaillage pour amener le maillage dans son etat conforme aux conventions.
 *
 * Precondition:
 *  Les membres suivants doivent etre valides :
 *  - Memes membres que creer_facettes_virtuelles
 *  - desc_facettes_ (remarque)
 *(remarque) desc_facettes_ doit etre un descripteur valide (correspondance entre
 *           elements distants et elements virtuels). En revanche, on suppose que
 *           la facette peut etre reelle sur n'importe quel processeur (pas forcement
 *           le proprietaire du premier sommet).
 * Postcondition:
 *  - Toutes les conditions qui definissent l'etat MINIMAL sont remplies.
 */
void Maillage_FT_Disc::corriger_proprietaires_facettes()
{
  //Process::Journal()<<"Maillage_FT_Disc::corriger_proprietaires_facettes nb_som="<<nb_sommets()<<finl;
  if (Comm_Group::check_enabled())
    check_mesh(1, 1 /* ne pas tester proprietaire facette */);
  assert(statut_ >= MINIMAL);

  static ArrOfIntFT facette_pe;
  static ArrOfIntFT liste_facettes;
  static ArrOfIntFT liste_PE;

  if (Objet_U::dimension == 3)
    {
      ordonner_sommets_facettes(facettes_,
                                sommet_PE_owner_,
                                sommet_num_owner_);
    }
  long i;
  const long moi = me();
  {
    const long nbfacettes = facettes_.dimension(0);
    // Quel est le proprietaire actuel des facettes ?
    // (d'apres le descripteur des facettes)
    facette_pe.resize_array(nbfacettes);
    desc_facettes_.remplir_element_pe(facette_pe);

    // A qui devraient etre les facettes ? (d'apres le numero du proprietaire
    // du premier sommet).
    // On construit facette_pe = -1 si le proprietaire est legitime
    //              facette_pe = numero du proprietaire legitime sinon.
    liste_facettes.resize_array(0);
    liste_PE.resize_array(0);
    for (i = 0; i < nbfacettes; i++)
      {
        long pe_actuel = facette_pe[i];
        long premier_sommet = facettes_(i, 0);
        long pe_legitime = sommet_PE_owner_[premier_sommet];
        if (pe_actuel == moi && pe_actuel != pe_legitime)
          {
            liste_facettes.append_array(i);
            liste_PE.append_array(pe_legitime);
            facette_pe[i] = pe_legitime;
          }
        else
          {
            facette_pe[i] = -1;
          }
      }
  }

  // Creation des facettes virtuelles sur le vrai proprietaire.
  // Le schema de comm est "les procs ayant des sommets virtuels
  // parlent aux procs ayant les sommets reels" (en effet,
  // si je suis le proprietaire illegitime d'une facette,
  // le sommet 0 est virtuel pour moi et reel pour le proprietaire
  // de la facette).
  {
    const Schema_Comm_FT& comm = desc_sommets_.schema_comm_inverse();
    creer_facettes_virtuelles(liste_facettes, liste_PE,
                              comm.get_send_pe_list(),
                              comm.get_recv_pe_list());
  }
  {
    // Simple resize, on n'a fait qu'ajouter des elements virtuels qui seront
    // ecrases par echange_espace_virtuel
    const long nbfacettes = facettes_.dimension(0);
    facette_pe.resize_array(nbfacettes);
    desc_facettes_.echange_espace_virtuel(facette_pe);
    // Mise a jour des espaces distants et virtuels (c'est la qu'on change
    // reelement de proprietaire)
    desc_facettes_.echanger_elements(facette_pe);
    // Mise a jour de facette_num_owner_
    facette_num_owner_.resize_array(nbfacettes);
    for (i = 0; i < nbfacettes; i++)
      facette_num_owner_[i] = i;
    desc_facettes_.echange_espace_virtuel(facette_num_owner_);
  }

  if (Comm_Group::check_enabled()) check_mesh();
  //Process::Journal()<<"FIN Maillage_FT_Disc::corriger_proprietaires_facettes nb_som="<<nb_sommets()<<finl;
}

/*! @brief Creation de facettes virtuelles sur le pe specifie.
 *
 * liste_facettes = liste de numeros de facettes reelles a envoyer sur le
 *                   processeur distant (la liste peut comporter des doublons)
 *                   Le processeur distant ne doit pas etre moi !
 *  liste_pe = numero du pe a qui il faut envoyer chaque facette
 *  comm = un schema ou send_pe_list contient les PEs mentionnes dans la liste.
 *
 *  Algo: Il faut d'abord creer les sommets des facettes, s'ils n'existent
 *  pas encore, puis creer les facettes.
 *  Cas general, le processeur A possede la facette (au sens du descripteur des
 *  facettes), le processeur B possede un sommet 's' de la facette, on veut envoyer
 *  la facette au processeur C. Il faut:
 *  A envoie a B le numero du sommet 's' et le numero du processeur C
 *  B cree sur C le sommet virtuel 's'
 *  A envoie a C les numeros des sommets 's' de la facette, C peut maintenant creer la facette.
 *
 */
void Maillage_FT_Disc::creer_facettes_virtuelles(const ArrOfInt& liste_facettes,
                                                 const ArrOfInt& liste_facettes_pe,
                                                 const ArrOfInt& facettes_send_pe_list,
                                                 const ArrOfInt& facettes_recv_pe_list)
{
  // Lorsque l'on arrive par Maillage_FT_IJK::creer_facettes_virtuelles
  // le tableau des compo_connexes doit etre a jour. On peut donc utiliser
  // le check_mesh == Maillage_FT_IJK::check_mesh.
  if (Comm_Group::check_enabled())
    check_mesh(1, 1 /* ne pas tester le proprietaire de la facette */);
  assert(statut_ >= MINIMAL);

  // On envoie des paquets contenant deux types de messages :
  // Messages envoyes au destinataire de la facette (pe0)
  //  token DATA_FACETTE
  //  numero du PE possedant le sommet0 (pe1)
  //  numero du sommet 0 de la facette sur le pe0
  //  numero du PE possedant le sommet1 (pe1)
  //  numero du sommet 1 de la facette sur pe1
  //[ numero du PE possedant le sommet2 (pe2)
  //  numero du sommet 2 de la facette sur pe2 ]  (en dimension 3)
  //
  // Messages envoyes au pe1 et pe2
  // pour leur dire de creer ces sommets sur le pe0
  //  token DATA_NOEUD
  //  numero du sommet
  //  numero du pe0
  static const long DATA_FACETTE = 0;
  static const long DATA_NOEUD = 1;
  long data_facette[8];
  long data_noeud[3];
  data_facette[0] = DATA_FACETTE;
  data_noeud[0] = DATA_NOEUD;

  const long n_proc = nproc();
  ArrOfIntFT liste_sommets;   // liste des sommets a envoyer
  ArrOfIntFT liste_sommets_pe;// pour chaque sommet, numero du pe a qui envoyer
  ArrOfIntFT recv_pe_flags;   // pour chaque pe, recoit-t-il des sommets ?
  ArrOfIntFT recv_pe_list;    // liste des pe qui recoivent des sommets
  ArrOfIntFT send_pe_flags;   // pour chaque pe, envoie-t-il des sommets ?
  ArrOfIntFT send_pe_list;    // liste des pe qui envoient des sommets

  long i;
  const long moi = me();

  liste_sommets.resize_array(0);
  liste_sommets_pe.resize_array(0);
  recv_pe_flags.resize_array(n_proc);
  recv_pe_flags = 0;
  send_pe_flags.resize_array(n_proc);
  send_pe_flags = 0;

  static Schema_Comm_FT comm;
  {
    // Construction d'un schema de comm contenant a la fois les destinataires
    // des facettes (liste_facettes_pe) et les proprietaires des sommets
    // des facettes.
    const Schema_Comm_FT& comm_sommets = desc_sommets_.schema_comm_inverse();
    send_pe_list = comm_sommets.get_send_pe_list();
    append_array_to_array(send_pe_list, facettes_send_pe_list);
    array_trier_retirer_doublons(send_pe_list);
    recv_pe_list = comm_sommets.get_recv_pe_list();
    append_array_to_array(recv_pe_list, facettes_recv_pe_list);
    array_trier_retirer_doublons(recv_pe_list);
    comm.set_send_recv_pe_list(send_pe_list, recv_pe_list);
  }
  comm.begin_comm();
  {
    const long nbfacettes = liste_facettes.size_array();
    Descripteur_FT& espace_distant = desc_facettes_.espace_distant();

    for (i = 0; i < nbfacettes; i++)
      {
        const long facette = liste_facettes[i];
        const long pe_destination = liste_facettes_pe[i];
        assert(pe_destination != moi);

        // Si la facette est deja virtuelle sur le pe_destination, inutile
        // de l'envoyer.
        if (espace_distant.contient_element(pe_destination, facette))
          continue;

        // Ajout de la facette a l'espace distant
        espace_distant.ajoute_element(pe_destination, facette);

        long n_data = 1;
        data_facette[n_data++] = facette; // Numero de la facette chez le proprietaire

        for (long j = 0; j < dimension; j++)
          {
            const long sommet = facettes_(facette, j);
            const long le_sommet_num_owner = sommet_num_owner_[sommet];
            const long sommet_pe = sommet_PE_owner_[sommet];

            // Donnees de la facette pour le destinataire de la facette
            data_facette[n_data++] = le_sommet_num_owner;
            data_facette[n_data++] = sommet_pe;
            // Demande d'envoi du sommet au proprietaire du sommet
            if (sommet_pe != pe_destination)
              {
                if (sommet_pe == moi)
                  {
                    // Si le sommet est a moi, j'empile tout de suite le sommet a envoyer
                    liste_sommets.append_array(le_sommet_num_owner);
                    liste_sommets_pe.append_array(pe_destination);
                    send_pe_flags[pe_destination] = 1;
                  }
                else
                  {
                    // Sinon, j'envoie la demande au pe qui possede le sommet
                    data_noeud[1] = le_sommet_num_owner;
                    data_noeud[2] = pe_destination;
                    comm.send_buffer(sommet_pe).put(data_noeud, 3);
                  }
              }
          }
        comm.send_buffer(pe_destination).put(data_facette, n_data);
      }
  }
  comm.echange_taille_et_messages();

  static ArrOfIntFT liste_facettes_pe_source;
  static ArrOfIntFT liste_facettes_num_owner;
  static ArrOfIntFT liste_facettes_sommets;
  static ArrOfIntFT liste_facettes_sommets_pe;
  liste_facettes_pe_source.resize_array(0);
  liste_facettes_num_owner.resize_array(0);
  liste_facettes_sommets.resize_array(0);
  liste_facettes_sommets_pe.resize_array(0);

  // Reception des donnees : liste des sommets a envoyer, liste
  // des facettes et des sommets a recevoir. On construit send/recv_pe_flags,
  // qui indiquent si on recoit ou si on envoie des sommets a chaque proc.
  // Cela servira a construire le schema de comm pour l'envoi des sommets
  // virtuels.
  {
    const ArrOfInt& pe_voisins = comm.get_recv_pe_list();
    const long nb_voisins = pe_voisins.size_array();
    for (i = 0; i < nb_voisins; i++)
      {
        const long pe_source = pe_voisins[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long token;
            buffer.get(& token, 1);
            if (buffer.eof())
              break;

            switch(token)
              {
              case DATA_FACETTE:
                {
                  buffer.get(data_facette, 2*dimension+1);
                  liste_facettes_pe_source.append_array(pe_source); // pe proprietaire de la facette
                  const long num_facette = data_facette[0];
                  liste_facettes_num_owner.append_array(num_facette);
                  for (long j = 0; j < dimension; j++)
                    {
                      // Numero du sommet sur le pe proprietaire
                      long sommet_num = data_facette[j*2+1];
                      // Numero du pe proprietaire du sommet
                      long sommet_pe = data_facette[j*2+2];
                      if (sommet_pe != moi)
                        recv_pe_flags[sommet_pe] = 1;
                      liste_facettes_sommets.append_array(sommet_num);
                      liste_facettes_sommets_pe.append_array(sommet_pe);
                    }
                  break;
                }
              case DATA_NOEUD:
                {
                  buffer.get(data_noeud, 2);
                  liste_sommets.append_array(data_noeud[0]); // numero du sommet a envoyer
                  liste_sommets_pe.append_array(data_noeud[1]);      // numero du pe destinataire
                  send_pe_flags[data_noeud[1]] = 1;
                  break;
                }
              default:
                assert(0);
              }
          }
      }
  }
  comm.end_comm();

  // Creation des sommets virtuels requis
  {
    // Construction du schema de comm. a partir de la liste des processeurs
    // a qui on envoie des donnees et de qui on recoit.
    static Schema_Comm_FT comm2;
    send_pe_list.resize_array(0);
    for (i = 0; i < n_proc; i++)
      if (send_pe_flags[i])
        send_pe_list.append_array(i);
    recv_pe_list.resize_array(0);
    for (i = 0; i < n_proc; i++)
      if (recv_pe_flags[i])
        recv_pe_list.append_array(i);
    comm2.set_send_recv_pe_list(send_pe_list, recv_pe_list);
    creer_sommets_virtuels(liste_sommets, liste_sommets_pe, comm2);
  }

  // Creation des facettes virtuelles
  {
    const long dimension3 = (dimension == 3);
    // Conversion des numeros distants des sommets en numeros locaux,
    // le resultat remplace le contenu de liste_facettes_sommets.
    convertir_numero_distant_local(desc_sommets_, sommet_num_owner_,
                                   liste_facettes_sommets,
                                   liste_facettes_sommets_pe,
                                   liste_facettes_sommets);
    const long nb_new_facettes = liste_facettes_pe_source.size_array();
    long nbfacettes = facettes_.dimension(0);
    facettes_.resize(nbfacettes + nb_new_facettes, dimension);
    facette_num_owner_.resize_array(nbfacettes + nb_new_facettes);
    for (i = 0; i < nb_new_facettes; i++)
      {
        facettes_(nbfacettes, 0) = liste_facettes_sommets[i*dimension];
        facettes_(nbfacettes, 1) = liste_facettes_sommets[i*dimension+1];
        if (dimension3)
          facettes_(nbfacettes, 2) = liste_facettes_sommets[i*dimension+2];
        desc_facettes_.espace_virtuel().ajoute_element(liste_facettes_pe_source[i],
                                                       nbfacettes);
        facette_num_owner_[nbfacettes] = liste_facettes_num_owner[i];
        nbfacettes++;
      }
  }
  desc_facettes_.espace_virtuel().calcul_liste_pe_voisins();
  desc_facettes_.espace_distant().calcul_liste_pe_voisins();
  desc_facettes_.calcul_schema_comm(facettes_.dimension(0));

  statut_ = MINIMAL;

  // Si cette methode est appelee par Maillage_FT_IJK::creer_facettes_virtuelles,
  // a ce stade, on a cree des facettes sans mettre a jour le tableau compo_connexe.
  // Le check_mesh a appeler ne doit donc pas verifier les compo_connexes.
  // La verification des compo_connexes est prevue a la fin de
  // Maillage_FT_IJK::creer_facettes_virtuelles
  if (Comm_Group::check_enabled())
    Maillage_FT_Disc::check_mesh(1, 1 /* ne pas tester le proprietaire de la facette */);
}

/*! @brief Echange des facettes dont les numeros sont fournis dans "liste_facettes" : Pour chaque facette a ajouter,
 *
 *  * On determine le processeur destinataire de la facette (c'est le
 *    processeur qui possede l'element_arrivee)
 *  * On ajoute la facette sur ce processeur si elle n'existe pas encore
 *    (creation d'une facette virtuelle et des sommets de la facette).
 *  * On envoie a ce processeur le numero de l'element_arrivee. Le numero est
 *    converti en numero local d'un element reel sur ce processeur.
 *  Cette methode est essentiellement utilisee dans le parcours de l'interface.
 *
 * Precondition: Le maillage doit etre dans l'etat minimal (en particulier,
 *  il doit respecter la convention "proprietaire facette=proprietaire 1er sommet").
 *
 * Trois categories de processeurs vont se parler:
 *  le processeur A qui connait le numero de la facette a envoyer,
 *  le processeur B qui possede la facette,
 *  le processeur C a qui on veut envoyer la facette.
 *
 * @param (liste_facettes) un tableau contenant des numeros de facettes reelles ou virtuelles, eventuellement avec doublons.
 * @param (liste_elem_arrivee) tableau de taille identique a liste_facettes, contenant pour chacune un numero d'element virtuel.
 * @param (facettes_recues_numfacettes) tableau rempli lors de l'echange avec la liste des facettes qui ont ete recues (meme si elles existaient deja et avec les doublons eventuels). Le contenu initial du tableau est efface. Aliasing ave liste_facettes autorise.
 * @param (facettes_recues_numelement) tableau rempli lors de l'echange: pour chaque facette recue, numero de l'element d'arrivee (c'est la traduction du numero d'element virtuel "liste_elem_arrivee" en numero d'element reel sur le processeur qui a recu la facette. Aliasing avec liste_elem_arrivee autorise.
 */
void Maillage_FT_Disc::echanger_facettes(const ArrOfInt& liste_facettes,
                                         const ArrOfInt& liste_elem_arrivee,
                                         ArrOfInt& facettes_recues_numfacettes,
                                         ArrOfInt& facettes_recues_numelement)
{
  if (Comm_Group::check_enabled()) check_mesh();
  const long moi = me();
  long i;
  const long nb_facettes_envoi = liste_facettes.size_array();

  // Numeros des elements euleriens associes a chaque facette (numero local
  // de l'element sur le pe qui possede cet element)
  static ArrOfIntFT liste_elem_arrivee_local;
  // Pour chaque element d'arrivee, determination du PE destination
  // et conversion du numero de l'element virtuel en numero local sur ce pe.
  static ArrOfIntFT liste_pe_dest;

  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine& le_dom = domaine_dis.domaine();
  const long nb_elem = le_dom.nb_elem(); // Nombre d'elements reels
  {
    liste_pe_dest.resize_array(nb_facettes_envoi);
    liste_elem_arrivee_local.resize_array(nb_facettes_envoi);

    const IntTab& elem_virt_pe_num = le_dom.elem_virt_pe_num();

    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long elem_voisin = liste_elem_arrivee[i];
        // Calcul de l'index de l'element dans elem_virt_pe_num (l'indice 0
        // correspond au premier element virtuel, soit l'element nb_elem)
        const long index = elem_voisin - nb_elem;
        // Numero du pe qui possede l'element voisin
        const long pe = elem_virt_pe_num(index, 0);
        // Numero local de l'element sur ce pe.
        const long elem_local = elem_virt_pe_num(index, 1);
        liste_pe_dest[i] = pe;
        liste_elem_arrivee_local[i] = elem_local;
      }
  }
  // =======================================================================
  //                         Communication A -> B
  //
  // Envoi de la liste des facettes a transmettre au processeur proprietaire de la
  // facette (A envoie les numeros de facettes a B)
  // On remplit facettes_to_send (liste des facettes que B doit envoyer a C)
  //            facettes_pe_dest (numero du processeur C a qui il faut envoyer)
  //            BtoC_send_pe_flags (drapeaux des processeurs a qui on envoie)
  // Attention : si B=C, on n'envoie pas la facette.
  static ArrOfIntFT facettes_to_send;
  static ArrOfIntFT facettes_pe_dest;
  // Lors de la communication finale entre les procs B et les procs C,
  // ces drapeaux indiquent si on recoit des donnees et si on en envoie
  // a chacun des processeurs.
  static long nbproc = Process::nproc();
  static ArrOfIntFT BtoC_send_pe_flags(nbproc);
  static ArrOfIntFT BtoC_recv_pe_flags(nbproc);
  BtoC_send_pe_flags = 0;
  facettes_to_send.resize_array(0);
  facettes_pe_dest.resize_array(0);
  {
    // A et B sont voisins au sens des espaces distants/virtuels des facettes
    // (un processeur chez qui la facette est virtuelle envoie des donnees au
    //  processeur chez qui elle est reelle => schema_comm_inverse)
    const Schema_Comm_FT& comm = desc_facettes_.schema_comm_inverse();
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = liste_facettes[i];
        const long PE_destinataire = liste_pe_dest[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        if (PE_proprietaire == moi)
          {
            facettes_to_send.append_array(facette_num_owner);
            facettes_pe_dest.append_array(PE_destinataire);
            BtoC_send_pe_flags[PE_destinataire] = 1;
          }
        else
          {
            if (PE_destinataire != PE_proprietaire)
              comm.send_buffer(PE_proprietaire) << facette_num_owner << PE_destinataire;
          }
      }
    comm.echange_taille_et_messages();
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette, PE_destinataire;
            buffer >> facette >> PE_destinataire;
            if (buffer.eof())
              break;
            facettes_to_send.append_array(facette);
            facettes_pe_dest.append_array(PE_destinataire);
            BtoC_send_pe_flags[PE_destinataire] = 1;
          }
      }
    comm.end_comm();
  }

  // =======================================================================
  //                         Communication A -> C
  //
  // A envoie a C la liste la liste des facettes, l'expediteur et le numero
  // de l'element d'entree.
  // On remplit facettes_recues_numfacettes (numerotation temporaire)
  //            facettes_recues_numelement
  //            BtoC_recv_pe_flags (drapeaux des procs de qui on recoit des facettes).
  // Attention : si B=C on ne recoit pas de donnees de B

  // Pour chaque facette recue, numero du PE proprietaire de la facette:
  static ArrOfIntFT nouvelles_facettes_pe_proprietaire;
  {
    nouvelles_facettes_pe_proprietaire.resize_array(0);
    BtoC_recv_pe_flags = 0;

    // A et C sont voisins au sens du maillage eulerien
    const Schema_Comm_FT& comm = schema_comm_domaine_;
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = liste_facettes[i];
        const long PE_destinataire = liste_pe_dest[i];
        const long element_arrivee = liste_elem_arrivee_local[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        assert (PE_destinataire != moi);
        comm.send_buffer(PE_destinataire) << facette_num_owner
                                          << PE_proprietaire
                                          << element_arrivee;
      }
    comm.echange_taille_et_messages();
    // On efface facettes_recues_numfacettes et numelement.
    // Pour le cas ou il y aurait aliasing des parametres, il ne faut plus
    // utiliser liste_facettes et liste_elem_arrivee.
    facettes_recues_numfacettes.resize_array(0);
    facettes_recues_numelement.resize_array(0);
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette = -1, PE_proprietaire = -1, element_arrivee = -1;
            buffer >> facette >> PE_proprietaire >> element_arrivee;
            if (buffer.eof())
              break;
            assert(facette >= 0 && PE_proprietaire >= 0 && element_arrivee >= 0);
            assert(element_arrivee < nb_elem); // Ce doit etre un element reel
            // Enregistre le numero de la nouvelle facette et l'element d'arrivee.
            // Pour l'instant c'est le numero de la facette sur le proprietaire.
            facettes_recues_numfacettes.append_array(facette);
            nouvelles_facettes_pe_proprietaire.append_array(PE_proprietaire);
            facettes_recues_numelement.append_array(element_arrivee);
            // Enregistre que l'on va recevoir des donnees du PE_proprietaire lors
            // de la derniere communication B->C
            if (PE_proprietaire != moi)
              BtoC_recv_pe_flags[PE_proprietaire] = 1;
          }
      }
    comm.end_comm();
  }

  // =======================================================================
  //                         Communication B -> C
  //
  // Construction du schema de communication pour envoyer les facettes reelles
  // au destinataire (B envoie les facettes a C) : il faut determiner
  // a qui j'envoie et de qui de recois.
  {
    static ArrOfIntFT send_pe_list;
    static ArrOfIntFT recv_pe_list;
    send_pe_list.resize_array(0);
    recv_pe_list.resize_array(0);
    for (i = 0; i < nbproc; i++)
      if (BtoC_send_pe_flags[i])
        send_pe_list.append_array(i);
    for (i = 0; i < nbproc; i++)
      if (BtoC_recv_pe_flags[i])
        recv_pe_list.append_array(i);
    // Envoi des facettes
    creer_facettes_virtuelles(facettes_to_send, facettes_pe_dest,
                              send_pe_list, recv_pe_list);
  }
  // Conversion du numero de la facette en numero local dans facettes_recues_numfacettes
  convertir_numero_distant_local(desc_facettes_,
                                 facette_num_owner_,
                                 facettes_recues_numfacettes,
                                 nouvelles_facettes_pe_proprietaire,
                                 facettes_recues_numfacettes);
  if (Comm_Group::check_enabled()) check_mesh();
}
// Description:
// idem echanger_facettes mais pour les faces
void Maillage_FT_Disc::echanger_facettes_face_x(const ArrOfInt& liste_facettes,
                                                const ArrOfInt& liste_face_arrivee,
                                                ArrOfInt& facettes_recues_numfacettes,
                                                ArrOfInt& facettes_recues_numface)
{
  if (Comm_Group::check_enabled()) check_mesh();
  const long moi = me();
  long i=0;
  long indice=0;
  long nb_facettes_envoi = liste_facettes.size_array();
  const long nb_facettes_envoi_init=liste_facettes.size_array();
  // Numeros des faces euleriennes associes a chaque facette (numero local
  // de la face sur le pe qui possede cet element)
  static ArrOfIntFT liste_facex_arrivee_local;
  ArrOfIntFT la_liste_facettes;
  la_liste_facettes.resize_array(nb_facettes_envoi);
  // Pour chaque face d'arrivee, determination du PE destination
  // et conversion du numero de la face virtuelle en numero local sur ce pe.
  static ArrOfIntFT liste_pe_dest_x;
  const Desc_Structure_FT& desc_facettes_const = desc_facettes();
  const Descripteur_FT& espace_distant = desc_facettes_const.espace_distant();
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& le_domaine_vf = ref_cast(Domaine_VF,domaine_dis.valeur());
  const ArrOfInt& faces_doubles = le_domaine_vf.faces_doubles();
  const IntTab& faces_doubles_pe_num = le_domaine_vf.faces_doubles_pe_num();
  const IntTab& faces_doubles_virt_pe_num = le_domaine_vf.faces_doubles_virt_pe_num();
  const long nb_face = le_domaine_vf.nb_faces(); // Nombre de faces reelles
  {
    liste_pe_dest_x.resize_array(nb_facettes_envoi);
    liste_facex_arrivee_local.resize_array(nb_facettes_envoi);

    const IntTab& face_virt_pe_num = le_domaine_vf.face_virt_pe_num();

    while (i<nb_facettes_envoi_init)
      {
        const long face_voisine = liste_face_arrivee[i];
        const long face_voisine_double = (faces_doubles_virt_pe_num(face_voisine,0)>0);
        if (faces_doubles(face_voisine) && face_voisine < nb_face)
          {
            liste_pe_dest_x[indice] = faces_doubles_pe_num(face_voisine,0);
            liste_facex_arrivee_local[indice] = faces_doubles_pe_num(face_voisine,1);
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;
          }
        else if (face_voisine_double && face_voisine >= nb_face)
          {

            // si la fa7 n'existe pas, on ne l'envoie pas
            // on est oblige de faire ca car complique d'identifier le pe auquel on veut envoyer la fa7 dans le cas suivant :

            //			 | a droite de cette zone, proc 2 uniquement  (et proc 1 a gauche)
            //    proc1  | proc 2
            //	-----------------------------------------------
            //	   |  -> |	 /  |     |
            //	------------/-----------------------------------  en dessous de cette zone : proc 0 uniquement
            //     |  -->| /    |     |     proc 0
            //  ------------------------------------------------   / : fa7, --> face parcourue, -> face que l'on veut parcourir
            // La face "->" que l'on veut parcourir appartient a la fois au proc 1 et au proc 2, c'est une face double
            // On a simplement l'info que la fa7 traverse cette face, mais on ne sait pas de quel cote ==> on ne sait pas a qui l'envoyer
            // On doit donc checker si le proc 1 et le proc 2 contiennent la fa7
            // Lors du parcours des faces ON NE CREE PAS DE NOUVELLES FA7, moralement, on demande simplement au proc de calculer l'intersection entre une fa7 connue
            // et une de ses faces
            // Ainsi, dans le cas present, le proc 1 ne connait pas la fa7 et le proc 2 la connait
            // c'est donc uniquement au proc 2 que l'on envoie la fa7

            const long fa7=liste_facettes[i];
            const long pe1=faces_doubles_virt_pe_num(face_voisine,0);
            const long face_1=faces_doubles_virt_pe_num(face_voisine,1);
            const long pe2=faces_doubles_virt_pe_num(face_voisine,2);
            const long face_2=faces_doubles_virt_pe_num(face_voisine,3);

            long fa7_reelle_pour_pe1=1;
            fa7_reelle_pour_pe1 = (sommet_PE_owner_(facettes_(fa7,0))==pe1) ? 1 : 0;
            long fa7_reelle_pour_pe2=1;
            fa7_reelle_pour_pe2 = (sommet_PE_owner_(facettes_(fa7,0))==pe2) ? 1 : 0;
            const long fa7_existante_pour_pe1= fa7_reelle_pour_pe1 || espace_distant.contient_element(pe1, fa7);
            const long fa7_existante_pour_pe2= fa7_reelle_pour_pe2 || espace_distant.contient_element(pe2, fa7);
            if (fa7_existante_pour_pe1)
              {
                liste_pe_dest_x[indice] = pe1;
                liste_facex_arrivee_local[indice] = face_1;
                la_liste_facettes[indice]=fa7;
                indice++;
              }
            i++;
            if (fa7_existante_pour_pe2)
              {
                liste_pe_dest_x[indice] = pe2;
                liste_facex_arrivee_local[indice] = face_2;
                la_liste_facettes[indice]=fa7;
                indice++;
              }
          }

        else
          {
            // Calcul de l'index de la face dans face_virt_pe_num (l'indice 0
            // correspond a la premiere face virtuelle, soit la face nb_face)
            const long index = face_voisine - nb_face;
            const long pe = face_virt_pe_num(index, 0);
            const long face_locale = face_virt_pe_num(index, 1);
            liste_pe_dest_x[indice] = pe;
            liste_facex_arrivee_local[indice] = face_locale;
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;
          }
        i++;
      }

    liste_pe_dest_x.resize_array(indice);
    liste_facex_arrivee_local.resize_array(indice);
    la_liste_facettes.resize_array(indice);
  }
  nb_facettes_envoi=indice;
  // =======================================================================
  //                         Communication A -> B
  //
  // Envoi de la liste des facettes a transmettre au processeur proprietaire de la
  // facette (A envoie les numeros de facettes a B)
  // On remplit facettes_to_send (liste des facettes que B doit envoyer a C)
  //            facettes_pe_dest (numero du processeur C a qui il faut envoyer)
  //            BtoC_send_pe_flags (drapeaux des processeurs a qui on envoie)
  // Attention : si B=C, on n'envoie pas la facette.
  static ArrOfIntFT facettes_to_send_x;
  static ArrOfIntFT facettes_pe_dest_x;
  // Lors de la communication finale entre les procs B et les procs C,
  // ces drapeaux indiquent si on recoit des donnees et si on en envoie
  // a chacun des processeurs.
  static long nbproc = Process::nproc();
  static ArrOfIntFT BtoC_send_pe_flags_x(nbproc);
  static ArrOfIntFT BtoC_recv_pe_flags_x(nbproc);
  BtoC_send_pe_flags_x = 0;
  facettes_to_send_x.resize_array(0);
  facettes_pe_dest_x.resize_array(0);
  {
    // A et B sont voisins au sens des espaces distants/virtuels des facettes
    // (un processeur chez qui la facette est virtuelle envoie des donnees au
    //  processeur chez qui elle est reelle => schema_comm_inverse)
    const Schema_Comm_FT& comm = desc_facettes_.schema_comm_inverse();
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_x[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        if (PE_proprietaire == moi)
          {
            facettes_to_send_x.append_array(facette_num_owner);
            facettes_pe_dest_x.append_array(PE_destinataire);
            BtoC_send_pe_flags_x[PE_destinataire] = 1;
          }
        else
          {
            if (PE_destinataire != PE_proprietaire)
              comm.send_buffer(PE_proprietaire) << facette_num_owner << PE_destinataire;
          }
      }
    comm.echange_taille_et_messages();
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette, PE_destinataire;
            buffer >> facette >> PE_destinataire;
            if (buffer.eof())
              break;
            facettes_to_send_x.append_array(facette);
            facettes_pe_dest_x.append_array(PE_destinataire);
            BtoC_send_pe_flags_x[PE_destinataire] = 1;
          }
      }
    comm.end_comm();
  }
  //Cerr << "Fin com A to B" << finl;
  // =======================================================================
  //                         Communication A -> C
  //
  // A envoie a C la liste la liste des facettes, l'expediteur et le numero
  // de l'element d'entree.
  // On remplit facettes_recues_numfacettes (numerotation temporaire)
  //            facettes_recues_numelement
  //            BtoC_recv_pe_flags (drapeaux des procs de qui on recoit des facettes).
  // Attention : si B=C on ne recoit pas de donnees de B

  // Pour chaque facette recue, numero du PE proprietaire de la facette:
  static ArrOfIntFT nouvelles_facettes_pe_proprietaire_x;
  {
    nouvelles_facettes_pe_proprietaire_x.resize_array(0);
    BtoC_recv_pe_flags_x = 0;

    // A et C sont voisins au sens du maillage eulerien
    const Schema_Comm_FT& comm = schema_comm_domaine_;
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_x[i];
        const long face_arrivee = liste_facex_arrivee_local[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        assert (PE_destinataire != moi);
        comm.send_buffer(PE_destinataire) << facette_num_owner
                                          << PE_proprietaire
                                          << face_arrivee;
      }
    comm.echange_taille_et_messages();
    // On efface facettes_recues_numfacettes et numelement.
    // Pour le cas ou il y aurait aliasing des parametres, il ne faut plus
    // utiliser liste_facettes et liste_elem_arrivee.
    facettes_recues_numfacettes.resize_array(0);
    facettes_recues_numface.resize_array(0);
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette = -1, PE_proprietaire = -1, face_arrivee = -1;
            buffer >> facette >> PE_proprietaire >> face_arrivee;
            if (buffer.eof())
              break;
            assert(facette >= 0 && PE_proprietaire >= 0 && face_arrivee >= 0);
            assert(face_arrivee < nb_face); // Ce doit etre un element reel
            // Enregistre le numero de la nouvelle facette et l'element d'arrivee.
            // Pour l'instant c'est le numero de la facette sur le proprietaire.
            facettes_recues_numfacettes.append_array(facette);
            nouvelles_facettes_pe_proprietaire_x.append_array(PE_proprietaire);
            facettes_recues_numface.append_array(face_arrivee);
            // Enregistre que l'on va recevoir des donnees du PE_proprietaire lors
            // de la derniere communication B->C
            if (PE_proprietaire != moi)
              BtoC_recv_pe_flags_x[PE_proprietaire] = 1;
          }
      }
    comm.end_comm();
  }
  //Cerr << "Fin com A to C" << finl;
  // =====
  // ==================================================================
  //                         Communication B -> C
  //
  // Construction du schema de communication pour envoyer les facettes reelles
  // au destinataire (B envoie les facettes a C) : il faut determiner
  // a qui j'envoie et de qui de recois.
  {
    static ArrOfIntFT send_pe_list_x;
    static ArrOfIntFT recv_pe_list_x;
    send_pe_list_x.resize_array(0);
    recv_pe_list_x.resize_array(0);
    for (i = 0; i < nbproc; i++)
      if (BtoC_send_pe_flags_x[i])
        send_pe_list_x.append_array(i);
    for (i = 0; i < nbproc; i++)
      if (BtoC_recv_pe_flags_x[i])
        recv_pe_list_x.append_array(i);
    // Envoi des facettes
    creer_facettes_virtuelles(facettes_to_send_x, facettes_pe_dest_x,
                              send_pe_list_x, recv_pe_list_x);
  }
  //Cerr << "Fin com B to C" << finl;
  // Conversion du numero de la facette en numero local dans facettes_recues_numfacettes
  convertir_numero_distant_local(desc_facettes_,
                                 facette_num_owner_,
                                 facettes_recues_numfacettes,
                                 nouvelles_facettes_pe_proprietaire_x,
                                 facettes_recues_numfacettes);
  if (Comm_Group::check_enabled()) check_mesh();
}

// Description:
// idem echanger_facettes mais pour les faces
void Maillage_FT_Disc::echanger_facettes_face_y(const ArrOfInt& liste_facettes,
                                                const ArrOfInt& liste_face_arrivee,
                                                ArrOfInt& facettes_recues_numfacettes,
                                                ArrOfInt& facettes_recues_numface)
{
  if (Comm_Group::check_enabled()) check_mesh();
  const long moi = me();
  long i=0;
  long indice=0;
  long nb_facettes_envoi = liste_facettes.size_array();
  const long nb_facettes_envoi_init=liste_facettes.size_array();
  // Numeros des faces euleriennes associes a chaque facette (numero local
  // de la face sur le pe qui possede cet element)
  static ArrOfIntFT liste_facey_arrivee_local;
  ArrOfIntFT la_liste_facettes;
  la_liste_facettes.resize_array(nb_facettes_envoi);
  // Pour chaque face d'arrivee, determination du PE destination
  // et conversion du numero de la face virtuelle en numero local sur ce pe.
  static ArrOfIntFT liste_pe_dest_y;
  const Desc_Structure_FT& desc_facettes_const = desc_facettes();
  const Descripteur_FT& espace_distant = desc_facettes_const.espace_distant();
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& le_domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());

  const ArrOfInt& faces_doubles = le_domaine_vf.faces_doubles();
  const IntTab& faces_doubles_pe_num = le_domaine_vf.faces_doubles_pe_num();
  const IntTab& faces_doubles_virt_pe_num = le_domaine_vf.faces_doubles_virt_pe_num();

  const long nb_face = le_domaine_vf.nb_faces(); // Nombre de faces reelles
  {
    liste_pe_dest_y.resize_array(nb_facettes_envoi);
    liste_facey_arrivee_local.resize_array(nb_facettes_envoi);

    const IntTab& face_virt_pe_num = le_domaine_vf.face_virt_pe_num();
    //Cerr << "liste_face_arrivee.size_array() " << liste_face_arrivee.size_array() << finl;
    while (i < nb_facettes_envoi_init)
      {
        const long face_voisine = liste_face_arrivee[i];
        const long face_voisine_double = (faces_doubles_virt_pe_num(face_voisine,0)>0);
        if (faces_doubles(face_voisine) && face_voisine < nb_face)
          {
            liste_pe_dest_y[indice] = faces_doubles_pe_num(face_voisine,0);
            liste_facey_arrivee_local[indice] = faces_doubles_pe_num(face_voisine,1);
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;
          }

        else if (face_voisine_double && face_voisine >= nb_face)
          {
            const long fa7=liste_facettes[i];
            const long pe1=faces_doubles_virt_pe_num(face_voisine,0);
            const long face_1=faces_doubles_virt_pe_num(face_voisine,1);
            const long pe2=faces_doubles_virt_pe_num(face_voisine,2);
            const long face_2=faces_doubles_virt_pe_num(face_voisine,3);

            long fa7_reelle_pour_pe1=1;
            fa7_reelle_pour_pe1 = (sommet_PE_owner_(facettes_(fa7,0))==pe1) ? 1 : 0;
            long fa7_reelle_pour_pe2=1;
            fa7_reelle_pour_pe2 = (sommet_PE_owner_(facettes_(fa7,0))==pe2) ? 1 : 0;
            const long fa7_existante_pour_pe1= fa7_reelle_pour_pe1 || espace_distant.contient_element(pe1, fa7);
            const long fa7_existante_pour_pe2= fa7_reelle_pour_pe2 || espace_distant.contient_element(pe2, fa7);

            if (fa7_existante_pour_pe1)
              {
                liste_pe_dest_y[indice] = pe1;
                liste_facey_arrivee_local[indice] = face_1;
                la_liste_facettes[indice]=fa7;
                indice++;
              }
            i++;
            if (fa7_existante_pour_pe2)
              {
                liste_pe_dest_y[indice] = pe2;
                liste_facey_arrivee_local[indice] = face_2;
                la_liste_facettes[indice]=fa7;
                indice++;
              }

          }

        else
          {
            // Calcul de l'index de la face dans face_virt_pe_num (l'indice 0
            // correspond a la premiere face virtuelle, soit la face nb_face)
            const long index = face_voisine - nb_face;
            // Numero du pe qui possede la face voisine
            const long pe = face_virt_pe_num(index, 0);
            // Numero local de la face sur ce pe.
            const long face_locale = face_virt_pe_num(index, 1);
            liste_pe_dest_y[indice] = pe;
            liste_facey_arrivee_local[indice] = face_locale;
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;
          }
        i++;
      }
    liste_pe_dest_y.resize_array(indice);
    liste_facey_arrivee_local.resize_array(indice);
    la_liste_facettes.resize_array(indice);
    //Cerr << "liste_facey_arrivee_local.size_array() " << liste_facey_arrivee_local.size_array() << finl;
  }
  nb_facettes_envoi=indice;
  // ===================================
  // =======================================================================
  //                         Communication A -> B
  //
  // Envoi de la liste des facettes a transmettre au processeur proprietaire de la
  // facette (A envoie les numeros de facettes a B)
  // On remplit facettes_to_send (liste des facettes que B doit envoyer a C)
  //            facettes_pe_dest (numero du processeur C a qui il faut envoyer)
  //            BtoC_send_pe_flags (drapeaux des processeurs a qui on envoie)
  // Attention : si B=C, on n'envoie pas la facette.
  static ArrOfIntFT facettes_to_send_y;
  static ArrOfIntFT facettes_pe_dest_y;
  // Lors de la communication finale entre les procs B et les procs C,
  // ces drapeaux indiquent si on recoit des donnees et si on en envoie
  // a chacun des processeurs.
  static long nbproc = Process::nproc();
  static ArrOfIntFT BtoC_send_pe_flags_y(nbproc);
  static ArrOfIntFT BtoC_recv_pe_flags_y(nbproc);
  BtoC_send_pe_flags_y = 0;
  facettes_to_send_y.resize_array(0);
  facettes_pe_dest_y.resize_array(0);
  {
    // A et B sont voisins au sens des espaces distants/virtuels des facettes
    // (un processeur chez qui la facette est virtuelle envoie des donnees au
    //  processeur chez qui elle est reelle => schema_comm_inverse)
    const Schema_Comm_FT& comm = desc_facettes_.schema_comm_inverse();
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_y[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        if (PE_proprietaire == moi)
          {
            facettes_to_send_y.append_array(facette_num_owner);
            facettes_pe_dest_y.append_array(PE_destinataire);
            BtoC_send_pe_flags_y[PE_destinataire] = 1;
          }
        else
          {
            if (PE_destinataire != PE_proprietaire)
              comm.send_buffer(PE_proprietaire) << facette_num_owner << PE_destinataire;
          }
      }
    comm.echange_taille_et_messages();
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette, PE_destinataire;
            buffer >> facette >> PE_destinataire;
            if (buffer.eof())
              break;
            facettes_to_send_y.append_array(facette);
            facettes_pe_dest_y.append_array(PE_destinataire);
            BtoC_send_pe_flags_y[PE_destinataire] = 1;
          }
      }
    comm.end_comm();
  }
  // =======================================================================
  //                         Communication A -> C
  //
  // A envoie a C la liste la liste des facettes, l'expediteur et le numero
  // de l'element d'entree.
  // On remplit facettes_recues_numfacettes (numerotation temporaire)
  //            facettes_recues_numelement
  //            BtoC_recv_pe_flags (drapeaux des procs de qui on recoit des facettes).
  // Attention : si B=C on ne recoit pas de donnees de B

  // Pour chaque facette recue, numero du PE proprietaire de la facette:
  static ArrOfIntFT nouvelles_facettes_pe_proprietaire_y;
  {
    nouvelles_facettes_pe_proprietaire_y.resize_array(0);
    BtoC_recv_pe_flags_y = 0;

    // A et C sont voisins au sens du maillage eulerien
    const Schema_Comm_FT& comm = schema_comm_domaine_;
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_y[i];
        const long face_arrivee = liste_facey_arrivee_local[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        assert (PE_destinataire != moi);
        comm.send_buffer(PE_destinataire) << facette_num_owner
                                          << PE_proprietaire
                                          << face_arrivee;
      }
    comm.echange_taille_et_messages();
    // On efface facettes_recues_numfacettes et numelement.
    // Pour le cas ou il y aurait aliasing des parametres, il ne faut plus
    // utiliser liste_facettes et liste_elem_arrivee.
    facettes_recues_numfacettes.resize_array(0);
    facettes_recues_numface.resize_array(0);
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette = -1, PE_proprietaire = -1, face_arrivee = -1;
            buffer >> facette >> PE_proprietaire >> face_arrivee;
            if (buffer.eof())
              break;
            assert(facette >= 0 && PE_proprietaire >= 0 && face_arrivee >= 0);
            assert(face_arrivee < nb_face); // Ce doit etre un element reel
            // Enregistre le numero de la nouvelle facette et l'element d'arrivee.
            // Pour l'instant c'est le numero de la facette sur le proprietaire.
            facettes_recues_numfacettes.append_array(facette);
            nouvelles_facettes_pe_proprietaire_y.append_array(PE_proprietaire);
            facettes_recues_numface.append_array(face_arrivee);
            // Enregistre que l'on va recevoir des donnees du PE_proprietaire lors
            // de la derniere communication B->C
            if (PE_proprietaire != moi)
              BtoC_recv_pe_flags_y[PE_proprietaire] = 1;
          }
      }
    comm.end_comm();
  }
  // =======================================================================
  //                         Communication B -> C
  //
  // Construction du schema de communication pour envoyer les facettes reelles
  // au destinataire (B envoie les facettes a C) : il faut determiner
  // a qui j'envoie et de qui de recois.
  {
    static ArrOfIntFT send_pe_list_y;
    static ArrOfIntFT recv_pe_list_y;
    send_pe_list_y.resize_array(0);
    recv_pe_list_y.resize_array(0);
    for (i = 0; i < nbproc; i++)
      if (BtoC_send_pe_flags_y[i])
        send_pe_list_y.append_array(i);
    for (i = 0; i < nbproc; i++)
      if (BtoC_recv_pe_flags_y[i])
        recv_pe_list_y.append_array(i);
    // Envoi des facettes
    creer_facettes_virtuelles(facettes_to_send_y, facettes_pe_dest_y,
                              send_pe_list_y, recv_pe_list_y);
  }
  // Conversion du numero de la facette en numero local dans facettes_recues_numfacettes
  convertir_numero_distant_local(desc_facettes_,
                                 facette_num_owner_,
                                 facettes_recues_numfacettes,
                                 nouvelles_facettes_pe_proprietaire_y,
                                 facettes_recues_numfacettes);
  if (Comm_Group::check_enabled()) check_mesh();
}

// Description:
// idem echanger_facettes mais pour les faces
void Maillage_FT_Disc::echanger_facettes_face_z(const ArrOfInt& liste_facettes,
                                                const ArrOfInt& liste_face_arrivee,
                                                ArrOfInt& facettes_recues_numfacettes,
                                                ArrOfInt& facettes_recues_numface)
{
  if (Comm_Group::check_enabled()) check_mesh();
  const long moi = me();
  long i=0;
  long indice=0;
  long nb_facettes_envoi = liste_facettes.size_array();
  const long nb_facettes_envoi_init=liste_facettes.size_array();
  // Numeros des faces euleriennes associes a chaque facette (numero local
  // de la face sur le pe qui possede cet element)
  static ArrOfIntFT liste_facez_arrivee_local;
  ArrOfIntFT la_liste_facettes;
  la_liste_facettes.resize_array(nb_facettes_envoi);
  // Pour chaque face d'arrivee, determination du PE destination
  // et conversion du numero de la face virtuelle en numero local sur ce pe.
  static ArrOfIntFT liste_pe_dest_z;
  const Desc_Structure_FT& desc_facettes_const = desc_facettes();
  const Descripteur_FT& espace_distant = desc_facettes_const.espace_distant();
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& le_domaine_vf = ref_cast(Domaine_VF,domaine_dis.valeur());

  const ArrOfInt& faces_doubles = le_domaine_vf.faces_doubles();
  const IntTab& faces_doubles_pe_num = le_domaine_vf.faces_doubles_pe_num();
  const IntTab& faces_doubles_virt_pe_num = le_domaine_vf.faces_doubles_virt_pe_num();

  const long nb_face = le_domaine_vf.nb_faces(); // Nombre de faces reelles
  {
    liste_pe_dest_z.resize_array(nb_facettes_envoi);
    liste_facez_arrivee_local.resize_array(nb_facettes_envoi);

    const IntTab& face_virt_pe_num = le_domaine_vf.face_virt_pe_num();
    while (i < nb_facettes_envoi_init)
      {
        const long face_voisine = liste_face_arrivee[i];
        const long face_voisine_double = (faces_doubles_virt_pe_num(face_voisine,0)>0);
        if (faces_doubles(face_voisine) && face_voisine < nb_face)
          {
            //Cerr << "face_double " << le_domaine_vf.xv(face_voisine,0) << " " << le_domaine_vf.xv(face_voisine,1) << " " << le_domaine_vf.xv(face_voisine,2) << finl;
            liste_pe_dest_z[indice] = faces_doubles_pe_num(face_voisine,0);
            liste_facez_arrivee_local[indice] = faces_doubles_pe_num(face_voisine,1);
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;
          }

        else if (face_voisine_double && face_voisine >= nb_face)
          {
            const long fa7=liste_facettes[i];
            const long pe1=faces_doubles_virt_pe_num(face_voisine,0);
            const long face_1=faces_doubles_virt_pe_num(face_voisine,1);
            const long pe2=faces_doubles_virt_pe_num(face_voisine,2);
            const long face_2=faces_doubles_virt_pe_num(face_voisine,3);

            long fa7_reelle_pour_pe1=1;
            fa7_reelle_pour_pe1 = (sommet_PE_owner_(facettes_(fa7,0))==pe1) ? 1 : 0;
            long fa7_reelle_pour_pe2=1;
            fa7_reelle_pour_pe2 = (sommet_PE_owner_(facettes_(fa7,0))==pe2) ? 1 : 0;
            const long fa7_existante_pour_pe1= fa7_reelle_pour_pe1 || espace_distant.contient_element(pe1, fa7);
            const long fa7_existante_pour_pe2= fa7_reelle_pour_pe2 || espace_distant.contient_element(pe2, fa7);

            if (fa7_existante_pour_pe1)
              {
                liste_pe_dest_z[indice] = pe1;
                liste_facez_arrivee_local[indice] = face_1;
                la_liste_facettes[indice]=fa7;
                indice++;
              }
            i++;
            if (fa7_existante_pour_pe2)
              {
                liste_pe_dest_z[indice] = pe2;
                liste_facez_arrivee_local[indice] = face_2;
                la_liste_facettes[indice]=fa7;
                indice++;
              }
          }
        else
          {
            // Calcul de l'index de la face dans face_virt_pe_num (l'indice 0
            // correspond a la premiere face virtuelle, soit la face nb_face)
            const long index = face_voisine - nb_face;
            // Numero du pe qui possede la face voisine
            const long pe = face_virt_pe_num(index, 0);
            // Numero local de la face sur ce pe.
            const long face_locale = face_virt_pe_num(index, 1);
            liste_pe_dest_z[indice] = pe;
            liste_facez_arrivee_local[indice] = face_locale;
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;
          }
        i++;
      }

    liste_pe_dest_z.resize_array(indice);
    liste_facez_arrivee_local.resize_array(indice);
    la_liste_facettes.resize_array(indice);
  }
  nb_facettes_envoi=indice;
  // =======================================================================
  //                         Communication A -> B
  //
  // Envoi de la liste des facettes a transmettre au processeur proprietaire de la
  // facette (A envoie les numeros de facettes a B)
  // On remplit facettes_to_send (liste des facettes que B doit envoyer a C)
  //            facettes_pe_dest (numero du processeur C a qui il faut envoyer)
  //            BtoC_send_pe_flags (drapeaux des processeurs a qui on envoie)
  // Attention : si B=C, on n'envoie pas la facette.
  static ArrOfIntFT facettes_to_send_z;
  static ArrOfIntFT facettes_pe_dest_z;
  // Lors de la communication finale entre les procs B et les procs C,
  // ces drapeaux indiquent si on recoit des donnees et si on en envoie
  // a chacun des processeurs.
  static long nbproc = Process::nproc();
  static ArrOfIntFT BtoC_send_pe_flags_z(nbproc);
  static ArrOfIntFT BtoC_recv_pe_flags_z(nbproc);
  BtoC_send_pe_flags_z = 0;
  facettes_to_send_z.resize_array(0);
  facettes_pe_dest_z.resize_array(0);
  {
    // A et B sont voisins au sens des espaces distants/virtuels des facettes
    // (un processeur chez qui la facette est virtuelle envoie des donnees au
    //  processeur chez qui elle est reelle => schema_comm_inverse)
    const Schema_Comm_FT& comm = desc_facettes_.schema_comm_inverse();
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_z[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        if (PE_proprietaire == moi)
          {
            facettes_to_send_z.append_array(facette_num_owner);
            facettes_pe_dest_z.append_array(PE_destinataire);
            BtoC_send_pe_flags_z[PE_destinataire] = 1;
          }
        else
          {
            if (PE_destinataire != PE_proprietaire)
              comm.send_buffer(PE_proprietaire) << facette_num_owner << PE_destinataire;
          }
      }
    comm.echange_taille_et_messages();
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette, PE_destinataire;
            buffer >> facette >> PE_destinataire;
            if (buffer.eof())
              break;
            facettes_to_send_z.append_array(facette);
            facettes_pe_dest_z.append_array(PE_destinataire);
            BtoC_send_pe_flags_z[PE_destinataire] = 1;
          }
      }
    comm.end_comm();
  }

  // =======================================================================
  //                         Communication A -> C
  //
  // A envoie a C la liste la liste des facettes, l'expediteur et le numero
  // de l'element d'entree.
  // On remplit facettes_recues_numfacettes (numerotation temporaire)
  //            facettes_recues_numelement
  //            BtoC_recv_pe_flags (drapeaux des procs de qui on recoit des facettes).
  // Attention : si B=C on ne recoit pas de donnees de B

  // Pour chaque facette recue, numero du PE proprietaire de la facette:
  static ArrOfIntFT nouvelles_facettes_pe_proprietaire_z;
  {
    nouvelles_facettes_pe_proprietaire_z.resize_array(0);
    BtoC_recv_pe_flags_z = 0;

    // A et C sont voisins au sens du maillage eulerien
    const Schema_Comm_FT& comm = schema_comm_domaine_;
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_z[i];
        const long face_arrivee = liste_facez_arrivee_local[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        assert (PE_destinataire != moi);
        comm.send_buffer(PE_destinataire) << facette_num_owner
                                          << PE_proprietaire
                                          << face_arrivee;
      }
    comm.echange_taille_et_messages();
    // On efface facettes_recues_numfacettes et numelement.
    // Pour le cas ou il y aurait aliasing des parametres, il ne faut plus
    // utiliser liste_facettes et liste_elem_arrivee.
    facettes_recues_numfacettes.resize_array(0);
    facettes_recues_numface.resize_array(0);
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette = -1, PE_proprietaire = -1, face_arrivee = -1;
            buffer >> facette >> PE_proprietaire >> face_arrivee;
            if (buffer.eof())
              break;
            assert(facette >= 0 && PE_proprietaire >= 0 && face_arrivee >= 0);
            assert(face_arrivee < nb_face); // Ce doit etre un element reel
            // Enregistre le numero de la nouvelle facette et l'element d'arrivee.
            // Pour l'instant c'est le numero de la facette sur le proprietaire.
            facettes_recues_numfacettes.append_array(facette);
            nouvelles_facettes_pe_proprietaire_z.append_array(PE_proprietaire);
            facettes_recues_numface.append_array(face_arrivee);
            // Enregistre que l'on va recevoir des donnees du PE_proprietaire lors
            // de la derniere communication B->C
            if (PE_proprietaire != moi)
              BtoC_recv_pe_flags_z[PE_proprietaire] = 1;
          }
      }
    comm.end_comm();
  }

  // =======================================================================
  //                         Communication B -> C
  //
  // Construction du schema de communication pour envoyer les facettes reelles
  // au destinataire (B envoie les facettes a C) : il faut determiner
  // a qui j'envoie et de qui de recois.
  {
    static ArrOfIntFT send_pe_list_z;
    static ArrOfIntFT recv_pe_list_z;
    send_pe_list_z.resize_array(0);
    recv_pe_list_z.resize_array(0);
    for (i = 0; i < nbproc; i++)
      if (BtoC_send_pe_flags_z[i])
        send_pe_list_z.append_array(i);
    for (i = 0; i < nbproc; i++)
      if (BtoC_recv_pe_flags_z[i])
        recv_pe_list_z.append_array(i);
    // Envoi des facettes
    creer_facettes_virtuelles(facettes_to_send_z, facettes_pe_dest_z,
                              send_pe_list_z, recv_pe_list_z);
  }
  // Conversion du numero de la facette en numero local dans facettes_recues_numfacettes
  convertir_numero_distant_local(desc_facettes_,
                                 facette_num_owner_,
                                 facettes_recues_numfacettes,
                                 nouvelles_facettes_pe_proprietaire_z,
                                 facettes_recues_numfacettes);
  if (Comm_Group::check_enabled()) check_mesh();
}

void Maillage_FT_Disc::echanger_facettes_arete_x(const ArrOfInt& liste_facettes,
                                                 const ArrOfInt& liste_arete_arrivee,
                                                 ArrOfInt& facettes_recues_numfacettes,
                                                 ArrOfInt& facettes_recues_numarete)
{
  if (Comm_Group::check_enabled()) check_mesh();
  const long moi = me();
  long i=0;
  long indice=0;
  long nb_facettes_envoi = liste_facettes.size_array();
  const long nb_facettes_envoi_init=liste_facettes.size_array();
  // Numeros des faces euleriennes associes a chaque facette (numero local
  // de la face sur le pe qui possede cet element)
  static ArrOfIntFT liste_aretex_arrivee_local;
  ArrOfIntFT la_liste_facettes;
  la_liste_facettes.resize_array(nb_facettes_envoi);
  // Pour chaque face d'arrivee, determination du PE destination
  // et conversion du numero de la face virtuelle en numero local sur ce pe.
  static ArrOfIntFT liste_pe_dest_aretex;
  const Desc_Structure_FT& desc_facettes_const = desc_facettes();
  const Descripteur_FT& espace_distant = desc_facettes_const.espace_distant();
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& le_domaine_vf = ref_cast(Domaine_VF,domaine_dis.valeur());
  const IntVect& aretes_multiples = le_domaine_vf.aretes_multiples();
  const IntTab& arete_virt_pe_num = le_domaine_vf.arete_virt_pe_num();
  //const IntTab& aretes_multiples_pe_num = le_domaine_vf.aretes_multiples_pe_num();
  const IntTab& aretes_multiples_virt_pe_num = le_domaine_vf.aretes_multiples_virt_pe_num();

  const long nb_aretes_reelles = le_domaine_vf.domaine().nb_aretes(); // Nombre de faces reelles
  {
    liste_pe_dest_aretex.resize_array(nb_facettes_envoi);
    liste_aretex_arrivee_local.resize_array(nb_facettes_envoi);

    while (i<nb_facettes_envoi_init)
      {
        const long arete_voisine = liste_arete_arrivee[i];
        const long arete_voisine_multiple = (aretes_multiples(arete_voisine)>0);
        const long arete_mult=aretes_multiples(arete_voisine);

        if (arete_voisine_multiple)
          {
            const long fa7=liste_facettes[i];
            const long pe1=aretes_multiples_virt_pe_num(arete_voisine,0);
            const long arete_1=aretes_multiples_virt_pe_num(arete_voisine,1);
            const long pe2=aretes_multiples_virt_pe_num(arete_voisine,2);
            const long arete_2=aretes_multiples_virt_pe_num(arete_voisine,3);
            const long pe3=aretes_multiples_virt_pe_num(arete_voisine,4);
            const long arete_3=aretes_multiples_virt_pe_num(arete_voisine,5);
            const long pe4=aretes_multiples_virt_pe_num(arete_voisine,6);
            const long arete_4=aretes_multiples_virt_pe_num(arete_voisine,7);
            long fa7_reelle_pour_pe1=1;
            fa7_reelle_pour_pe1 = (sommet_PE_owner_(facettes_(fa7,0))==pe1) ? 1 : 0;
            long fa7_reelle_pour_pe2=1;
            fa7_reelle_pour_pe2 = (sommet_PE_owner_(facettes_(fa7,0))==pe2) ? 1 : 0;
            long fa7_reelle_pour_pe3=1;
            fa7_reelle_pour_pe3 = (sommet_PE_owner_(facettes_(fa7,0))==pe3) ? 1 : 0;
            long fa7_reelle_pour_pe4=1;
            fa7_reelle_pour_pe4 = (sommet_PE_owner_(facettes_(fa7,0))==pe4) ? 1 : 0;

            const long fa7_existante_pour_pe1= (pe1>=0) ? (fa7_reelle_pour_pe1 || espace_distant.contient_element(pe1, fa7)) : 0;
            const long fa7_existante_pour_pe2= (pe2>=0) ? (fa7_reelle_pour_pe2 || espace_distant.contient_element(pe2, fa7)) : 0;
            const long fa7_existante_pour_pe3= (pe3>=0) ? (fa7_reelle_pour_pe3 || espace_distant.contient_element(pe3, fa7)) : 0;
            const long fa7_existante_pour_pe4= (pe4>=0) ? (fa7_reelle_pour_pe4 || espace_distant.contient_element(pe4, fa7)) : 0;

            if (fa7_existante_pour_pe1 && pe1!=moi)
              {
                liste_pe_dest_aretex[indice] = pe1;
                liste_aretex_arrivee_local[indice] = arete_1;
                la_liste_facettes[indice]=fa7;
                indice++;
              }
            if (arete_mult>1)
              {
                i++;
                if (fa7_existante_pour_pe2 && pe2!=moi)
                  {
                    liste_pe_dest_aretex[indice] = pe2;
                    liste_aretex_arrivee_local[indice] = arete_2;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }
            if (arete_mult>2)
              {
                i++;
                if (fa7_existante_pour_pe3 && pe3!=moi)
                  {
                    liste_pe_dest_aretex[indice] = pe3;
                    liste_aretex_arrivee_local[indice] = arete_3;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }
            if (1)
              {
                i++;
                if (fa7_existante_pour_pe4 && pe4!=moi)
                  {
                    liste_pe_dest_aretex[indice] = pe4;
                    liste_aretex_arrivee_local[indice] = arete_4;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }
          }
        else
          {
            const long index = arete_voisine - nb_aretes_reelles;
            const long pe = arete_virt_pe_num(index, 0);
            const long arete_locale = arete_virt_pe_num(index, 1);
            liste_pe_dest_aretex[indice] = pe;
            liste_aretex_arrivee_local[indice] = arete_locale;
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;

          }

        i++;
      }

    liste_pe_dest_aretex.resize_array(indice);
    liste_aretex_arrivee_local.resize_array(indice);
    la_liste_facettes.resize_array(indice);
  }
  nb_facettes_envoi=indice;

  // =======================================================================
  //                         Communication A -> B
  //
  // Envoi de la liste des facettes a transmettre au processeur proprietaire de la
  // facette (A envoie les numeros de facettes a B)
  // On remplit facettes_to_send (liste des facettes que B doit envoyer a C)
  //            facettes_pe_dest (numero du processeur C a qui il faut envoyer)
  //            BtoC_send_pe_flags (drapeaux des processeurs a qui on envoie)
  // Attention : si B=C, on n'envoie pas la facette.
  static ArrOfIntFT facettes_to_send_aretex;
  static ArrOfIntFT facettes_pe_dest_aretex;
  // Lors de la communication finale entre les procs B et les procs C,
  // ces drapeaux indiquent si on recoit des donnees et si on en envoie
  // a chacun des processeurs.
  static long nbproc = Process::nproc();
  static ArrOfIntFT BtoC_send_pe_flags_x(nbproc);
  static ArrOfIntFT BtoC_recv_pe_flags_x(nbproc);
  BtoC_send_pe_flags_x = 0;
  facettes_to_send_aretex.resize_array(0);
  facettes_pe_dest_aretex.resize_array(0);
  {
    // A et B sont voisins au sens des espaces distants/virtuels des facettes
    // (un processeur chez qui la facette est virtuelle envoie des donnees au
    //  processeur chez qui elle est reelle => schema_comm_inverse)
    const Schema_Comm_FT& comm = desc_facettes_.schema_comm_inverse();
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_aretex[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        if (PE_proprietaire == moi)
          {
            facettes_to_send_aretex.append_array(facette_num_owner);
            facettes_pe_dest_aretex.append_array(PE_destinataire);
            BtoC_send_pe_flags_x[PE_destinataire] = 1;
          }
        else
          {
            if (PE_destinataire != PE_proprietaire)
              comm.send_buffer(PE_proprietaire) << facette_num_owner << PE_destinataire;
          }
      }
    comm.echange_taille_et_messages();
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette, PE_destinataire;
            buffer >> facette >> PE_destinataire;
            if (buffer.eof())
              break;
            facettes_to_send_aretex.append_array(facette);
            facettes_pe_dest_aretex.append_array(PE_destinataire);
            BtoC_send_pe_flags_x[PE_destinataire] = 1;
          }
      }
    comm.end_comm();
  }
  //Cerr << "Fin com A to B" << finl;
  // =======================================================================
  //                         Communication A -> C
  //
  // A envoie a C la liste la liste des facettes, l'expediteur et le numero
  // de l'element d'entree.
  // On remplit facettes_recues_numfacettes (numerotation temporaire)
  //            facettes_recues_numelement
  //            BtoC_recv_pe_flags (drapeaux des procs de qui on recoit des facettes).
  // Attention : si B=C on ne recoit pas de donnees de B

  // Pour chaque facette recue, numero du PE proprietaire de la facette:
  static ArrOfIntFT nouvelles_facettes_pe_proprietaire_aretex;
  {
    nouvelles_facettes_pe_proprietaire_aretex.resize_array(0);
    BtoC_recv_pe_flags_x = 0;

    // A et C sont voisins au sens du maillage eulerien
    const Schema_Comm_FT& comm = schema_comm_domaine_;
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_aretex[i];
        const long arete_arrivee = liste_aretex_arrivee_local[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        assert (PE_destinataire != moi);
        comm.send_buffer(PE_destinataire) << facette_num_owner
                                          << PE_proprietaire
                                          << arete_arrivee;
      }
    comm.echange_taille_et_messages();
    // On efface facettes_recues_numfacettes et numelement.
    // Pour le cas ou il y aurait aliasing des parametres, il ne faut plus
    // utiliser liste_facettes et liste_elem_arrivee.
    facettes_recues_numfacettes.resize_array(0);
    facettes_recues_numarete.resize_array(0);
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette = -1, PE_proprietaire = -1, arete_arrivee = -1;
            buffer >> facette >> PE_proprietaire >> arete_arrivee;
            if (buffer.eof())
              break;
            assert(facette >= 0 && PE_proprietaire >= 0 && arete_arrivee >= 0);
            // assert(arete_arrivee < nb_aretes_reelles); // Ce doit etre un element reel
            // Enregistre le numero de la nouvelle facette et l'element d'arrivee.
            // Pour l'instant c'est le numero de la facette sur le proprietaire.
            facettes_recues_numfacettes.append_array(facette);
            nouvelles_facettes_pe_proprietaire_aretex.append_array(PE_proprietaire);
            facettes_recues_numarete.append_array(arete_arrivee);
            // Enregistre que l'on va recevoir des donnees du PE_proprietaire lors
            // de la derniere communication B->C
            if (PE_proprietaire != moi)
              BtoC_recv_pe_flags_x[PE_proprietaire] = 1;
          }
      }
    comm.end_comm();
  }
  //Cerr << "Fin com A to C" << finl;
  // =====
  // ==================================================================
  //                         Communication B -> C
  //
  // Construction du schema de communication pour envoyer les facettes reelles
  // au destinataire (B envoie les facettes a C) : il faut determiner
  // a qui j'envoie et de qui de recois.
  {
    static ArrOfIntFT send_pe_list_x;
    static ArrOfIntFT recv_pe_list_x;
    send_pe_list_x.resize_array(0);
    recv_pe_list_x.resize_array(0);
    for (i = 0; i < nbproc; i++)
      if (BtoC_send_pe_flags_x[i])
        send_pe_list_x.append_array(i);
    for (i = 0; i < nbproc; i++)
      if (BtoC_recv_pe_flags_x[i])
        recv_pe_list_x.append_array(i);
    // Envoi des facettes
    creer_facettes_virtuelles(facettes_to_send_aretex, facettes_pe_dest_aretex,
                              send_pe_list_x, recv_pe_list_x);
  }
  //Cerr << "Fin com B to C" << finl;
  // Conversion du numero de la facette en numero local dans facettes_recues_numfacettes
  convertir_numero_distant_local(desc_facettes_,
                                 facette_num_owner_,
                                 facettes_recues_numfacettes,
                                 nouvelles_facettes_pe_proprietaire_aretex,
                                 facettes_recues_numfacettes);
  if (Comm_Group::check_enabled()) check_mesh();
}
void Maillage_FT_Disc::echanger_facettes_arete_y(const ArrOfInt& liste_facettes,
                                                 const ArrOfInt& liste_arete_arrivee,
                                                 ArrOfInt& facettes_recues_numfacettes,
                                                 ArrOfInt& facettes_recues_numarete)
{
  if (Comm_Group::check_enabled()) check_mesh();
  const long moi = me();
  long i=0;
  long indice=0;
  long nb_facettes_envoi = liste_facettes.size_array();
  const long nb_facettes_envoi_init=liste_facettes.size_array();

  // Numeros des faces euleriennes associes a chaque facette (numero local
  // de la face sur le pe qui possede cet element)
  static ArrOfIntFT liste_aretey_arrivee_local;
  ArrOfIntFT la_liste_facettes;
  la_liste_facettes.resize_array(nb_facettes_envoi);
  // Pour chaque face d'arrivee, determination du PE destination
  // et conversion du numero de la face virtuelle en numero local sur ce pe.
  static ArrOfIntFT liste_pe_dest_aretey;
  const Desc_Structure_FT& desc_facettes_const = desc_facettes();
  const Descripteur_FT& espace_distant = desc_facettes_const.espace_distant();
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& le_domaine_vf = ref_cast(Domaine_VF,domaine_dis.valeur());
  const IntVect& aretes_multiples = le_domaine_vf.aretes_multiples();
  const IntTab& arete_virt_pe_num = le_domaine_vf.arete_virt_pe_num();
  //const IntTab& aretes_multiples_pe_num = le_domaine_vf.aretes_multiples_pe_num();
  const IntTab& aretes_multiples_virt_pe_num = le_domaine_vf.aretes_multiples_virt_pe_num();
  const long nb_aretes_reelles = le_domaine_vf.domaine().nb_aretes(); // Nombre de faces reelles
  {
    liste_pe_dest_aretey.resize_array(nb_facettes_envoi);
    liste_aretey_arrivee_local.resize_array(nb_facettes_envoi);

    while (i<nb_facettes_envoi_init)
      {
        const long arete_voisine = liste_arete_arrivee[i];
        const long arete_voisine_multiple = (aretes_multiples(arete_voisine)>0);
        const long arete_mult=aretes_multiples(arete_voisine);

        if (arete_voisine_multiple)
          {
            const long fa7=liste_facettes[i];
            const long pe1=aretes_multiples_virt_pe_num(arete_voisine,0);
            const long arete_1=aretes_multiples_virt_pe_num(arete_voisine,1);
            const long pe2=aretes_multiples_virt_pe_num(arete_voisine,2);
            const long arete_2=aretes_multiples_virt_pe_num(arete_voisine,3);
            const long pe3=aretes_multiples_virt_pe_num(arete_voisine,4);
            const long arete_3=aretes_multiples_virt_pe_num(arete_voisine,5);
            const long pe4=aretes_multiples_virt_pe_num(arete_voisine,6);
            const long arete_4=aretes_multiples_virt_pe_num(arete_voisine,7);
            long fa7_reelle_pour_pe1=1;
            fa7_reelle_pour_pe1 = (sommet_PE_owner_(facettes_(fa7,0))==pe1) ? 1 : 0;
            long fa7_reelle_pour_pe2=1;
            fa7_reelle_pour_pe2 = (sommet_PE_owner_(facettes_(fa7,0))==pe2) ? 1 : 0;
            long fa7_reelle_pour_pe3=1;
            fa7_reelle_pour_pe3 = (sommet_PE_owner_(facettes_(fa7,0))==pe3) ? 1 : 0;
            long fa7_reelle_pour_pe4=1;
            fa7_reelle_pour_pe4 = (sommet_PE_owner_(facettes_(fa7,0))==pe4) ? 1 : 0;

            const long fa7_existante_pour_pe1= (pe1>=0) ? fa7_reelle_pour_pe1 || espace_distant.contient_element(pe1, fa7) : 0;
            const long fa7_existante_pour_pe2= (pe2>=0) ? fa7_reelle_pour_pe2 || espace_distant.contient_element(pe2, fa7) : 0;
            const long fa7_existante_pour_pe3= (pe3>=0) ? fa7_reelle_pour_pe3 || espace_distant.contient_element(pe3, fa7) : 0;
            const long fa7_existante_pour_pe4= (pe4>=0) ? fa7_reelle_pour_pe4 || espace_distant.contient_element(pe4, fa7) : 0;

            if (fa7_existante_pour_pe1 && pe1!=moi)
              {
                liste_pe_dest_aretey[indice] = pe1;
                liste_aretey_arrivee_local[indice] = arete_1;
                la_liste_facettes[indice]=fa7;
                indice++;
              }
            if (arete_mult>1)
              {
                i++;
                if (fa7_existante_pour_pe2 && pe2!=moi)
                  {
                    liste_pe_dest_aretey[indice] = pe2;
                    liste_aretey_arrivee_local[indice] = arete_2;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }
            if (arete_mult>2)
              {
                i++;
                if (fa7_existante_pour_pe3 && pe3!=moi)
                  {
                    liste_pe_dest_aretey[indice] = pe3;
                    liste_aretey_arrivee_local[indice] = arete_3;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }

            if (1) // (arete_voisine>nb_aretes_reelles)
              {
                i++;
                if (fa7_existante_pour_pe4 && pe4!=moi)
                  {
                    liste_pe_dest_aretey[indice] = pe4;
                    liste_aretey_arrivee_local[indice] = arete_4;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }
          }
        else
          {
            const long index = arete_voisine - nb_aretes_reelles;
            const long pe = arete_virt_pe_num(index, 0);
            const long arete_locale = arete_virt_pe_num(index, 1);
            liste_pe_dest_aretey[indice] = pe;
            liste_aretey_arrivee_local[indice] = arete_locale;
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;
          }

        i++;
      }

    liste_pe_dest_aretey.resize_array(indice);
    liste_aretey_arrivee_local.resize_array(indice);
    la_liste_facettes.resize_array(indice);
  }
  nb_facettes_envoi=indice;

  // =======================================================================
  //                         Communication A -> B
  //
  // Envoi de la liste des facettes a transmettre au processeur proprietaire de la
  // facette (A envoie les numeros de facettes a B)
  // On remplit facettes_to_send (liste des facettes que B doit envoyer a C)
  //            facettes_pe_dest (numero du processeur C a qui il faut envoyer)
  //            BtoC_send_pe_flags (drapeaux des processeurs a qui on envoie)
  // Attention : si B=C, on n'envoie pas la facette.
  static ArrOfIntFT facettes_to_send_aretey;
  static ArrOfIntFT facettes_pe_dest_aretey;
  // Lors de la communication finale entre les procs B et les procs C,
  // ces drapeaux indiquent si on recoit des donnees et si on en envoie
  // a chacun des processeurs.
  static long nbproc = Process::nproc();
  static ArrOfIntFT BtoC_send_pe_flags_aretey(nbproc);
  static ArrOfIntFT BtoC_recv_pe_flags_aretey(nbproc);
  BtoC_send_pe_flags_aretey = 0;
  facettes_to_send_aretey.resize_array(0);
  facettes_pe_dest_aretey.resize_array(0);
  {
    // A et B sont voisins au sens des espaces distants/virtuels des facettes
    // (un processeur chez qui la facette est virtuelle envoie des donnees au
    //  processeur chez qui elle est reelle => schema_comm_inverse)
    const Schema_Comm_FT& comm = desc_facettes_.schema_comm_inverse();
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_aretey[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        if (PE_proprietaire == moi)
          {
            facettes_to_send_aretey.append_array(facette_num_owner);
            facettes_pe_dest_aretey.append_array(PE_destinataire);
            BtoC_send_pe_flags_aretey[PE_destinataire] = 1;
          }
        else
          {
            if (PE_destinataire != PE_proprietaire)
              comm.send_buffer(PE_proprietaire) << facette_num_owner << PE_destinataire;
          }
      }
    comm.echange_taille_et_messages();
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette, PE_destinataire;
            buffer >> facette >> PE_destinataire;
            if (buffer.eof())
              break;
            facettes_to_send_aretey.append_array(facette);
            facettes_pe_dest_aretey.append_array(PE_destinataire);
            BtoC_send_pe_flags_aretey[PE_destinataire] = 1;
          }
      }
    comm.end_comm();
  }
  //Cerr << "Fin com A to B" << finl;
  // =======================================================================
  //                         Communication A -> C
  //
  // A envoie a C la liste la liste des facettes, l'expediteur et le numero
  // de l'element d'entree.
  // On remplit facettes_recues_numfacettes (numerotation temporaire)
  //            facettes_recues_numelement
  //            BtoC_recv_pe_flags (drapeaux des procs de qui on recoit des facettes).
  // Attention : si B=C on ne recoit pas de donnees de B

  // Pour chaque facette recue, numero du PE proprietaire de la facette:
  static ArrOfIntFT nouvelles_facettes_pe_proprietaire_y;
  {
    nouvelles_facettes_pe_proprietaire_y.resize_array(0);
    BtoC_recv_pe_flags_aretey = 0;

    // A et C sont voisins au sens du maillage eulerien
    const Schema_Comm_FT& comm = schema_comm_domaine_;
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_aretey[i];
        const long arete_arrivee = liste_aretey_arrivee_local[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        assert (PE_destinataire != moi);
        comm.send_buffer(PE_destinataire) << facette_num_owner
                                          << PE_proprietaire
                                          << arete_arrivee;
      }
    comm.echange_taille_et_messages();
    // On efface facettes_recues_numfacettes et numelement.
    // Pour le cas ou il y aurait aliasing des parametres, il ne faut plus
    // utiliser liste_facettes et liste_elem_arrivee.
    facettes_recues_numfacettes.resize_array(0);
    facettes_recues_numarete.resize_array(0);
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette = -1, PE_proprietaire = -1, arete_arrivee = -1;
            buffer >> facette >> PE_proprietaire >> arete_arrivee;
            if (buffer.eof())
              break;
            assert(facette >= 0 && PE_proprietaire >= 0 && arete_arrivee >= 0);
            // assert(arete_arrivee < nb_aretes_reelles); // Ce doit etre un element reel
            // Enregistre le numero de la nouvelle facette et l'element d'arrivee.
            // Pour l'instant c'est le numero de la facette sur le proprietaire.
            facettes_recues_numfacettes.append_array(facette);
            nouvelles_facettes_pe_proprietaire_y.append_array(PE_proprietaire);
            facettes_recues_numarete.append_array(arete_arrivee);
            // Enregistre que l'on va recevoir des donnees du PE_proprietaire lors
            // de la derniere communication B->C
            if (PE_proprietaire != moi)
              BtoC_recv_pe_flags_aretey[PE_proprietaire] = 1;
          }
      }
    comm.end_comm();
  }
  //Cerr << "Fin com A to C" << finl;
  // =====
  // ==================================================================
  //                         Communication B -> C
  //
  // Construction du schema de communication pour envoyer les facettes reelles
  // au destinataire (B envoie les facettes a C) : il faut determiner
  // a qui j'envoie et de qui de recois.
  {
    static ArrOfIntFT send_pe_list_y;
    static ArrOfIntFT recv_pe_list_y;
    send_pe_list_y.resize_array(0);
    recv_pe_list_y.resize_array(0);
    for (i = 0; i < nbproc; i++)
      if (BtoC_send_pe_flags_aretey[i])
        send_pe_list_y.append_array(i);
    for (i = 0; i < nbproc; i++)
      if (BtoC_recv_pe_flags_aretey[i])
        recv_pe_list_y.append_array(i);
    // Envoi des facettes
    creer_facettes_virtuelles(facettes_to_send_aretey, facettes_pe_dest_aretey,
                              send_pe_list_y, recv_pe_list_y);
  }
  //Cerr << "Fin com B to C" << finl;
  // Conversion du numero de la facette en numero local dans facettes_recues_numfacettes
  convertir_numero_distant_local(desc_facettes_,
                                 facette_num_owner_,
                                 facettes_recues_numfacettes,
                                 nouvelles_facettes_pe_proprietaire_y,
                                 facettes_recues_numfacettes);
  if (Comm_Group::check_enabled()) check_mesh();
}
void Maillage_FT_Disc::echanger_facettes_arete_z(const ArrOfInt& liste_facettes,
                                                 const ArrOfInt& liste_arete_arrivee,
                                                 ArrOfInt& facettes_recues_numfacettes,
                                                 ArrOfInt& facettes_recues_numarete)
{
  if (Comm_Group::check_enabled()) check_mesh();
  const long moi = me();
  long i=0;
  long indice=0;
  long nb_facettes_envoi = liste_facettes.size_array();
  const long nb_facettes_envoi_init=liste_facettes.size_array();
  // Numeros des faces euleriennes associes a chaque facette (numero local
  // de la face sur le pe qui possede cet element)
  static ArrOfIntFT liste_aretez_arrivee_local;
  ArrOfIntFT la_liste_facettes;
  la_liste_facettes.resize_array(nb_facettes_envoi);
  // Pour chaque face d'arrivee, determination du PE destination
  // et conversion du numero de la face virtuelle en numero local sur ce pe.
  static ArrOfIntFT liste_pe_dest_aretez;
  const Desc_Structure_FT& desc_facettes_const = desc_facettes();
  const Descripteur_FT& espace_distant = desc_facettes_const.espace_distant();
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& le_domaine_vf = ref_cast(Domaine_VF,domaine_dis.valeur());
  const IntVect& aretes_multiples = le_domaine_vf.aretes_multiples();
  const IntTab& arete_virt_pe_num = le_domaine_vf.arete_virt_pe_num();
  //const IntTab& aretes_multiples_pe_num = le_domaine_vf.aretes_multiples_pe_num();
  const IntTab& aretes_multiples_virt_pe_num = le_domaine_vf.aretes_multiples_virt_pe_num();
  const long nb_aretes_reelles = le_domaine_vf.domaine().nb_aretes(); // Nombre de faces reelles
  {
    liste_pe_dest_aretez.resize_array(nb_facettes_envoi);
    liste_aretez_arrivee_local.resize_array(nb_facettes_envoi);

    while (i<nb_facettes_envoi_init)
      {
        const long arete_voisine = liste_arete_arrivee[i];
        const long arete_voisine_multiple = (aretes_multiples(arete_voisine)>0);
        const long arete_mult=aretes_multiples(arete_voisine);
        if (arete_voisine_multiple)
          {
            const long fa7=liste_facettes[i];
            const long pe1=aretes_multiples_virt_pe_num(arete_voisine,0);
            const long arete_1=aretes_multiples_virt_pe_num(arete_voisine,1);
            const long pe2=aretes_multiples_virt_pe_num(arete_voisine,2);
            const long arete_2=aretes_multiples_virt_pe_num(arete_voisine,3);
            const long pe3=aretes_multiples_virt_pe_num(arete_voisine,4);
            const long arete_3=aretes_multiples_virt_pe_num(arete_voisine,5);
            const long pe4=aretes_multiples_virt_pe_num(arete_voisine,6);
            const long arete_4=aretes_multiples_virt_pe_num(arete_voisine,7);

            long fa7_reelle_pour_pe1=1;
            fa7_reelle_pour_pe1 = (sommet_PE_owner_(facettes_(fa7,0))==pe1) ? 1 : 0;
            long fa7_reelle_pour_pe2=1;
            fa7_reelle_pour_pe2 = (sommet_PE_owner_(facettes_(fa7,0))==pe2) ? 1 : 0;
            long fa7_reelle_pour_pe3=1;
            fa7_reelle_pour_pe3 = (sommet_PE_owner_(facettes_(fa7,0))==pe3) ? 1 : 0;
            long fa7_reelle_pour_pe4=1;
            fa7_reelle_pour_pe4 = (sommet_PE_owner_(facettes_(fa7,0))==pe4) ? 1 : 0;

            const long fa7_existante_pour_pe1= (pe1>=0) ? fa7_reelle_pour_pe1 || espace_distant.contient_element(pe1, fa7) : 0;
            const long fa7_existante_pour_pe2= (pe2>=0) ? fa7_reelle_pour_pe2 || espace_distant.contient_element(pe2, fa7) : 0;
            const long fa7_existante_pour_pe3= (pe3>=0) ? fa7_reelle_pour_pe3 || espace_distant.contient_element(pe3, fa7) : 0;
            const long fa7_existante_pour_pe4= (pe4>=0) ? fa7_reelle_pour_pe4 || espace_distant.contient_element(pe4, fa7) : 0;

            /*
            Cerr << "arete_mult " << arete_mult <<
                 "\tpe_1 " << pe1 << " a_1 " << arete_1 << " " << fa7_reelle_pour_pe1 << " " << fa7_existante_pour_pe1 ;
            Cerr << "\tpe2 " << pe2 << " a_2 " << arete_2 << " " << fa7_reelle_pour_pe2 << " " << fa7_existante_pour_pe2 ;
            Cerr << "\tpe3 " << pe3 << " a_3 " << arete_3 << " " << fa7_reelle_pour_pe3 << " " << fa7_existante_pour_pe3;
            Cerr << "\tpe4 " << pe4 << " a_4 " << arete_4 << " " << fa7_reelle_pour_pe4 << " " << fa7_existante_pour_pe4 <<
                 "\tcg " << la_zone_vf.xa(arete_voisine,0) << " " << la_zone_vf.xa(arete_voisine,1) << " " << la_zone_vf.xa(arete_voisine,2) <<
                 "\tfa7 " << cg_fa7_(fa7,0) << " " << cg_fa7_(fa7,1) << " " <<cg_fa7_(fa7,2) << finl;

            */
            if (fa7_existante_pour_pe1 && pe1!=moi)
              {
                liste_pe_dest_aretez[indice] = pe1;
                liste_aretez_arrivee_local[indice] = arete_1;
                la_liste_facettes[indice]=fa7;
                indice++;
              }
            if (arete_mult>1)
              {
                i++;
                if (fa7_existante_pour_pe2 && pe2!=moi)
                  {
                    liste_pe_dest_aretez[indice] = pe2;
                    liste_aretez_arrivee_local[indice] = arete_2;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }
            if (arete_mult>2)
              {
                i++;
                if (fa7_existante_pour_pe3 && pe3!=moi)
                  {
                    liste_pe_dest_aretez[indice] = pe3;
                    liste_aretez_arrivee_local[indice] = arete_3;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }
            if(1) // (arete_voisine>nb_aretes_reelles)
              {
                i++;
                if (fa7_existante_pour_pe4 && pe4!=moi)
                  {
                    liste_pe_dest_aretez[indice] = pe4;
                    liste_aretez_arrivee_local[indice] = arete_4;
                    la_liste_facettes[indice]=fa7;
                    indice++;
                  }
              }
          }
        else
          {
            const long index = arete_voisine - nb_aretes_reelles;
            const long pe = arete_virt_pe_num(index, 0);
            const long arete_locale = arete_virt_pe_num(index, 1);
            liste_pe_dest_aretez[indice] = pe;
            liste_aretez_arrivee_local[indice] = arete_locale;
            la_liste_facettes[indice]=liste_facettes[i];
            indice++;

          }

        i++;
      }

    liste_pe_dest_aretez.resize_array(indice);
    liste_aretez_arrivee_local.resize_array(indice);
    la_liste_facettes.resize_array(indice);
  }
  nb_facettes_envoi=indice;
  /*Cerr << "APRES TRI" << finl;
  for (i =0; i<liste_aretez_arrivee_local.size_array(); i++)
    {
      long arete=liste_aretez_arrivee_local(i);
      Cerr << "arete " << arete << " pe_dist " << liste_pe_dest_aretez(i) << finl;
    }
    */
  // =======================================================================
  //                         Communication A -> B
  //
  // Envoi de la liste des facettes a transmettre au processeur proprietaire de la
  // facette (A envoie les numeros de facettes a B)
  // On remplit facettes_to_send (liste des facettes que B doit envoyer a C)
  //            facettes_pe_dest (numero du processeur C a qui il faut envoyer)
  //            BtoC_send_pe_flags (drapeaux des processeurs a qui on envoie)
  // Attention : si B=C, on n'envoie pas la facette.
  static ArrOfIntFT facettes_to_send_aretez;
  static ArrOfIntFT facettes_pe_dest_aretez;
  // Lors de la communication finale entre les procs B et les procs C,
  // ces drapeaux indiquent si on recoit des donnees et si on en envoie
  // a chacun des processeurs.
  static long nbproc = Process::nproc();
  static ArrOfIntFT BtoC_send_pe_flags_aretez(nbproc);
  static ArrOfIntFT BtoC_recv_pe_flags_aretez(nbproc);
  BtoC_send_pe_flags_aretez = 0;
  facettes_to_send_aretez.resize_array(0);
  facettes_pe_dest_aretez.resize_array(0);
  {
    // A et B sont voisins au sens des espaces distants/virtuels des facettes
    // (un processeur chez qui la facette est virtuelle envoie des donnees au
    //  processeur chez qui elle est reelle => schema_comm_inverse)
    const Schema_Comm_FT& comm = desc_facettes_.schema_comm_inverse();
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_aretez[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        if (PE_proprietaire == moi)
          {
            facettes_to_send_aretez.append_array(facette_num_owner);
            facettes_pe_dest_aretez.append_array(PE_destinataire);
            BtoC_send_pe_flags_aretez[PE_destinataire] = 1;
          }
        else
          {
            if (PE_destinataire != PE_proprietaire)
              comm.send_buffer(PE_proprietaire) << facette_num_owner << PE_destinataire;
          }
      }
    comm.echange_taille_et_messages();
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette, PE_destinataire;
            buffer >> facette >> PE_destinataire;
            if (buffer.eof())
              break;
            facettes_to_send_aretez.append_array(facette);
            facettes_pe_dest_aretez.append_array(PE_destinataire);
            BtoC_send_pe_flags_aretez[PE_destinataire] = 1;
          }
      }
    comm.end_comm();
  }
  //Cerr << "Fin com A to B" << finl;
  // =======================================================================
  //                         Communication A -> C
  //
  // A envoie a C la liste la liste des facettes, l'expediteur et le numero
  // de l'element d'entree.
  // On remplit facettes_recues_numfacettes (numerotation temporaire)
  //            facettes_recues_numelement
  //            BtoC_recv_pe_flags (drapeaux des procs de qui on recoit des facettes).
  // Attention : si B=C on ne recoit pas de donnees de B

  // Pour chaque facette recue, numero du PE proprietaire de la facette:
  static ArrOfIntFT nouvelles_facettes_pe_proprietaire_aretez;
  {
    nouvelles_facettes_pe_proprietaire_aretez.resize_array(0);
    BtoC_recv_pe_flags_aretez = 0;

    // A et C sont voisins au sens du maillage eulerien
    const Schema_Comm_FT& comm = schema_comm_domaine_;
    comm.begin_comm();
    for (i = 0; i < nb_facettes_envoi; i++)
      {
        const long facette = la_liste_facettes[i];
        const long PE_destinataire = liste_pe_dest_aretez[i];
        const long arete_arrivee = liste_aretez_arrivee_local[i];
        const long premier_sommet = facettes_(facette, 0);
        const long PE_proprietaire = sommet_PE_owner_[premier_sommet];
        const long facette_num_owner = facette_num_owner_[facette];
        assert (PE_destinataire != moi);
        comm.send_buffer(PE_destinataire) << facette_num_owner
                                          << PE_proprietaire
                                          << arete_arrivee;
      }
    comm.echange_taille_et_messages();
    // On efface facettes_recues_numfacettes et numelement.
    // Pour le cas ou il y aurait aliasing des parametres, il ne faut plus
    // utiliser liste_facettes et liste_elem_arrivee.
    facettes_recues_numfacettes.resize_array(0);
    facettes_recues_numarete.resize_array(0);
    const ArrOfInt& recv_pe_list = comm.get_recv_pe_list();
    const long nb_recv_pe = recv_pe_list.size_array();
    for (i = 0; i < nb_recv_pe; i++)
      {
        const long pe_source = recv_pe_list[i];
        Entree& buffer = comm.recv_buffer(pe_source);
        while (1)
          {
            long facette = -1, PE_proprietaire = -1, arete_arrivee = -1;
            buffer >> facette >> PE_proprietaire >> arete_arrivee;
            if (buffer.eof())
              break;
            assert(facette >= 0 && PE_proprietaire >= 0 && arete_arrivee >= 0);
            // assert(arete_arrivee < nb_aretes_reelles); // Ce doit etre un element reel
            // Enregistre le numero de la nouvelle facette et l'element d'arrivee.
            // Pour l'instant c'est le numero de la facette sur le proprietaire.
            facettes_recues_numfacettes.append_array(facette);
            nouvelles_facettes_pe_proprietaire_aretez.append_array(PE_proprietaire);
            facettes_recues_numarete.append_array(arete_arrivee);
            // Enregistre que l'on va recevoir des donnees du PE_proprietaire lors
            // de la derniere communication B->C
            if (PE_proprietaire != moi)
              BtoC_recv_pe_flags_aretez[PE_proprietaire] = 1;
          }
      }
    comm.end_comm();
  }
  //Cerr << "Fin com A to C" << finl;
  // =====
  // ==================================================================
  //                         Communication B -> C
  //
  // Construction du schema de communication pour envoyer les facettes reelles
  // au destinataire (B envoie les facettes a C) : il faut determiner
  // a qui j'envoie et de qui de recois.
  {
    static ArrOfIntFT send_pe_list_aretez;
    static ArrOfIntFT recv_pe_list_aretez;
    send_pe_list_aretez.resize_array(0);
    recv_pe_list_aretez.resize_array(0);
    for (i = 0; i < nbproc; i++)
      if (BtoC_send_pe_flags_aretez[i])
        send_pe_list_aretez.append_array(i);
    for (i = 0; i < nbproc; i++)
      if (BtoC_recv_pe_flags_aretez[i])
        recv_pe_list_aretez.append_array(i);
    // Envoi des facettes
    creer_facettes_virtuelles(facettes_to_send_aretez, facettes_pe_dest_aretez,
                              send_pe_list_aretez, recv_pe_list_aretez);
  }
  //Cerr << "Fin com B to C" << finl;
  // Conversion du numero de la facette en numero local dans facettes_recues_numfacettes
  convertir_numero_distant_local(desc_facettes_,
                                 facette_num_owner_,
                                 facettes_recues_numfacettes,
                                 nouvelles_facettes_pe_proprietaire_aretez,
                                 facettes_recues_numfacettes);
  if (Comm_Group::check_enabled()) check_mesh();
}
/*! @brief Conversion des couples (numeros_distants, pe) en numeros_locaux.
 *
 * Parametres:
 *   descripteur       = soit desc_sommets_, soit desc_facettes_, ou un autre
 *   element_num_owner = soit sommets_num_owner_, soit facettes_num_owner_, ou equivalent
 *   numeros_distants  = une liste de numeros de facettes ou sommets qui appartiennent au pe_distant
 *                      et qui sont virtuels chez moi. Le numero est le numero
 *                      de la facette ou du sommet sur le pe_distant
 *   pe_distant        = pour chaque numero de facette ou sommet, le pe proprietaire
 *
 * @param (numeros_locaux) tableau rempli par cette methode. Le contenu initial est efface. La taille est identique a celle de numeros_distants. On y met le numero local des elements (ce sont des numeros d'elements virtuels). Attention: numeros_locaux et numeros_distants peuvent referencer le meme tableau
 */
void Maillage_FT_Disc::convertir_numero_distant_local(const Desc_Structure_FT& descripteur,
                                                      const ArrOfInt& element_num_owner,
                                                      const ArrOfInt& numeros_distants,
                                                      const ArrOfInt& pe_distant,
                                                      ArrOfInt& numeros_locaux) const
{
  const Descripteur_FT& espace_virtuel = descripteur.espace_virtuel();
  const long nb_numeros = numeros_distants.size_array();
  numeros_locaux.resize_array(nb_numeros);
  const long moi = me();

  for (long index = 0; index < nb_numeros; index++)
    {
      const long numero_distant = numeros_distants[index];
      const long pe = pe_distant[index];
      long numero_local = -1;

      if (pe == moi)
        {
          numero_local = numero_distant;
        }
      else
        {
          // Recherche parmi les sommets virtuels du pe:
          const ArrOfInt& elements_virtuels = espace_virtuel.elements(pe);
          const long n = elements_virtuels.size_array();
          // Elem est un numero de facette ou de sommet
          long i;
          for (i = 0; i < n; i++)
            {
              const long elem = elements_virtuels[i];
              const long numero_distant_de_elem = element_num_owner[elem];
              if (numero_distant == numero_distant_de_elem)
                {
                  numero_local = elem;
                  break;
                }
            }
          assert(numero_local >= 0);
        }
      numeros_locaux[index] = numero_local;
    }
}
void Maillage_FT_Disc::convertir_numero_distant_local(const Desc_Structure_FT& descripteur,
                                                      const ArrOfInt& element_num_owner,
                                                      const long numero_distant,
                                                      const long pe_distant,
                                                      long& numero_local) const
{
  const Descripteur_FT& espace_virtuel = descripteur.espace_virtuel();
  const long moi = me();

  numero_local = -1;

  if (pe_distant == moi)
    {
      numero_local = numero_distant;
    }
  else
    {
      // Recherche parmi les sommets virtuels du pe:
      const ArrOfInt& elements_virtuels = espace_virtuel.elements(pe_distant);
      const long n = elements_virtuels.size_array();
      // Elem est un numero de facette ou de sommet
      long i;
      for (i = 0; i < n; i++)
        {
          const long elem = elements_virtuels[i];
          const long numero_distant_de_elem = element_num_owner[elem];
          if (numero_distant == numero_distant_de_elem)
            {
              numero_local = elem;
              break;
            }
        }
      //assert(numero_local >= 0);
    }
}

// Description :
//  Methode outil utilisee dans "deplacer_sommets".
//  (pour decouper le code en petits bouts...).
//  On modifie x, y, z, element, face_bord et eventuellement sommets_envoyes...
// x,y,z:    position initiale en entree, position finale en sortie
//           (le sommet est deplace en direction de la destination jusqu'a
//            ce qu'il sorte de l'"element")
// x1,y1,z1: destination
// element : en entree, numero de l'element ou se trouve le sommet,
//           en sortie, numero de l'element ou arrive le sommet,
// face_bord: idem pour la face de bord si le sommet est sur un bord,
//            -1 sinon.
// sommets_envoyes, ... : donnees pour les echanges entre processeurs si
//           le sommet traverse un joint.
// valeur de retour : 0 si x1,y1,z1 est dans l'"element" initial (on a fini le deplacement)
//                    1 sinon (la destination n'a pas encore ete atteinte).

//RQ : Fonction scindee en 2 :
// - deplacement du point (deplacer_un_point)
// - deplacement du sommet (deplacer_un_sommet) : appel a deplacer_un_point + gestion des transferts entre processeurs
//fonctions statiques (pour ne pas utiliser des membres propres au maillage
long Maillage_FT_Disc::deplacer_un_point(double& x, double& y, double& z,
                                         double x1, double y1, double z1,
                                         long& element,
                                         long& face_bord,
                                         const Parcours_interface& parcours,
                                         const Domaine_VF& domaine_vf,
                                         const IntTab& face_voisins,
                                         long skip_facettes)
{
  // Valeur de retour de la fonction :
  long continuer = -2;
  // Rappel : face_bord < -1  <=>  le sommet est sur une ligne de contact
  // La face de bord ou passe le sommet (si ligne contact) :
  long face_suivante = -2;
  // L'element ou passe le sommet :
  long element_suivant = -2;

  if (face_bord < 0)
    {
      // Le noeud n'est pas une ligne de contact
      // Si on ne trouve pas d'intersection, on recupere pos_intersection=1.
      double pos_intersection = 1.;
      const long face_sortie =
        parcours.calculer_face_sortie_element(domaine_vf, element,
                                              x, y, z,
                                              x1, y1, z1,
                                              pos_intersection);
      // Mettre a jour la position x,y,z
      x = x * (1. - pos_intersection) + x1 * pos_intersection;
      y = y * (1. - pos_intersection) + y1 * pos_intersection;
      z = z * (1. - pos_intersection) + z1 * pos_intersection;
      if (face_sortie >= 0)
        {
          // Le sommet sort de l'element, quel est l'element voisin ?
          const long elem0 = face_voisins(face_sortie, 0);
          const long elem1 = face_voisins(face_sortie, 1);
          element_suivant = elem0 + elem1 - element;
          if (element_suivant < 0)
            {
              // Le sommet touche un bord => il devient ligne de contact
              element_suivant = element;
              face_suivante = face_sortie;
              continuer = 1;
            }
          else
            {
              // Le sommet passe dans un autre element
              face_suivante = -1;
              continuer = 1;
            }
        }
      else
        {
          // Le sommet ne change pas d'element
          element_suivant = element;
          face_suivante = -1;
          continuer = 0;
        }
    }
  else
    {
      //On ne modifie pas face_bord et elem dans le cas des marqueurss (skip_facettes)
      if (!skip_facettes)
        {
          // Le sommet est sur un bord => ligne de contact
          // Le sommet se deplace sur la face de bord. Par ou sort-il ?
          const long nouvelle_face_bord =
            parcours.calculer_sortie_face_bord(face_bord, element,
                                               x, y, z,
                                               x1, y1, z1,
                                               x, y, z); // Nouvelle position
          if (nouvelle_face_bord >= 0)
            {
              // Le sommet change de face de bord :
              // Numero de la face suivante :
              face_suivante = nouvelle_face_bord;
              // Numero de l'element adjacent a cette face
              const long elem0 = face_voisins(face_suivante, 0);
              const long elem1 = face_voisins(face_suivante, 1);
              element_suivant = elem0 + elem1 + 1;
              assert((elem0 == -1 || elem1 == -1)
                     && (elem0 >= 0 || elem1 >= 0));
              continuer = 1;
            }
          else
            {
              // Le sommet reste sur la meme face de bord
              face_suivante = face_bord;
              element_suivant = element;
              continuer = 0;
            }
        } //fin !skip_facettes
      else
        {
          element_suivant = element;
          face_suivante = face_bord;
          continuer = 0;
        }
    }

  element = element_suivant;
  face_bord = face_suivante;



  return continuer;
}

long Maillage_FT_Disc::deplacer_un_sommet(double& x, double& y, double& z,
                                          double x1, double y1, double z1,
                                          long& element,
                                          long& face_bord,
                                          const long num_sommet,
                                          const Parcours_interface& parcours,
                                          const Domaine_VF& domaine_vf,
                                          const IntTab& face_voisins,
                                          ArrOfInt& sommets_envoyes,
                                          ArrOfInt& element_virtuel_arrivee,
                                          ArrOfInt& face_virtuelle_arrivee,
                                          DoubleTab& deplacement_restant,
                                          long skip_facettes)

{
  // Valeur de retour de la fonction :
  long continuer = -2;
  // Rappel : face_bord < -1  <=>  le sommet est sur une ligne de contact
  // La face de bord ou passe le sommet (si ligne contact) :
  long face_suivante = -2;
  // L'element ou passe le sommet :
  long element_suivant = -2;

  // -------------------------------------------------------------
  // PREMIERE ETAPE : calcul de "face_suivante", "element_suivant" et "continuer" :
  // -------------------------------------------------------------

  element_suivant = element;
  face_suivante = face_bord;
  continuer = Maillage_FT_Disc::deplacer_un_point(x,y,z, x1,y1,z1,
                                                  element_suivant, face_suivante,
                                                  parcours, domaine_vf, face_voisins,
                                                  skip_facettes);

  // On n'a rien oublie ?
  assert(continuer > -2 && face_suivante > -2 && element_suivant > -2);

  // -------------------------------------------------------------
  // DEUXIEME ETAPE : traitement des traversees de joints entre processeurs
  // -------------------------------------------------------------
  if (continuer)
    {
      const long nb_elem_reels = domaine_vf.nb_elem();
      if (element_suivant >= nb_elem_reels)
        {
          if (num_sommet>=0)
            {
              //cas general
              // L'element voisin est virtuel, on traverse un joint.
              // On empile les valeurs a transmettre.
              sommets_envoyes.append_array(num_sommet);
              element_virtuel_arrivee.append_array(element_suivant);
              face_virtuelle_arrivee.append_array(face_suivante);
              const long n = deplacement_restant.dimension(0);
              deplacement_restant.resize(n+1, 3);
              deplacement_restant(n,0) = x1 - x;
              deplacement_restant(n,1) = y1 - y;
              deplacement_restant(n,2) = z1 - z;
              // Le sommet reste pour l'instant dans l'element reel courant.
              // Il sera mis a jour dans "echanger_sommets_PE".
              // On arrete les iterations pour ce sommet, la suite sera
              // traitee par un autre processeur:
              continuer = 0;
            }
          else
            {
              //cas ou on deplace un point, pas un sommet
              continuer = -1;
            }
        }
      else
        {
          // Mise a jour de element et face_bord
          element = element_suivant;
          face_bord = face_suivante;
        }
    }
  return continuer;
}


/*! @brief Applique un vecteur deplacement aux noeuds dont le numero est dans "liste_noeud", puis echange les espaces virtuels.
 *
 *   Si un noeud traverse un joint, il change de proprietaire.
 *   Si un noeud rencontre une paroi, il s'arrete sur la paroi et on le
 *   transforme en noeud "ligne de contact".
 *   Si un noeud est sur une paroi (noeud "ligne de contact"), il longe
 *   la paroi (on deplace le noeud de face de bord en face de bord en
 *   minimisant la distance entre le noeud et le deplacement demande.
 *
 * @param (liste_sommets_initiale) une liste non redondante de noeuds REELS a deplacer
 * @param (deplacement_initial) le vecteur deplacement des noeuds de "liste_noeuds" ( dimension(0)==liste_noeuds.size_array() et dimension(1)==dimension )
 */
void Maillage_FT_Disc::deplacer_sommets(const ArrOfInt& liste_sommets_initiale,
                                        const DoubleTab& deplacement_initial,
                                        ArrOfInt& liste_sommets_sortis,
                                        ArrOfInt& numero_face_sortie,
                                        long skip_facettes)

{
  assert(deplacement_initial.dimension(0) == liste_sommets_initiale.size_array());
  assert(deplacement_initial.dimension(1) == Objet_U::dimension);

  if (Comm_Group::check_enabled()) check_mesh(1,0,skip_facettes);

  const long dimension3 = (Objet_U::dimension == 3);
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, refdomaine_dis_.valeur().valeur());
  const Parcours_interface& parcours = refparcours_interface_.valeur();
  const IntTab& face_voisins = domaine_vf.face_voisins();

  //const ArrOfInt& faces_doubles = zone_vf.faces_doubles();
  //Cerr << "Maillage_FT_Disc::deplacer_sommets : faces_doubles " << faces_doubles.size_array()<< finl;
  //Cerr << "nb_sommets : "<<sommets_.dimension(0)<<  "\tnb_sommets_elem : " << sommet_elem_.size_array()<< finl;
  //Cerr << "sommet_elem_\n" << sommet_elem_<< "\n"<<finl;
  //Cerr << "sommets_\n" << sommets_<< finl;
  //
  // ALGORITHME :
  // On deplace tous les noeuds jusqu'a ce qu'ils rencontrent un joint.
  // Puis on change les noeuds de processeurs et on continue le deplacement
  // des noeuds qui viennent d'arriver.
  //

  // Liste de sommets et deplacement a l'iteration courante
  static ArrOfIntFT  liste_sommets;
  static DoubleTabFT deplacement;
  // Liste des sommets qui atteignent un joint et donnees a envoyer
  static ArrOfIntFT  sommets_envoyes;
  static ArrOfIntFT  element_virtuel_arrivee;
  static ArrOfIntFT  face_virtuelle_arrivee;
  static DoubleTabFT deplacement_restant;

  //On reinitialise sommet_face_bord_  a -1 dans le cas des marqueurs (skip_facettes)
  if (is_solid_particle_) // HMS : suppression des marqueurs de bords // EB : is_solid_particle initialise dans Transport_Interfaces_FT_Disc::discretiser
    sommet_face_bord_ = -1;
  else if (skip_facettes)
    sommet_face_bord_ = -1;


  liste_sommets = liste_sommets_initiale;
  {
    // Remplissage du tableau "deplacement" qui a toujours 3 colonnes.
    const long n = deplacement_initial.dimension(0);
    deplacement.resize(n, 3);
    for (long i = 0; i < n; i++)
      {
        deplacement(i, 0) = deplacement_initial(i, 0);
        deplacement(i, 1) = deplacement_initial(i, 1);
        deplacement(i, 2) = (dimension3) ? deplacement_initial(i, 2) : 0.;
      }
  }

  // Boucle : "tant qu'il reste des noeuds a deplacer sur l'un des processeurs"
  long somme_nb_sommets_envoyes = 0;
  const long max_iterations_echange = 50;
  long nb_iterations_echange = 0;

  do
    {
      const long nbsommets = liste_sommets.size_array();
      sommets_envoyes.resize_array(0);
      element_virtuel_arrivee.resize_array(0);
      face_virtuelle_arrivee.resize_array(0);
      deplacement_restant.resize(0, 3);

      for (long i_sommet = 0; i_sommet < nbsommets; i_sommet++)
        {
          const long num_sommet = liste_sommets[i_sommet];
          // Coordonnee courante du sommet (on la deplace d'element en element)
          double x0  = sommets_(num_sommet, 0);
          double y0  = sommets_(num_sommet, 1);
          double z0  = (dimension3 ? sommets_(num_sommet, 2) : 0.);
          // Position du point d'arrivee
          const double x1 = x0 + deplacement(i_sommet, 0);
          const double y1 = y0 + deplacement(i_sommet, 1);
          const double z1 = (dimension3) ? (z0 + deplacement(i_sommet, 2)) : 0.;
          // Numero de l'element ou se trouve le sommet
          long element = sommet_elem_[num_sommet];
          // Numero de la face de bord ou se trouve le sommet (si ligne de contact)
          long face_bord = sommet_face_bord_[num_sommet];

          assert(element >= 0); // Pas de sommet virtuel dans la liste !
          long continuer = 0;
          long face_bord_m2, face_bord_m1 = -10;
          // Boucle de deplacement du sommet d'element en element

          do
            {
              face_bord_m2 = face_bord_m1;
              face_bord_m1 = face_bord;
              continuer = Maillage_FT_Disc::deplacer_un_sommet(x0, y0, z0,
                                                               x1, y1, z1,
                                                               element,
                                                               face_bord,
                                                               num_sommet,
                                                               parcours,
                                                               domaine_vf,
                                                               face_voisins,
                                                               sommets_envoyes,
                                                               element_virtuel_arrivee,
                                                               face_virtuelle_arrivee,
                                                               deplacement_restant,
                                                               skip_facettes);

              //ajout test pour limiter les allers-retours antre plsieurs faces de bord
              //pb des coins
              if (face_bord!=-1 && face_bord==face_bord_m2)
                {
                  //Process::Journal()<<"face_bord==face_bord_m2 : fin deplacement"<<finl;
                  continuer = 0;
                }
            }
          while (continuer);

          sommets_(num_sommet, 0) = x0;
          sommets_(num_sommet, 1) = y0;
          if (dimension3)
            sommets_(num_sommet, 2) = z0;
          sommet_elem_[num_sommet] = element;
          sommet_face_bord_[num_sommet] = face_bord;
          ///////////////////////////////////////////////////////////////////////////
          //debut  EB
          /*
          IntVect check_som_in_face(2*dimension);
          //const long nb_faces_reelles=zone_vf.nb_faces(); // EB
          for (long dim=0; dim<dimension; dim++) sommet_face_(num_sommet,dim)=-1;

          long mon_elem=-1;
          long elem_virt=-1;
          mon_elem=sommet_elem_(num_sommet); // Si le sommet m'appartient, alors je recupere l'element qui le contient -> je remonte aux faces par rapport a sa position relative au cg de l'element
          //Cerr << "mon_elem " << mon_elem << finl;
          if (mon_elem>=0) check_som_in_face=1;
          else // Le sommet ne m'appartient pas. Je ne connais pas l'element qui le contient -> je le cherche
            {
              elem_virt=zone_vf.zone().chercher_elements(sommets_(num_sommet,0),sommets_(num_sommet,1),sommets_(num_sommet,2)); // element eulerien contenant le sommet virtuel
              assert(elem_virt>=0);
              // Si le sommet est dans le volume de controle d'une face double, on met check_som_in_face a 1.
              check_som_in_face=0;
              for (long dim=0; dim<2*dimension; dim++) // on parcourt les faces du volume de controle de l'element. Si une des faces est reelles, alors le sommet est dans le volume de controle de la premiere couche de joint
                {
                  long la_face=zone_vf.elem_faces(elem_virt,dim);

                  if (faces_doubles(la_face))  // la face est reelle mais l'element est virtuel -> face_double
                    {
                      check_som_in_face(dim)=1;
                      mon_elem=elem_virt;
                    }
                }
            }
          if (mon_elem>=0)
            {

              IntVect faces_elem_eulerien(2*dimension);
              for (long dim=0; dim<2*dimension; dim++) faces_elem_eulerien(dim) = zone_vf.elem_faces(mon_elem,dim); // on recupere les faces de l'element eulerien
              double pos_sommet, coord_elem;
              for (long dim=0; dim<dimension; dim++)
                {
                  pos_sommet=sommets_(num_sommet,dim); // on recupere la position du sommet suivant l'axe "dim"
                  coord_elem=zone_vf.xp(mon_elem,dim);
                  if (pos_sommet < coord_elem && check_som_in_face(dim)) sommet_face_(num_sommet,dim)=faces_elem_eulerien(dim);
                  else if (pos_sommet >= coord_elem && check_som_in_face(dim+dimension)) sommet_face_(num_sommet,dim)=faces_elem_eulerien(dim+dimension);
                  assert(sommet_face_(num_sommet,dim)>=0);
                }
            }
          */
          ////////////////////////////////////////////////////////////////////////////
          // fin EB
        }
      //Cerr << "sommet_elem_\n" << sommet_elem_<< "\n"<<finl;
      //Cerr << "sommets_\n" << sommets_<< finl;
      //Cerr <<"\n" << finl;
      // Mise a jour des espaces virtuels des coordonnees des sommets
      desc_sommets_.echange_espace_virtuel(sommets_);
      // Mise a jour des espaces virtuels des face_bord (si le sommet est ligne
      // de contact, sommet_face_bord_ doit etre >= 0 meme si c'est un sommet virtuel)
      desc_sommets_.echange_espace_virtuel(sommet_face_bord_);

      const long nb_sommets_envoyes = sommets_envoyes.size_array();
      somme_nb_sommets_envoyes = mp_sum(nb_sommets_envoyes);
      if (somme_nb_sommets_envoyes > 0)
        {
          // Echange des sommets et des deplacements restants
          // On recupere une nouvelle "liste_sommets" et un nouveau "deplacement"

          echanger_sommets_PE(sommets_envoyes,
                              element_virtuel_arrivee,
                              face_virtuelle_arrivee,
                              deplacement_restant,
                              liste_sommets,
                              deplacement,
                              skip_facettes);
        }
      nb_iterations_echange++;
    }
  while (somme_nb_sommets_envoyes > 0 && nb_iterations_echange < max_iterations_echange);
  if (nb_iterations_echange == max_iterations_echange)
    {
      if (je_suis_maitre())
        Cerr << "Maillage_FT_Disc.cpp::deplacer_sommets max_iterations_echange atteint." << finl;
    }

  if (Comm_Group::check_enabled()) check_sommets();

  liste_sommets.resize_array(0);
  deplacement.resize(0);
  sommets_envoyes.resize_array(0);
  element_virtuel_arrivee.resize_array(0);
  face_virtuelle_arrivee.resize_array(0);
  deplacement_restant.resize(0);
}

// Description :
//  Verification de la coherence de la structure de donnees des sommets.
long Maillage_FT_Disc::check_sommets(long error_is_fatal) const
{
  const double invalid_value = DMAXFLOAT*0.9;
  const long dimension3 = (Objet_U::dimension == 3);
  const long moi = Process::me();
  const Domaine_dis& domaine_dis = refdomaine_dis_.valeur();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis.valeur());
  const long nb_elements_reels = domaine_vf.nb_elem();
  long i, j;

  if (statut_ == RESET)
    {
      long ok = (nb_sommets() == 0);
      ok = ok && (sommet_elem_.size_array() == 0);
      ok = ok && (sommet_face_.size_array() == 0); // EB
      ok = ok && (sommet_face_bord_.size_array() == 0);
      ok = ok && (sommet_PE_owner_.size_array() == 0);
      ok = ok && (sommet_num_owner_.size_array() == 0);
      ok = ok && (desc_sommets_.espace_virtuel().pe_voisins().size_array() == 0);
      ok = ok && (desc_sommets_.espace_distant().pe_voisins().size_array() == 0);
      if (!ok && error_is_fatal)
        {
          Journal() << "Erreur Maillage_FT_Disc::check_sommets : maillage RESET invalide" << finl;
          assert(0);
          exit();
        }
      return !ok;
    }
  //Journal() << "Entree dans Maillage_FT_Disc::check_sommets" << finl;
  // Verification que les espaces distants et virtuels sont coherents:
  desc_sommets_.check();
  // Verification des tailles des tableaux:
  const long nbsommets = sommets_.dimension(0);
  if (sommets_.dimension(1) != Objet_U::dimension)
    Journal() << "Erreur sommets_.dimension(1) != Objet_U::dimension" << finl;
  // Verification de sommet_PE_owner_ : on le recalcule et on compare
  {
    if (sommet_PE_owner_.size_array() != nbsommets)
      {
        if (error_is_fatal)
          {
            assert(0);
            exit();
          }
        Cerr << "" << finl;
      }

    ArrOfIntFT pe_owner(nbsommets);
    desc_sommets_.remplir_element_pe(pe_owner);
    for (i = 0; i < nbsommets; i++)
      {
        if (sommet_PE_owner_[i] != pe_owner[i])
          {
            if (error_is_fatal)
              {
                assert(0);
                exit();
              }
            break;
          }
      }
    if (i < nbsommets)
      Cerr << "Erreur Verification de sommet_PE_owner_" << finl;
  }
  // Verification de sommet_num_owner_ : on le recalcule et on compare
  {
    if (sommet_num_owner_.size_array() != nbsommets)
      {
        Journal() << "Erreur sommet_num_owner_.size_array() != nb_sommets" << finl;
        if (error_is_fatal)
          {
            assert(0);
            exit();
          }
      }
    ArrOfIntFT num_owner(nbsommets);
    for (i = 0; i < nbsommets; i++)
      num_owner[i] = i;
    desc_sommets_.echange_espace_virtuel(num_owner);
    for (i = 0; i < nbsommets; i++)
      {
        if (num_owner[i] != sommet_num_owner_[i])
          {
            Journal() << "Erreur num_owner[" << i << "] = " << num_owner[i];
            Journal() << "  sommet_num_owner_ = " << sommet_num_owner_[i] << finl;
            if (error_is_fatal)
              {
                assert(0);
                exit();
              }
          }
      }
  }

  // Verification des coordonnees des sommets virtuels : comparaison avec une copie echangee
  {
    DoubleTabFT copie_sommets = sommets_;
    // On invalide les elements virtuels. Cela permet de detecter le cas
    // ou un element se trouverait a la fois dans l'espace distant et dans
    // l'espace virtuel.
    const Descripteur_FT& espace_virtuel = desc_sommets_.espace_virtuel();
    const ArrOfInt& pe_voisins = espace_virtuel.pe_voisins();
    const long nb_pe_voisins = pe_voisins.size_array();
    for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
      {
        const long pe = pe_voisins[indice_pe];
        const ArrOfInt& elements = espace_virtuel.elements(pe);
        const long n = elements.size_array();
        for (i = 0; i < n; i++)
          {
            const long num_sommet = elements[i];
            copie_sommets(num_sommet, 0) = invalid_value;
            copie_sommets(num_sommet, 1) = invalid_value;
            if (dimension3)
              copie_sommets(num_sommet, 2) = invalid_value;
          }
      }
    // Echange espace virtuel sur la copie
    desc_sommets_.echange_espace_virtuel(copie_sommets);
    // Comparaison de la copie et du tableau des sommets.
    for (i = 0; i < nbsommets; i++)
      {
        for (j = 0; j < Objet_U::dimension; j++)
          {
            if (copie_sommets(i,j) == invalid_value)
              {
                Journal() << "Erreur copie_sommets(" << i << ",";
                Journal() << j << ") == DMAX_FLOAT" << finl;
                if (error_is_fatal)
                  {
                    assert(0);
                    exit();
                  }
              }
            if (copie_sommets(i,j) != sommets_(i,j))
              {
                Journal() << "Erreur copie_sommets(" << i << ",";
                Journal() << j << ") = " << copie_sommets(i,j);
                Journal() << " sommets_ = " << sommets_(i,j) << finl;
                if (error_is_fatal)
                  {
                    assert(0);
                    exit();
                  }
              }
          }
      }
  }

  // Verification de sommet_elem_ :
  {
    if (sommet_elem_.size_array() != nbsommets)
      {
        Journal() << "Erreur sommet_elem_.size_array() != nb_sommets" << finl;
        if (error_is_fatal)
          {
            assert(0);
            exit();
          }
      }
    const Parcours_interface& parcours =refparcours_interface_.valeur();
    const double epsilon = parcours.get_erreur_geometrique();

    for (i = 0; i < nbsommets; i++)
      {
        const long num_element = sommet_elem_[i];
        const long pe = sommet_PE_owner_[i];
        if (pe != moi && num_element >= 0)
          {
            Journal() << "Erreur sommet_num_owner_[" << i;
            Journal() << "] && sommet_elem_ >= 0" << finl;
            if (error_is_fatal)
              {
                assert(0);
                exit();
              }
          }
        // EB : on se passe de cette verif pour sommet_face car ce tableau est construit a partir de sommet_elem
        // Et avec les faces doubles, c'est possible d'avoir pe!=moi && num_face >=0 pour un sommet virtuel
        if (pe == moi)
          {
            if (num_element > nb_elements_reels)
              {
                Journal() << "Erreur i=" << i;
                Journal() << "   sommet_elem_[i]>nb_elements_reels" << finl;
                if (error_is_fatal)
                  {
                    assert(0);
                    exit();
                  }
              }
            else
              {
                // On verifie que le sommet est dans l'element a epsilon pres.
                // Le transport assure une distance < epsilon. On laisse un peu de marge,
                // on teste si distance < 2. * epsilon.
                double x, y, z = 0.;
                x = sommets_(i, 0);
                y = sommets_(i, 1);
                if (dimension3)
                  z = sommets_(i, 2);
                const double distance = parcours.distance_sommet_faces(domaine_vf, num_element,
                                                                       x, y, z);
                if (distance > 2. * epsilon)
                  {
                    Journal() << "Erreur i=" << i << "  num_elemen=" << num_element;
                    Journal() << "  distance_sommet_faces=" << distance << finl;
                    if (error_is_fatal)
                      {
                        Cerr << "Trouble when checking mesh in Maillage_FT_Disc::check_sommets." << finl;
                        Cerr << "Contact TRUST support." << finl;
                        exit();
                      }
                  }
              }
          }
      }
    // Verification de sommet_face_bord_ :
    {
      if (sommet_face_bord_.size_array() != nbsommets)
        {
          Journal() << "Erreur sommet_face_bord_.size_array() != nb_sommets" << finl;
          if (error_is_fatal)
            {
              assert(0);
              exit();
            }
        }
      const IntTab& face_voisins = domaine_vf.face_voisins();
      for (i = 0; i < nbsommets; i++)
        {
          const long face = sommet_face_bord_[i];
          const long element = sommet_elem_[i];
          if (face >= 0 && !sommet_virtuel(i))
            {
              const long voisin0 = face_voisins(face, 0);
              const long voisin1 = face_voisins(face, 1);
              if (voisin0 >= 0 && voisin1 >= 0)
                {
                  Journal() << "Erreur sommet_face_bord_[" << i;
                  Journal() << "]=" << face << " (la face n'est pas au bord)." << finl;
                  if (error_is_fatal)
                    {
                      assert(0);
                      exit();
                    }
                }
              if (voisin0 != element && voisin1 != element)
                {
                  Journal() << "Erreur sommet_face_bord_[" << i;
                  Journal() << "]=" << face;
                  Journal() << "\n   (la face n'est pas voisine de l'element ";
                  Journal() << element << finl;
                  if (error_is_fatal)
                    {
                      assert(0);
                      exit();
                    }
                }
            }
        }
    }
  }
  return 0;
}

static True_int fct_tri_facettes(const void *pt1, const void *pt2)
{
  const long *a = (const long *) pt1;
  const long *b = (const long *) pt2;

  long i, x = 0;
  const long dim = Objet_U::dimension;
  for (i = 0; i < dim; i++)
    {
      x = a[i] - b[i];
      if (x != 0)
        break;
    }
  return x;
}

// Description
//  Verifie la coherence des structures de donnees de l'interface (notamment
//  la distribution des donnees sur les differents processeurs, espaces
//  virtuels, etc). On verifie que toutes les conditions imposees pour avoir
//  un maillage dans l'etat MINIMAL sont remplies.
//
//  Si skip_facette_pe != 0, on ne verifie pas la condition "proprietaire facette
//  == proprietaire premier sommet".
//  Si skip_facettes != 0    on ne considere que check_sommets()

long Maillage_FT_Disc::check_mesh(long error_is_fatal, long skip_facette_pe, long skip_facettes) const
{
  static const Stat_Counter_Id stat_counter = statistiques().new_counter(3, "Check_mesh", "FrontTracking");
  statistiques().begin_count(stat_counter);

  const double invalid_value = DMAXFLOAT*0.9;
  long i, j;
  long return_code = -1;

  return_code = check_sommets(error_is_fatal);
  if (return_code == 0)
    return_code = -1;

  if (skip_facettes)
    return return_code;

  if (return_code < 0 && statut_ == RESET)
    {
      long ok = (nb_facettes() == 0);
      ok = ok && (facette_num_owner_.size_array() == 0);
      ok = ok && (desc_facettes_.espace_virtuel().pe_voisins().size_array() == 0);
      ok = ok && (desc_facettes_.espace_distant().pe_voisins().size_array() == 0);
      if (!ok && error_is_fatal)
        {
          Journal() << "Erreur Maillage_FT_Disc::check_mesh : maillage RESET invalide" << finl;
          assert(0);
          exit();
        }
      return_code = !ok;
    }

  if (return_code < 0)
    {
      desc_facettes_.check();
      // Verification des tailles des tableaux:
      const long nbsommets = sommets_.dimension(0);
      const long nbfacettes = facettes_.dimension(0);
      if (facettes_.dimension(1) != Objet_U::dimension)
        Journal() << "Erreur facettes_.dimension(1) != Objet_U::dimension" << finl;

      // Verification des facettes :

      // Les numeros de sommets existent-t-ils ?
      // Les espaces virtuels sont-ils coherents ?
      // On verifie que les coordonnees des sommets sont les memes sur tous les procs.
      // Pour ca, on cree un tableau contenant les coordonnees des sommets des facettes,
      // on echange l'espace virtuel du tableau et on compare au tableau facettes_ et sommets_
      const long dim_carre = Objet_U::dimension * Objet_U::dimension;
      DoubleTabFT coord_facettes(nbfacettes, dim_carre);

      for (i = 0; i < nbfacettes; i++)
        {
          long n = 0;
          for (j = 0; j < Objet_U::dimension; j++)
            {
              const long som = facettes_(i, j);
              if (som >= nbsommets)
                {
                  Journal() << "Erreur facettes_(" << i << "," << j;
                  Journal() << ") >= nb_sommets" << finl;
                  printFa7(i,1,Journal());
                  if (error_is_fatal)
                    {
                      assert(0);
                      exit();
                    }
                  return_code = 0;
                }
              for (long k = 0; k < Objet_U::dimension; k++)
                coord_facettes(i, n++) = sommets_(som, k);
            }
        }
      // On invalide les espaces virtuels:
      const Descripteur_FT& espace_virtuel = desc_facettes_.espace_virtuel();
      const ArrOfInt& pe_voisins = espace_virtuel.pe_voisins();
      const long nb_pe_voisins = pe_voisins.size_array();
      for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
        {
          const long pe = pe_voisins[indice_pe];
          const ArrOfInt& elements = espace_virtuel.elements(pe);
          const long n = elements.size_array();
          for (i = 0; i < n; i++)
            {
              const long num_facette = elements[i];
              for (j = 0; j < dim_carre; j++)
                coord_facettes(num_facette, j) = invalid_value;
            }
        }
      desc_facettes_.echange_espace_virtuel(coord_facettes);
      // Verification des coordonnees
      for (i = 0; i < nbfacettes; i++)
        {
          long n = 0;
          for (j = 0; j < Objet_U::dimension; j++)
            {
              const long som = facettes_(i, j);
              for (long k = 0; k < Objet_U::dimension; k++)
                {
                  double original = sommets_(som, k);
                  double copie = coord_facettes(i, n++);
                  if (original != copie)
                    {
                      Journal() << "Erreur facette " << i << " sommet " << j << " direction " << k;
                      Journal() << "   " << original << " != " << copie << finl;
                      printFa7(i,1,Journal());
                      Journal() <<" Fa7 orig ="<<finl;
                      if (Objet_U::dimension ==3)
                        {
                          Journal() << " som0= " << sommets_(i, 0) << " " << sommets_(i, 1) << " " << sommets_(i, 2) << finl;
                          Journal() << " som1= " << sommets_(i, 3) << " " << sommets_(i, 4) << " " << sommets_(i, 5) << finl;
                          Journal() << " som2= " << sommets_(i, 6) << " " << sommets_(i, 7) << " " << sommets_(i, 8) << finl;
                        }
                      else
                        {
                          Journal() << " som0= " << sommets_(i, 0) << " " << sommets_(i, 1) <<finl;
                          Journal() << " som1= " << sommets_(i, 2) << " " << sommets_(i, 3) << finl;
                        }
                      Journal() <<" Fa7 copie ="<<finl;
                      if (Objet_U::dimension ==3)
                        {
                          Journal() << " som0= " << coord_facettes(i, 0) << " " << coord_facettes(i, 1) << " " << coord_facettes(i, 2) << finl;
                          Journal() << " som1= " << coord_facettes(i, 3) << " " << coord_facettes(i, 4) << " " << coord_facettes(i, 5) << finl;
                          Journal() << " som2= " << coord_facettes(i, 6) << " " << coord_facettes(i, 7) << " " << coord_facettes(i, 8) << finl;
                        }
                      else
                        {
                          Journal() << " som0= " << coord_facettes(i, 0) << " " << coord_facettes(i, 1) <<finl;
                          Journal() << " som1= " << coord_facettes(i, 2) << " " << coord_facettes(i, 3) << finl;
                        }
                      if (error_is_fatal)
                        {
                          Process::exit();
                        }
                      return_code = 0;
                    }
                }
            }
        }

      // Verification qu'il n'existe pas deux fois la meme facette
      {
        IntTabFT copie_facettes = facettes_;
        // tri du tableau
        long * data = copie_facettes.addr();
        const long nbr_facettes = facettes_.dimension(0);
        assert(Objet_U::dimension == facettes_.dimension(1));
        qsort(data, nbr_facettes, facettes_.dimension(1)*sizeof(long),
              fct_tri_facettes);
        // recherche des doublons
        long ii, jj;
        long count = 0;
        const long nb_som_facettes = Objet_U::dimension;
        for (ii = 1; ii < nbr_facettes; ii++)
          {
            long facette_a_supprimer = (copie_facettes(ii,0) == copie_facettes(ii,1));
            if (! facette_a_supprimer)
              {
                for (jj = 0; jj < nb_som_facettes; jj++)
                  {
                    if (copie_facettes(ii,jj) != copie_facettes(ii-1,jj))
                      break;
                  }
                if (jj == nb_som_facettes)
                  {
                    count++;
                  }
              }
          }
        if (count > 0)
          {
            Cerr << "Erreur facette : " << count << " facettes identiques sur PE " << me() << finl;
          }
      }

      // Verification du proprietaire de la facette.
      // On reconstruit le tableau des proprietaires d'apres l'espace virtuel,
      // on verifie que c'est bien le proprietaire du premier sommet de la facette.
      // (code identique au debut de "corriger_proprietaire_facette").
      if (! skip_facette_pe)
        {
          // Quel est le proprietaire actuel des facettes ?
          // (d'apres le descripteur des facettes)
          ArrOfIntFT facette_pe(nbfacettes);
          desc_facettes_.remplir_element_pe(facette_pe);
          // On verifie:
          for (i = 0; i < nbfacettes; i++)
            {
              long pe_actuel = facette_pe[i];
              long premier_sommet = facettes_(i, 0);
              long pe_legitime = sommet_PE_owner_[premier_sommet];
              if (pe_actuel != pe_legitime)
                {
                  Journal() << "Erreur facette " << i;
                  Journal() << " : owner(d'apres descripteur)=" << pe_actuel;
                  Journal() << "   owner(d'apres 1er sommet)=" << pe_legitime << finl;
                  printFa7(i,1,Journal());
                  if (error_is_fatal)
                    {
                      assert(0);
                      exit();
                    }
                  return_code = 0;
                }
            }
          // Verification de facette_num_owner_
          ArrOfIntFT num_owner(nbfacettes);
          for (i = 0; i < nbfacettes; i++)
            num_owner[i] = i;
          desc_facettes().echange_espace_virtuel(num_owner);
          for (i = 0; i < nbfacettes; i++)
            if (num_owner[i] != facette_num_owner_[i])
              {
                Journal() << "Erreur facette " << i;
                Journal() << " : facette_num_owner[i]=" << facette_num_owner_[i];
                Journal() << " devrait valoir " << num_owner[i] << finl;
                if (error_is_fatal)
                  {
                    assert(0);
                    exit();
                  }
                return_code = 0;
              }
        }
    }
  statistiques().end_count(stat_counter);

  return return_code;
}

/*! @brief Retire toutes les facettes virtuelles et tous les sommets qui ne sont pas utilises.
 *
 */
void Maillage_FT_Disc::nettoyer_elements_virtuels()
{
  //lance le nettoyage
  nettoyer_maillage();
  // maillage_modifie() est fait dans nettoyer_maillage
}

/*! @brief Retire toutes les facettes virtuelles, toutes les facettes invalides (sommet0 == sommet1) et tous les sommets qui ne sont pas utilises.
 *
 */
void Maillage_FT_Disc::nettoyer_maillage()
{
  Cerr << "Maillage_FT_Disc::nettoyer_maillage " << finl;
  assert(statut_ >= MINIMAL);
  // Lorsque l'on arrive par Maillage_FT_IJK::nettoyer_maillage et qu'il y a
  // eu un nettoyage preliminaire de compo_connexe_facettes_, l'appel :
  // if (Comm_Group::check_enabled()) check_mesh(1,0,1);
  // va pointer vers Maillage_FT_IJK::check_mesh et va tester la taille du tableau
  // compo_connexe. A ce moment, le tableau compo_connexe_facettes_ a ete nettoye
  // de ces facettes virtuel alors que le maillage ne va l'etre que par la suite.
  // On force donc ici l'appel a Maillage_FT_Disc::check_mesh.
  // La verification par un Maillage_FT_IJK::check_mesh aura lieu a la fin de cette routine.
  if (Comm_Group::check_enabled()) Maillage_FT_Disc::check_mesh(1,0,1);

  const long dimension3 = (Objet_U::dimension==3);
  ArrOfIntFT new_elements;

  // On decale toutes les facettes reelles pour remplir les trous
  // dans les tableaux :
  // * facettes_
  // * facette_num_owner_
  {
    const long nbfacettes = facettes_.dimension(0);
    long n = 0;
    long i;
    for (i = 0; i < nbfacettes; i++)
      {
        const long invalide = (facettes_(i,0) == facettes_(i,1));
        const long virtuelle = facette_virtuelle(i);
        if (!invalide && !virtuelle)
          {
            facettes_(n, 0) = facettes_(i, 0);
            facettes_(n, 1) = facettes_(i, 1);
            if (dimension3)
              facettes_(n, 2) = facettes_(i, 2);
            n++;
          }
      }
    facettes_.resize(n, Objet_U::dimension);

    facette_num_owner_.resize_array(n);
    for (i = 0; i < n; i++)
      facette_num_owner_[i] = i;
    // Inutile d'echanger les espaces virtuels de facette_num_owner_:
    // il n'y a plus d'espaces virtuels.
  }

  // Il n'y a plus aucune facette virtuelle, on vide les espaces distants et
  // virtuels:
  desc_facettes_.reset();
  desc_facettes_.calcul_schema_comm(facettes_.dimension(0));

  // Marquage des sommets utilises:
  ArrOfIntFT sommets_utilises(nb_sommets());
  {
    ArrOfInt& tab = facettes_;  // vu comme unidimensionnel
    const long nb_sommets_facettes = tab.size_array();
    sommets_utilises = 0;
    // Boucle sur tous les sommets de toutes les facettes
    for (long i = 0; i < nb_sommets_facettes; i++)
      {
        const long sommet = tab[i];
        sommets_utilises[sommet]++;
      }
  }

  // Suppression des sommets virtuels inutilises
  // Pour cela, on parcours l'espace virtuel des sommets,
  // et si un sommet n'est pas utilise on le retire de l'espace virtuel
  // et de l'espace distant correspondant (on envoie au proprietaire
  // du sommet le rang du sommet a supprimer dans le tableau des elements
  // distants).
  // Le schema de comm a utiliser est "les proprietaires d'elements virtuels
  //  parlent aux proprietaires des elements reels" :
  {
    const Schema_Comm_FT& comm = desc_sommets_.schema_comm_inverse();
    comm.begin_comm();
    // Boucle sur les espaces virtuels de sommets
    {
      Descripteur_FT& espace_virtuel = desc_sommets_.espace_virtuel();
      const ArrOfInt& pe_voisins = espace_virtuel.pe_voisins();
      const long nb_pe_voisins = pe_voisins.size_array();
      for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
        {
          const long pe = pe_voisins[indice_pe];
          const ArrOfInt& elements = espace_virtuel.elements(pe); // EB : les elements sont les sommets lagrangiens. Voir Descipteur_FT.h
          const long nb_elements = elements.size_array(); // EB : nb_elements = nombre de sommets distants/virtuels du pe
          Sortie& buffer = comm.send_buffer(pe);
          new_elements.resize_array(0);
          for (long i = 0; i < nb_elements; i++)
            {
              const long sommet = elements[i];
              if (sommets_utilises[sommet])
                {
                  // Le sommet est utilise, on le laisse dans la liste
                  new_elements.append_array(sommet);
                }
              else
                {
                  // Le sommet n'est pas utilise, on envoie le rang du sommet a supprimer
                  // au proprietaire
                  buffer << i;
                }
            }
          // Mise a jour des elements virtuels
          espace_virtuel.set_elements(pe, new_elements);
        }
      espace_virtuel.calcul_liste_pe_voisins();
    }
    comm.echange_taille_et_messages();
    // On retire les elements des listes d'elements distants
    {
      Descripteur_FT& espace_distant = desc_sommets_.espace_distant();
      const ArrOfInt& pe_voisins = comm.get_recv_pe_list();
      const long nb_pe_voisins = pe_voisins.size_array();
      for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
        {
          const long pe = pe_voisins[indice_pe];
          const ArrOfInt& elements = espace_distant.elements(pe);
          const long nb_elements = elements.size_array();
          Entree& buffer = comm.recv_buffer(pe);
          long element_courant = 0;
          new_elements.resize_array(0);
          do
            {
              // rang du prochain sommet a supprimer dans "elements"
              long rang_sommet;
              buffer >> rang_sommet;
              if (buffer.eof())
                rang_sommet = nb_elements;
              else
                assert(rang_sommet < nb_elements);
              // On conserve les sommets jusqu'a rang_sommet exclu
              for (; element_courant < rang_sommet; element_courant++)
                {
                  const long sommet = elements[element_courant];
                  new_elements.append_array(sommet);
                  // On marque du sommet car il est utilise
                  // dans un espace distant
                  sommets_utilises[sommet]++;
                }
              // On passe le sommet a retirer
              element_courant++;
            }
          while (element_courant < nb_elements);
          // Mise a jour des elements distants
          espace_distant.set_elements(pe, new_elements);
        }
      espace_distant.calcul_liste_pe_voisins();
    }
    comm.end_comm();
  }
  // Suppression des sommets reels inutilises :
  // Mise a jour de :
  //  * sommets_
  //  * sommet_elem_
  //  * sommet_face_
  //  * sommet_face_bord_
  //  * sommet_PE_owner_
  //  * drapeaux_sommets_
  // On construit en meme temps le tableau renum_sommets:
  // si i est l'indice actuel d'un sommet, renum_sommets[i]
  // est l'indice du sommet apres suppression des sommets
  // inutilises.
  ArrOfIntFT renum_sommets(sommets_.dimension(0));
  {
    const long nbsommets = sommets_.dimension(0);
    // Compteur de sommets apres suppression
    long n = 0;
    for (long i = 0; i < nbsommets; i++)
      {
        if (sommets_utilises[i])
          {
            renum_sommets[i] = n;
            sommets_(n, 0) = sommets_(i, 0);
            sommets_(n, 1) = sommets_(i, 1);
            if (dimension3)
              sommets_(n, 2) = sommets_(i, 2);
            sommet_elem_[n] = sommet_elem_[i];
            if(i<sommet_face_.dimension(0) && n <sommet_face_.dimension(0))
              for (long dim=0; dim<dimension; dim++) sommet_face_(n,dim) = sommet_face_(i,dim); // EB
            sommet_face_bord_[n] = sommet_face_bord_[i];
            sommet_PE_owner_[n] = sommet_PE_owner_[i];
            drapeaux_sommets_[n] = drapeaux_sommets_[i];
            n++;
          }
        else
          {
            renum_sommets[i] = -1; // Le sommet i est supprime
          }
      }
    sommets_         .resize(n, Objet_U::dimension);
    sommet_elem_     .resize_array(n);
    sommet_face_     .resize(n,Objet_U::dimension); // EB
    sommet_face_bord_.resize_array(n);
    sommet_PE_owner_ .resize_array(n);
    drapeaux_sommets_.resize_array(n);
  }

  // Renumerotation des sommets dans desc_sommets
  {
    // Boucle sur les deux espaces : distant et virtuel
    for (long num_espace = 0; num_espace < 2; num_espace++)
      {
        // Espace_dv est soit l'ESPACE_D(istant), soit l'ESPACE_V(irtuel).
        Descripteur_FT& espace_dv =
          (num_espace == 0)
          ? desc_sommets_.espace_distant()
          : desc_sommets_.espace_virtuel();

        const ArrOfInt& pe_voisins = espace_dv.pe_voisins();
        const long nb_pe_voisins = pe_voisins.size_array();
        for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
          {
            const long pe = pe_voisins[indice_pe];
            const ArrOfInt& elements = espace_dv.elements(pe);
            const long nb_elements = elements.size_array();
            new_elements.resize_array(nb_elements);
            for (long i = 0; i < nb_elements; i++)
              {
                const long num_sommet = elements[i];
                const long new_num = renum_sommets[num_sommet];
                assert(new_num >= 0);
                new_elements[i] = new_num;
              }
            espace_dv.set_elements(pe, new_elements);
          }
        // Mise a jour du descripteur:
        espace_dv.calcul_liste_pe_voisins();
      }
    desc_sommets_.calcul_schema_comm(sommets_.dimension(0));
  }

  // Mise a jour de sommet_num_owner_:
  {
    const long nbsommets = sommets_.dimension(0);
    sommet_num_owner_.resize_array(nbsommets);
    for (long i = 0; i < nbsommets; i++)
      sommet_num_owner_[i] = i;
    desc_sommets_.echange_espace_virtuel(sommet_num_owner_);
  }

  // Renumerotation des sommets dans facettes_
  {
    ArrOfInt& tab = facettes_;  // vu comme unidimensionnel
    const long nb_sommets_facettes = tab.size_array();
    // Boucle sur tous les sommets de toutes les facettes
    for (long i = 0; i < nb_sommets_facettes; i++)
      {
        const long sommet = tab[i];
        const long new_num = renum_sommets[sommet];
        tab[i] = new_num;
      }
  }

  maillage_modifie(MINIMAL);

  if (Comm_Group::check_enabled()) check_mesh();
}

/*! @brief Supprime les facettes dont les indices locaux sont donnes en parametre.
 *
 * Le maillage est nettoye et retourne a l'etat MINIMAL.
 *
 */
void Maillage_FT_Disc::supprimer_facettes(const ArrOfInt& liste_facettes)
{
  const long n = liste_facettes.size_array();
  for (long i = 0; i < n; i++)
    {
      long j = liste_facettes[i];
      // On rend la facette j "invalide" pour nettoyer_maillage()
      facettes_(j, 1) = facettes_(j, 0);
    }
  nettoyer_maillage();
}

//Procedure de nettoyage des sommets que l on ne veut pas considerer
//a l etape suivante: sommets virtuels et sur des frontieres ouvertes
void Maillage_FT_Disc::nettoyer_noeuds_virtuels_et_frontieres()
{
  assert(statut_ >= MINIMAL);

  if (Comm_Group::check_enabled()) check_mesh(1,1,1);
  const long dimension3 = (Objet_U::dimension==3);
  ArrOfIntFT new_elements;
  ArrOfIntFT sommets_utilises(nb_sommets());
  sommets_utilises=0;

  //On detecte les sommets utilises
  //On ne retient pas les sommets virtuels et ceux situes sur des faces de frontiere ouverte
  for (long som=0; som<nb_sommets(); som++)
    {
      const Domaine_Cl_dis_base& zcl = equation_transport().get_probleme_base().equation(0).domaine_Cl_dis().valeur();
      long face_loc;
      long face_bord = sommet_face_bord_[som];
      long face_fr_ouverte = 0;
      if ((face_bord!=-1))
        {
          const Cond_lim_base& type_cl = zcl.condition_limite_de_la_face_reelle(face_bord,face_loc);
          if (((!sub_type(Dirichlet,type_cl)) && (!sub_type(Dirichlet_homogene,type_cl)) && (!sub_type(Symetrie,type_cl)))
              || (sub_type(Dirichlet_entree_fluide,type_cl)))
            face_fr_ouverte = 1;
        }

      if ((!sommet_virtuel(som)) && (!face_fr_ouverte))
        sommets_utilises[som] ++;
    }


  // Suppression des sommets virtuels inutilises
  // Pour cela, on parcours l'espace virtuel des sommets,
  // et si un sommet n'est pas utilise on le retire de l'espace virtuel
  // et de l'espace distant correspondant (on envoie au proprietaire
  // du sommet le rang du sommet a supprimer dans le tableau des elements
  // distants).
  // Le schema de comm a utiliser est "les proprietaires d'elements virtuels
  //  parlent aux proprietaires des elements reels" :
  {
    const Schema_Comm_FT& comm = desc_sommets_.schema_comm_inverse();
    comm.begin_comm();
    // Boucle sur les espaces virtuels de sommets
    {
      Descripteur_FT& espace_virtuel = desc_sommets_.espace_virtuel();
      const ArrOfInt& pe_voisins = espace_virtuel.pe_voisins();
      const long nb_pe_voisins = pe_voisins.size_array();
      for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
        {
          const long pe = pe_voisins[indice_pe];
          const ArrOfInt& elements = espace_virtuel.elements(pe);
          const long nb_elements = elements.size_array();
          Sortie& buffer = comm.send_buffer(pe);
          new_elements.resize_array(0);
          for (long i = 0; i < nb_elements; i++)
            {
              const long sommet = elements[i];
              if (sommets_utilises[sommet])
                {
                  // Le sommet est utilise, on le laisse dans la liste
                  new_elements.append_array(sommet);
                }
              else
                {
                  // Le sommet n'est pas utilise, on envoie le rang du sommet a supprimer
                  // au proprietaire
                  buffer << i;
                }
            }
          // Mise a jour des elements virtuels
          espace_virtuel.set_elements(pe, new_elements);
        }
      espace_virtuel.calcul_liste_pe_voisins();
    }
    comm.echange_taille_et_messages();
    // On retire les elements des listes d'elements distants
    {
      Descripteur_FT& espace_distant = desc_sommets_.espace_distant();
      const ArrOfInt& pe_voisins = comm.get_recv_pe_list();
      const long nb_pe_voisins = pe_voisins.size_array();
      for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
        {
          const long pe = pe_voisins[indice_pe];
          const ArrOfInt& elements = espace_distant.elements(pe);
          const long nb_elements = elements.size_array();
          Entree& buffer = comm.recv_buffer(pe);
          long element_courant = 0;
          new_elements.resize_array(0);
          do
            {
              // rang du prochain sommet a supprimer dans "elements"
              long rang_sommet;
              buffer >> rang_sommet;
              if (buffer.eof())
                rang_sommet = nb_elements;
              else
                assert(rang_sommet < nb_elements);
              // On conserve les sommets jusqu'a rang_sommet exclu
              for (; element_courant < rang_sommet; element_courant++)
                {
                  const long sommet = elements[element_courant];
                  new_elements.append_array(sommet);
                  // On marque du sommet car il est utilise
                  // dans un espace distant
                  sommets_utilises[sommet]++;
                }
              // On passe le sommet a retirer
              element_courant++;
            }
          while (element_courant < nb_elements);
          // Mise a jour des elements distants
          espace_distant.set_elements(pe, new_elements);
        }
      espace_distant.calcul_liste_pe_voisins();
    }
    comm.end_comm();
  }

  // Suppression des sommets reels inutilises
  // (dans notre cas ceux sont des sommets sur des faces de frontieres ouvertes):
  // Mise a jour de :
  //  * sommets_
  //  * sommet_elem_
  //  * sommet_face_  // EB : verifier qu'il ne faille pas definir un tableau sommets_utilises_face
  //  * sommet_face_bord_
  //  * sommet_PE_owner_
  //  * drapeaux_sommets_
  // On construit en meme temps le tableau renum_sommets:
  // si i est l'indice actuel d'un sommet, renum_sommets[i]
  // est l'indice du sommet apres suppression des sommets
  // inutilises.
  ArrOfIntFT renum_sommets(sommets_.dimension(0));
  {
    const long nbsommets = sommets_.dimension(0);
    // Compteur de sommets apres suppression
    long n = 0;
    for (long i = 0; i < nbsommets; i++)
      {
        if (sommets_utilises[i])
          {
            renum_sommets[i] = n;
            sommets_(n, 0) = sommets_(i, 0);
            sommets_(n, 1) = sommets_(i, 1);
            if (dimension3)
              sommets_(n, 2) = sommets_(i, 2);
            sommet_elem_[n] = sommet_elem_[i];
            for (long dim=0; dim<dimension; dim++) sommet_face_(n,dim) = sommet_face_(i,dim); // EB
            sommet_face_bord_[n] = sommet_face_bord_[i];
            sommet_PE_owner_[n] = sommet_PE_owner_[i];
            drapeaux_sommets_[n] = drapeaux_sommets_[i];
            n++;
          }
        else
          {
            renum_sommets[i] = -1; // Le sommet i est supprime
          }
      }
    sommets_         .resize(n, Objet_U::dimension);
    sommet_elem_     .resize_array(n);
    sommet_face_     .resize(n, Objet_U::dimension); // EB
    sommet_face_bord_.resize_array(n);
    sommet_PE_owner_ .resize_array(n);
    drapeaux_sommets_.resize_array(n);

    equation_transport().nettoyer_proprietes_particules(sommets_utilises);

  }

  //On supprime ce qui est lie aux espaces distants et virtuel

  // Renumerotation des sommets dans desc_sommets
  {
    // Boucle sur les deux espaces : distant et virtuel
    for (long num_espace = 0; num_espace < 2; num_espace++)
      {
        // Espace_dv est soit l'ESPACE_D(istant), soit l'ESPACE_V(irtuel).
        Descripteur_FT& espace_dv =
          (num_espace == 0)
          ? desc_sommets_.espace_distant()
          : desc_sommets_.espace_virtuel();

        const ArrOfInt& pe_voisins = espace_dv.pe_voisins();
        const long nb_pe_voisins = pe_voisins.size_array();
        for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
          {
            const long pe = pe_voisins[indice_pe];
            const ArrOfInt& elements = espace_dv.elements(pe);
            const long nb_elements = elements.size_array();
            new_elements.resize_array(nb_elements);
            for (long i = 0; i < nb_elements; i++)
              {
                const long num_sommet = elements[i];
                const long new_num = renum_sommets[num_sommet];
                assert(new_num >= 0);
                new_elements[i] = new_num;
              }
            espace_dv.set_elements(pe, new_elements);
          }
        // Mise a jour du descripteur:
        espace_dv.calcul_liste_pe_voisins();
      }
    desc_sommets_.calcul_schema_comm(sommets_.dimension(0));
  }

  // Mise a jour de sommet_num_owner_:
  {
    const long nbsommets = sommets_.dimension(0);
    sommet_num_owner_.resize_array(nbsommets);
    for (long i = 0; i < nbsommets; i++)
      sommet_num_owner_[i] = i;
    desc_sommets_.echange_espace_virtuel(sommet_num_owner_);
  }

  maillage_modifie(MINIMAL);
  if (Comm_Group::check_enabled()) check_mesh(1,1,1);
}

void Maillage_FT_Disc::nettoyer_phase(const Nom& nom_eq, const long phase)
{
  assert(statut_ >= MINIMAL);

  if (Comm_Group::check_enabled()) check_mesh(1,1,1);
  const long dimension3 = (Objet_U::dimension==3);
  ArrOfIntFT new_elements;
  ArrOfIntFT sommets_utilises(nb_sommets());
  sommets_utilises=0;
  const Equation_base& eq = equation_transport().probleme().get_equation_by_name(nom_eq);
  Transport_Interfaces_FT_Disc& eq_interf = ref_cast_non_const(Transport_Interfaces_FT_Disc,eq);
  const DoubleTab& indic =  eq_interf.inconnue().valeurs();
  double phase_reelle = double(phase);
  long elem;

  //Les sommets utilises sont ceux places dans la phase marquee
  for (long som=0; som<nb_sommets(); som++)
    {
      elem = sommet_elem_[som];
      if (est_egal(indic(elem),phase_reelle))
        sommets_utilises[som] ++;
    }


  // Suppression des sommets virtuels inutilises
  // Pour cela, on parcours l'espace virtuel des sommets,
  // et si un sommet n'est pas utilise on le retire de l'espace virtuel
  // et de l'espace distant correspondant (on envoie au proprietaire
  // du sommet le rang du sommet a supprimer dans le tableau des elements
  // distants).
  // Le schema de comm a utiliser est "les proprietaires d'elements virtuels
  //  parlent aux proprietaires des elements reels" :
  {
    const Schema_Comm_FT& comm = desc_sommets_.schema_comm_inverse();
    comm.begin_comm();
    // Boucle sur les espaces virtuels de sommets
    {
      Descripteur_FT& espace_virtuel = desc_sommets_.espace_virtuel();
      const ArrOfInt& pe_voisins = espace_virtuel.pe_voisins();
      const long nb_pe_voisins = pe_voisins.size_array();
      for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
        {
          const long pe = pe_voisins[indice_pe];
          const ArrOfInt& elements = espace_virtuel.elements(pe);
          const long nb_elements = elements.size_array();
          Sortie& buffer = comm.send_buffer(pe);
          new_elements.resize_array(0);
          for (long i = 0; i < nb_elements; i++)
            {
              const long sommet = elements[i];
              if (sommets_utilises[sommet])
                {
                  // Le sommet est utilise, on le laisse dans la liste
                  new_elements.append_array(sommet);
                }
              else
                {
                  // Le sommet n'est pas utilise, on envoie le rang du sommet a supprimer
                  // au proprietaire
                  buffer << i;
                }
            }
          // Mise a jour des elements virtuels
          espace_virtuel.set_elements(pe, new_elements);
        }
      espace_virtuel.calcul_liste_pe_voisins();
    }
    comm.echange_taille_et_messages();
    // On retire les elements des listes d'elements distants
    {
      Descripteur_FT& espace_distant = desc_sommets_.espace_distant();
      const ArrOfInt& pe_voisins = comm.get_recv_pe_list();
      const long nb_pe_voisins = pe_voisins.size_array();
      for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
        {
          const long pe = pe_voisins[indice_pe];
          const ArrOfInt& elements = espace_distant.elements(pe);
          const long nb_elements = elements.size_array();
          Entree& buffer = comm.recv_buffer(pe);
          long element_courant = 0;
          new_elements.resize_array(0);
          do
            {
              // rang du prochain sommet a supprimer dans "elements"
              long rang_sommet;
              buffer >> rang_sommet;
              if (buffer.eof())
                rang_sommet = nb_elements;
              else
                assert(rang_sommet < nb_elements);
              // On conserve les sommets jusqu'a rang_sommet exclu
              for (; element_courant < rang_sommet; element_courant++)
                {
                  const long sommet = elements[element_courant];
                  new_elements.append_array(sommet);
                  // On marque du sommet car il est utilise
                  // dans un espace distant
                  sommets_utilises[sommet]++;
                }
              // On passe le sommet a retirer
              element_courant++;
            }
          while (element_courant < nb_elements);
          // Mise a jour des elements distants
          espace_distant.set_elements(pe, new_elements);
        }
      espace_distant.calcul_liste_pe_voisins();
    }
    comm.end_comm();
  }

  // Suppression des sommets reels inutilises
  // (dans notre cas ceux sont des sommets sur des faces de frontieres ouvertes):
  // Mise a jour de :
  //  * sommets_
  //  * sommet_elem_
  //  * sommet_face_bord_
  //  * sommet_PE_owner_
  //  * drapeaux_sommets_
  // On construit en meme temps le tableau renum_sommets:
  // si i est l'indice actuel d'un sommet, renum_sommets[i]
  // est l'indice du sommet apres suppression des sommets
  // inutilises.
  ArrOfIntFT renum_sommets(sommets_.dimension(0));
  {
    const long nbsommets = sommets_.dimension(0);
    // Compteur de sommets apres suppression
    long n = 0;
    for (long i = 0; i < nbsommets; i++)
      {
        if (sommets_utilises[i])
          {
            renum_sommets[i] = n;
            sommets_(n, 0) = sommets_(i, 0);
            sommets_(n, 1) = sommets_(i, 1);
            if (dimension3)
              sommets_(n, 2) = sommets_(i, 2);
            sommet_elem_[n] = sommet_elem_[i];
            for (long dim=0; dim<dimension; dim++) sommet_face_(n,dim) = sommet_face_(i,dim); // EB
            sommet_face_bord_[n] = sommet_face_bord_[i];
            sommet_PE_owner_[n] = sommet_PE_owner_[i];
            drapeaux_sommets_[n] = drapeaux_sommets_[i];
            n++;
          }
        else
          {
            renum_sommets[i] = -1; // Le sommet i est supprime
          }
      }
    sommets_         .resize(n, Objet_U::dimension);
    sommet_elem_     .resize_array(n);
    sommet_face_     .resize(n, Objet_U::dimension); // EB
    sommet_face_bord_.resize_array(n);
    sommet_PE_owner_ .resize_array(n);
    drapeaux_sommets_.resize_array(n);

    equation_transport().nettoyer_proprietes_particules(sommets_utilises);

  }

  //On supprime ce qui est lie aux espaces distants et virtuel

  // Renumerotation des sommets dans desc_sommets
  {
    // Boucle sur les deux espaces : distant et virtuel
    for (long num_espace = 0; num_espace < 2; num_espace++)
      {
        // Espace_dv est soit l'ESPACE_D(istant), soit l'ESPACE_V(irtuel).
        Descripteur_FT& espace_dv =
          (num_espace == 0)
          ? desc_sommets_.espace_distant()
          : desc_sommets_.espace_virtuel();

        const ArrOfInt& pe_voisins = espace_dv.pe_voisins();
        const long nb_pe_voisins = pe_voisins.size_array();
        for (long indice_pe = 0; indice_pe < nb_pe_voisins; indice_pe++)
          {
            const long pe = pe_voisins[indice_pe];
            const ArrOfInt& elements = espace_dv.elements(pe);
            const long nb_elements = elements.size_array();
            new_elements.resize_array(nb_elements);
            for (long i = 0; i < nb_elements; i++)
              {
                const long num_sommet = elements[i];
                const long new_num = renum_sommets[num_sommet];
                assert(new_num >= 0);
                new_elements[i] = new_num;
              }
            espace_dv.set_elements(pe, new_elements);
          }
        // Mise a jour du descripteur:
        espace_dv.calcul_liste_pe_voisins();
      }
    desc_sommets_.calcul_schema_comm(sommets_.dimension(0));
  }

  // Mise a jour de sommet_num_owner_:
  {
    const long nbsommets = sommets_.dimension(0);
    sommet_num_owner_.resize_array(nbsommets);
    for (long i = 0; i < nbsommets; i++)
      sommet_num_owner_[i] = i;
    desc_sommets_.echange_espace_virtuel(sommet_num_owner_);
  }

  maillage_modifie(MINIMAL);
  if (Comm_Group::check_enabled()) check_mesh(1,1,1);
}

#if (defined(PATCH_HYSTERESIS_V2) || defined(PATCH_HYSTERESIS_V3))
// From Remaillage_FT : FTd_vecteur3 to ArrOfDouble...
inline double produit_scalaire(const ArrOfDouble& a, const ArrOfDouble& b)
{
  assert(a.size_array()==3);
  assert(b.size_array()==3);
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void produit_vectoriel(const ArrOfDouble& a, const ArrOfDouble& b, ArrOfDouble& resu)
{
  assert(a.size_array()==3);
  assert(b.size_array()==3);
  assert(resu.size_array()==3);
  resu[0] = a[1]*b[2] - a[2]*b[1];
  resu[1] = a[2]*b[0] - a[0]*b[2];
  resu[2] = a[0]*b[1] - a[1]*b[0];
}

void Maillage_FT_Disc::pre_lissage_courbure(ArrOfDouble& store_courbure_sommets, const long nitera) const
{
  // Un peu de lissage :
  const long nsom = nb_sommets();
  const ArrOfDouble& surface = get_update_surface_facettes();
  ArrOfDouble tmp_courbure_sommets(store_courbure_sommets);
  ArrOfDouble surface_autour_sommets(store_courbure_sommets);

  long i_facette;
  for (long iter = 0; iter < nitera; iter++)
    {
      tmp_courbure_sommets = store_courbure_sommets;
      store_courbure_sommets = 0.;
      surface_autour_sommets = 0.;
      for (i_facette = 0; i_facette < nb_facettes(); i_facette++)
        {
          // Ne pas calculer de flux pour les facettes virtuelles:
          if (facette_virtuelle(i_facette))
            continue;

          //const long s0 = facettes_(i_facette, 0);
          double moy = 0.;
          long count =0;
          for (long i = 0; i < 3; i++)
            {
              const long si = facettes_(i_facette, i);
              if (!sommet_ligne_contact(si))
                {
                  moy += tmp_courbure_sommets[si];
                  count++;
                }
            }
          if (count>0)
            moy /= count;
          for (long i = 0; i < 3; i++)
            {
              const long s = facettes_(i_facette, i);
              const double surf = surface(i_facette);
              store_courbure_sommets[s] += moy*surf;
              surface_autour_sommets[s] +=surf;
            }
        }
      for (long isom = 0; isom < nsom; isom++)
        if (surface_autour_sommets[isom]>0)
          store_courbure_sommets[isom] /=surface_autour_sommets[isom];
    }
}

// Correction a la contact line pour prendre en compte que la tangente au cercle n'est pas rigoureusement donnee par la facette
// Il faut faire une petite correction  (d'angle eps) pour evaluer l'angle de contact au sommet...
// Pour evaluer cette correction, on a besoin de connaitre la courbure normale kappa_n=1/r...
// Il faut donc une methode iterative (a 2 passe au moins).
//    1ere pass : On ne connait pas la courbure, on ne fait pas de correction.
//    2eme pass (et potentiellement les suivantes) : On utilise l'approximation de la courbure
//                                                   locale issue de la passe precedente.
void Maillage_FT_Disc::correction_costheta(const double c, const long s0, const long facette,
                                           /* const ArrOfDouble& s0s1, const ArrOfDouble& s0s2, */
                                           double ps) const
{
  // Soit l=norme(s0G) ou G est le cdg de la facette.
  // s0G = 1./3. * (s0s1+s0s2)
  double l2=0, r=0;
  /*  for (long j = 0; j < Objet_U::dimension; j++)
      {
        double tmp= (s0s1[j]+s0s2[j]);
        l2 += tmp*tmp;
      }
    l2 *= 1./3.;
  */
  // Un estimation de la taille de la facette est donnee par sa surface :
  //const ArrOfDouble& surface = get_update_surface_facettes();
  //l2 = surface(facette);
#if 1
  long ii=0;
  for (ii=0; ii<3; ii++)
    {
      if (facettes_(facette,ii) == s0)
        {
          break;
        }

    }
  assert(facettes_(facette,ii)==s0);
  const long s1 = facettes_(facette,(ii+1)%3);
  const long s2 = facettes_(facette,(ii+2)%3);

  ArrOfDouble som0(3),som1(3),som2(3), s0s1(3), s0s2(3);
  for (long j=0; j<Objet_U::dimension; j++)
    {
      som0[j] = sommets_(s0,j);
      som1[j] = sommets_(s1,j);
      som2[j] = sommets_(s2,j);
      s0s1[j] = som1[j] - som0[j];
      s0s2[j] = som2[j] - som0[j];
    }

  /*
   double tmp1=0,tmp2=0;
  for (long j = 0; j < Objet_U::dimension; j++)
    {
      tmp1 += (s0s1[j]*s0s1[j]);
      tmp2 += (s0s2[j]*s0s2[j]);
    }
  l2 = 0.5*(tmp1+tmp2);
  */
  for (long j = 0; j < Objet_U::dimension; j++)
    {
      double tmp= (s0s1[j]+s0s2[j]);
      l2 += tmp*tmp;
    }
  l2 *= 1./3.;
#endif

  // Pour r, on utilise la courbure... Donc iteratif...
  if (c!=0.)
    {
      r = correction_contact_courbure_coeff_/c;
      //r = 2./c; // Hypothese : localement isotrope : kappa_n=kappa_t=0.5kappa
      //r = 0.;   // Hypothese : kappa_n=0
      //r = 1./c; // Hypothese : courbure uniquement normale : kappa_n=kappa
      double r2 = r*r;
      if (l2>=r2)
        {
          Cerr << "Probleme lors du calcul de la courbure a la ligne de contact. "<< finl;
          Cerr << "Il semblerait que la courbure au sommet " << s0 << " sur le process "
               << Process::me() << " soit tres grande (ie un rayon de courbure plus petit"
               << " que la longueur caracteristique de l'element du Front. " << finl;
          Cerr << "Une erreur possible est un angle de contact initialement trop loin de la "
               << "gamme d'angles de contacts fournis dans le jeu de donnees. " << finl;
          Cerr << "Verifiez votre JDD puis contactez TRUST support." << finl;
          Process::exit();
        }
      //assert(l2<r2);
      // theta_reel = theta_apparent + eps // x=theta_apparent; eps+alpha= pi/2
      //            alpha est l'angle entre le plan de la facette et la normale au cercle (vers l'interieur).
      // cos(x+eps) = cos(x)*cos(eps)-sin(x)sin(eps)
      // avec :
      //  o cos(x) = ps
      //  o sin(x) = sqrt(1-ps^2)
      //  o cos(eps) = sin(alpha)=sin(pi/2-eps) = sqrt((r^2-l^2)/r^2)
      //  o sin(eps) = cos(alpha)=cos(pi/2-eps)=l/r
      ps = ps*sqrt((r2-l2)/r2) - sqrt((1.-ps*ps)*l2/r2);
    }
}

// Methode permettant de calculer lors d'une hysterisis l'angle de contact
// dans la plage des angles autorises qui soit le plus proche de l'angle de
// contact reel. C'est la premiere fois qu'on a besoin d'un angle de contact
// reelement mesure sur le maillage (et pas fourni en CL)
double Maillage_FT_Disc::calculer_costheta_objectif(const long s0, const long facette, const long call, const double c,
                                                    const DoubleTabFT& tab_cos_theta, ArrOfBit& drapeau_angle_in_range) const
{
  const Parcours_interface& parcours = refparcours_interface_.valeur();
  const DoubleTab& normale = get_update_normale_facettes();
  const long numface0 = sommet_face_bord_[s0];
  long j; /*
    ArrOfDouble som0(3),som1(3),som2(3), s0s1(3), s0s2(3);
    for (j=0; j<Objet_U::dimension; j++)
      {
        som0[j] = sommets_(s0,j);
        som1[j] = sommets_(s1,j);
        som2[j] = sommets_(s2,j);
        s0s1[j] = som1[j] - som0[j];
        s0s2[j] = som2[j] - som0[j];
      } */
  ArrOfDouble nprime(3), vect_base(3), ntheta(3);
  // Si on prend ntheta egale a nfacette
  // 				   au lieu d'essayer de determiner le vrai ntheta (qui correspond a l'objectif fixee)
  // ntheta[0] = normale(facette, 0) ;
  // ntheta[1] = normale(facette, 1) ;
  // ntheta[2] = normale(facette, 2) ;
  //                   Avec cette approximation, on ne tient pas compte de la CL et
  // 				   C'est comme si la courbure (dans la direction normale au bord) etait nulle!
  //      			   Dans le cas d'un cylindre comme la courbure tangeante est nulle, on obtient
  //                   Une courbure totale evaluee a 0 et donc un equilibre impossible quel que soit la CL...
  //
  // Ce que l'on veut plutot faire, c'est que ntheta forme un angle theta avec la normale au bord.
  //                   (ou theta doit etre calcule pour etre la projection de la valeur theta reele sur
  //                    l'intervalle d'angles authorises)
  //
  // Pour calculer la contribution en s0, on a besoin de connaitre l'angle de contact en s0
  // Au lieu de l'angle reel, on veut connaitre l'angle objectif, cad celui qui
  // doit etre atteint (si on est hors equilibbre) pour satisfaire la ligne de contact.

  // Pour cela, il faut ici calculer l'angle de contact reel a la paroi et soit :
  // 1. il est dans la gamme d'angles de contact souhaitee : C'est alors cet angle qu'on utilise
  //        pour le calcul de la courbure.
  // 2. il n'est pas dans la gamme des angles authorises : Il faut alors rechercher la borne de
  // l'intervalle la plus proche de l'angle reel et utiliser celle-ci pour calculer la courbure.
  // Dans ce second cas, cela aura pour effet de creer une courbure locale (aux sommets de la ligne de contact)
  // qui n'est pas homogene a celle a l'interieure du dom, ce qui va creer un terme source a divergence non nulle
  // pour la qdm, qui, via navier stokes va mettre le fluide en mouvement, qui, une fois la vitesse interpolee,
  // va deplacer le noeud en paroi dans la bonne direction.
  // C'est donc dans le cas 2 une methode tres indirecte d'imposer l'angle de contact.
  // Mais qui marche? Voir fiche hysterisis...

  const double costheta0 = tab_cos_theta(s0, 0);
  const double costheta1 = tab_cos_theta(s0, 1);
  /*  if (std::fabs(costheta0-costheta1)<1e-5) {
  	// Pas besoin de se casser la tete s'il n'y a pas d'hysterisis...
      return costheta0;
      }
  */

  // Normale unitaire a la face de bord (vers l'interieur)
  double nf0[3] = {0., 0., 0.};
  parcours.calculer_normale_face_bord(numface0,
                                      sommets_(s0,0), sommets_(s0,1),  sommets_(s0,2),
                                      nf0[0], nf0[1], nf0[2]);

  // Calcul du produit scalaire entre la normale au bord et la normale unitaire a la facette
  // On a besoin ici de celle orientee vers la vapeur pour etre coherent avec la definition de
  // l'angle de contact (qui est pris du cote liquide). Donc il faut -normale dans le ps :
  double ps = 0;
  for (j = 0; j < Objet_U::dimension; j++)
    ps -= nf0[j]* normale(facette, j);
  // Le produit scalaire entre les deux normales, c'est comme celui entre les 2 tangentes...
  // On a donc ps=cos(theta)

  // Correction de l'angle de contact (prise en compte de la tangente au cercle
  // qui n'est pas rigoureusement donnee par la facette
  if ((call> 1) && (correction_contact_courbure_coeff_!=0.))
    correction_costheta(c, s0, facette, ps);

  double costheta=0.;
  if ((ps-costheta0)*(ps-costheta1) <=0 )
    {
      // L'angle reel est dans l'intervalle authorise. On l'utilise pour le calcul de la courbure.
      costheta = ps;
    }
  else
    {
      drapeau_angle_in_range.clearbit(s0);
      // L'angle reel est en dehors de l'intervalle. On retient donc la borne de l'intervalle la
      // plus proche comme valeur pour utiliser dans le calcul de la courbure.
      if (std::fabs(ps-costheta0) <std::fabs(ps-costheta1))
        {
          // On est proche de la borne costheta0. C'est elle qu'on va utiliser.
          costheta = costheta0;
        }
      else
        {
          costheta = costheta1;
        }
    }
  return costheta;
}
#endif

#ifdef PATCH_HYSTERESIS_V2
// A: Le point a reflechir.
// nprime: une normale au plan quelconque! -> devient unitaire en sortie!
// Ar: La reflexion de A
static void miroir(const ArrOfDouble& A, ArrOfDouble& nprime,const ArrOfDouble& O, ArrOfDouble& Ar)
{
  // Rendre la normale unitaire :
  const double l = sqrt(nprime[0] * nprime[0] + nprime[1] * nprime[1] + nprime[2] * nprime[2]);
  if (l != 0.)
    {
      double inv_l = 1. / l;
      nprime[0] *= inv_l;
      nprime[1] *= inv_l;
      nprime[2] *= inv_l;
    }

  const long m=A.size_array();
  assert(m==3);
  ArrOfDouble OA(A);
  OA-=O;// OA=A-O;
  double ps=produit_scalaire(OA,nprime);
  //OAr=Ar-O=OA-2.*ps*n => Ar=A-2.*ps*n
  for (long i=0; i<m; i++)
    {
      Ar[i] =A[i]-2.*ps*nprime[i];
    }
}

static void normalize(ArrOfDouble& nprime)
{
  // Rendre la normale unitaire :
  const double l = sqrt(nprime[0] * nprime[0] + nprime[1] * nprime[1] + nprime[2] * nprime[2]);
  if (l != 0.)
    {
      double inv_l = 1. / l;
      nprime[0] *= inv_l;
      nprime[1] *= inv_l;
      nprime[2] *= inv_l;
    }
}
#endif

/*! @brief Calcul de la courbure discrete du maillage aux sommets.
 *
 * Methode de calcul : voir these B. Mathieu paragraphe 3.3.3 page 97
 *   La courbure est egale a la differentielle de la surface d'interface par rapport
 *    au deplacement des sommets, divisee par la differentielle du volume.
 *
 *
 * @param (courbure_sommets) Tableau dans lequel on veut stocker la courbure aux sommets. La valeur initiale du tableau est perdue. L'espace virtuel du tableau resultat est a jour.
 */
void Maillage_FT_Disc::calcul_courbure_sommets(ArrOfDouble& courbure_sommets, const long call) const
{
  const DoubleTab& normale = get_update_normale_facettes();
  const ArrOfDouble& surface = get_update_surface_facettes();

  const long nsom = nb_sommets();
  const long dim = Objet_U::dimension;
  const long nfaces = nb_facettes();

#if (defined(PATCH_HYSTERESIS_V2) || defined(PATCH_HYSTERESIS_V3))
  ArrOfBit drapeau_angle_in_range;
  drapeau_angle_in_range.resize_array(nsom);
  drapeau_angle_in_range=1; // par defaut on presume que oui. Si pour une seule facette on est en dehors,
  //                                                          on le met a 0 pour le sommet en question.
  // coef default value is 2.
  Cerr << "Your choice of parameters : correction coefficient " << correction_contact_courbure_coeff_
       << " and iter pre-lissage : " << niter_pre_lissage_ << finl;
  ArrOfDouble store_courbure_sommets(courbure_sommets);
  if (call> 1)
    {
      if (courbure_sommets.size_array() != nsom)
        {
          Cerr << "Erreur dans la seconde passe du calcul de la courbure... Dimensions des tableaux differentes!" << finl;
          Process::exit();
        }
      pre_lissage_courbure(store_courbure_sommets, niter_pre_lissage_);
    }
  else
    {
      // Lors du premier passage, il faut dimensionner correctement le tableau de courbure...
      courbure_sommets.resize_array(nsom);
      store_courbure_sommets.resize_array(nsom); // aaa
    }
#else
  courbure_sommets.resize_array(nsom);
#endif

  // Differentielle de la surface d'interface par rapport au deplacement de chaque sommet
  DoubleTab d_surface(nsom, dim);
  // Differentielle du volume de la phase 1 par rapport au deplacement de chaque sommet
  DoubleTab d_volume(nsom, dim);

  double n_s[3] = {0.,0.,0.};
  const double inverse_dimension = 1. / (double) dimension;
  const double un_tiers = 1. / 3.;
  const double un_sixieme = 1. / 6.;

  // This is based on the angle given in the data file (that is the micros/Young contact angle)
  DoubleTabFT tab_cos_theta;
  calculer_costheta_minmax(tab_cos_theta);

  long i, j, facette;

  // Cette classe sert pour promener les sommets sur le bord du domaine.
  const Parcours_interface& parcours = refparcours_interface_.valeur();
#if TCL_MODEL
  long flag_tcl = 0;
  double l_v = 0.;
  // interfacial velocity
  DoubleTab vit(nsom, dim);
  if (refequation_transport_.non_nul())
    {
      const Transport_Interfaces_FT_Disc& eq_interfaces = refequation_transport_.valeur();
      const Probleme_base& pb = eq_interfaces.get_probleme_base();
      Probleme_FT_Disc_gen& pb_ft = ref_cast_non_const(Probleme_FT_Disc_gen, pb);
      Triple_Line_Model_FT_Disc& tcl = pb_ft.tcl();
      if (tcl.is_activated() && tcl.is_capillary_activated())
        {
          flag_tcl = 1;
          l_v = tcl.get_lv();
          Postraitement_base::Localisation loc = Postraitement_base::SOMMETS;
          Motcle nom_du_champ = "vitesse";
          Cerr << "Validation and checking required in Maillage_FT_Disc::calcul_courbure_sommets" << finl;
          Process::exit();
          // Useless init to 0:
          // for (long ii=0; ii< nsom; ii++)
          //  for (long jj=0; jj< dim; jj++)
          //    vit(ii,jj) = 0.;
          eq_interfaces.get_champ_post_FT(nom_du_champ, loc, &vit); // HACK !!!! (warning, try debug to make sure it works if you want to remove it!!)

          //  const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF, domaine_cl);
          // It relies on the classical assumption in the FT module that the first equation is NS.
        }
      //Cerr <<"TCL Eq. 0 is " <<  equation_transport().get_probleme_base().equation(0).que_suis_je() << finl;
    }
#endif

  for (facette = 0; facette < nfaces; facette++)
    {
      // Si la facette est reelle, on ajoute la contribution a la differentielle
      // des deux ou trois sommets
      if (! facette_virtuelle(facette))
        {
          if (!bidim_axi)
            {
              long ii;
              const double surface_sur_dim = surface[facette] * inverse_dimension;
              for (ii = 0; ii < dim; ii++)
                {
                  // La differentielle de volume par rapport a un deplacement v est
                  // (v scalaire n) * surface_facette / dimension
                  // avec n la normale a la facette dirigee vers la phase 1.
                  n_s[ii] = normale(facette, ii) * surface_sur_dim;
                }
              for (ii = 0; ii < dim; ii++)
                {
                  const long sommet = facettes_(facette, ii);
                  for (j = 0; j < dim; j++)
                    d_volume(sommet, j) += n_s[j];
                }
              // Calcul de la differentielle de surface :
              if (dim == 2)
                {
                  long som[2];
                  som[0] = facettes_(facette, 0);
                  som[1] = facettes_(facette, 1);
                  double n[2];
                  n[0] = normale(facette, 0);
                  n[1] = normale(facette, 1);
                  // La differentielle est orthogonale a la normale :
                  d_surface(som[0], 0) -= n[1];
                  d_surface(som[0], 1) += n[0];
                  d_surface(som[1], 0) += n[1];
                  d_surface(som[1], 1) -= n[0];

                  // Traitement des lignes de contact: ajout d'une contribution
                  // cos(theta) * differentielle de la surface de bord mouillee par
                  // la phase 0
                  for (long i2 = 0; i2 < 2; i2++)
                    {
                      const long isom2 = som[i2];
                      const long face = sommet_face_bord_[isom2];
                      if (face < 0) // pas une face de bord
                        continue;

#ifndef PATCH_HYSTERESIS_V2
                      const double costheta = tab_cos_theta(isom2, 0);
#else
                      const double costheta0 = tab_cos_theta(som[i2], 0);
                      const double costheta1 = tab_cos_theta(som[i2], 1);
                      // Normale unitaire a la face de bord (vers l'interieur)
                      double nf1[3] = {0., 0., 0.};
                      parcours.calculer_normale_face_bord(face, som[i2], som[i2], 0,
                                                          nf1[0], nf1[1], nf1[2]);

                      // Calcul du produit scalaire entre la normale unitaire a la facette
                      double ps = 0;
                      for (j = 0; j < dim; j++)
                        {
                          ps += nf1[j]* n[j];
                        }
                      // Le produit scalaire entre les deux normales, c'est comme celui entre les 2 tangentes...
                      // On a donc cos(theta)
                      // A priori, si l'hysterisis fonctionne bien, il est dans l'intervalle
                      Nom st;
                      st = (ps-costheta0)*(ps-costheta1) <=0 ? "YES":"NO";
                      Cerr << "GB_CALCUL_COURBURE 2D PLAN : ps "
                           << ps << " in costheta ["<< costheta0 << " ; " << costheta1 << "]? " << st << finl;
                      const double costheta = ps;
                      Cerr << "Quelle est le sommet a choisir? Est-ce que i2 " << i2 << "est bien le sommet du bord?" << finl;
                      Cerr << "GB_CALCUL_COURBURE 2D PLAN: "
                           << "A valider depuis hysterisis. Le developpement fait en 3D n'a pas ete applique ici!" << finl;
                      Process::exit();
#endif

                      // Normale unitaire au bord
                      double nx, ny, nz;
                      parcours.calculer_normale_face_bord(face,
                                                          sommets_(som[i2],0), sommets_(som[i2],1), 0.,
                                                          nx, ny, nz);
                      // Ajout de la contribution de la surface:
                      double signe = (i2==0) ? -1. : 1.;
                      d_surface(som[i2], 0) += signe * ny * costheta;
                      d_surface(som[i2], 1) -= signe * nx * costheta;
                    }

                }
              else
                {
                  long sommets_loc[3];
                  sommets_loc[0] = facettes_(facette, 0);
                  sommets_loc[1] = facettes_(facette, 1);
                  sommets_loc[2] = facettes_(facette, 2);
                  ArrOfDouble n(3);// GB mod double[3] to ArrOfDouble
                  n[0] = normale(facette, 0) * 0.5;
                  n[1] = normale(facette, 1) * 0.5;
                  n[2] = normale(facette, 2) * 0.5;
                  for (ii = 0; ii < 3; ii++)
                    {
                      // La differentielle de surface pour un deplacement du sommet i
                      // est le produit vectoriel de la normale par le vecteur
                      // s2s1 = (sommet[(i+1)%3] - sommet[(i+2)%3]) * 0.5
                      // (vecteur de norme "base du triangle * 0.5" et de direction
                      //  la hauteur du triangle)
                      const long s0 = sommets_loc[ii];
                      const long s1 = sommets_loc[ (ii+1)%3 ];
                      const long s2 = sommets_loc[ (ii+2)%3 ];
                      ArrOfDouble s2s1(3); // GB mod double[3] to ArrOfDouble
                      s2s1[0] = sommets_(s1,0) - sommets_(s2,0);
                      s2s1[1] = sommets_(s1,1) - sommets_(s2,1);
                      s2s1[2] = sommets_(s1,2) - sommets_(s2,2);
                      // On calcul la courbure partout de la meme maniere
                      // meme si on est sur une ligne de contact (ie dans un premier temps, sans tenir compte de
                      // l'effet de la ligne de contact). On corrige l'effet de la ligne de contact par la suite...
                      d_surface(s0,0) += s2s1[1] * n[2] - s2s1[2] * n[1];
                      d_surface(s0,1) += s2s1[2] * n[0] - s2s1[0] * n[2];
                      d_surface(s0,2) += s2s1[0] * n[1] - s2s1[1] * n[0];
#ifdef PATCH_HYSTERESIS_V2
                      if (methode_calcul_courbure_contact_line_ == MIRROR)
                        {
                          const long numface0 = sommet_face_bord_[s0];
                          if (numface0>=0)
                            {
                              // Le sommet s0 est sur une ligne de contact.

                              // Voila la suite :
                              ArrOfDouble som0(3),som1(3),som2(3), s0s1(3), s0s2(3);
                              for (j=0; j<dim; j++)
                                {
                                  som0[j] = sommets_(s0,j);
                                  som1[j] = sommets_(s1,j);
                                  som2[j] = sommets_(s2,j);
                                  s0s1[j] = som1[j] - som0[j];
                                  s0s2[j] = som2[j] - som0[j];
                                }
                              ArrOfDouble nprime(3), vect_base(3), ntheta(3);
                              // Si on prend ntheta egale a nfacette
                              // 				   au lieu d'essayer de determiner le vrai ntheta (qui correspond a l'objectif fixee)
                              // ntheta[0] = normale(facette, 0) ;
                              // ntheta[1] = normale(facette, 1) ;
                              // ntheta[2] = normale(facette, 2) ;
                              //                   Avec cette approximation, on ne tient pas compte de la CL et
                              // 				   C'est comme si la courbure (dans la direction normale au bord) etait nulle!
                              //      			   Dans le cas d'un cylindre comme la courbure tangeante est nulle, on obtient
                              //                   Une courbure totale evaluee a 0 et donc un equilibre impossible quel que soit la CL...
                              //
                              // Ce que l'on veut plutot faire, c'est que ntheta forme un angle theta avec la normale au bord.
                              //                   (ou theta doit etre calcule pour etre la projection de la valeur theta reele sur
                              //                    l'intervalle d'angles authorises)
                              //
                              // Pour calculer la contribution en s0, on a besoin de connaitre l'angle de contact en s0
                              // Au lieu de l'angle reel, on veut connaitre l'angle objectif, cad celui qui
                              // doit etre atteint (si on est hors equilibbre) pour satisfaire la ligne de contact.

                              // Pour cela, il faut ici calculer l'angle de contact reel a la paroi et soit :
                              // 1. il est dans la gamme d'angles de contact souhaitee : C'est alors cet angle qu'on utilise
                              //        pour le calcul de la courbure.
                              // 2. il n'est pas dans la gamme des angles authorises : Il faut alors rechercher la borne de
                              // l'intervalle la plus proche de l'angle reel et utiliser celle-ci pour calculer la courbure.
                              // Dans ce second cas, cela aura pour effet de creer une courbure locale (aux sommets de la ligne de contact)
                              // qui n'est pas homogene a celle a l'interieure du dom, ce qui va creer un terme source a divergence non nulle
                              // pour la qdm, qui, via navier stokes va mettre le fluide en mouvement, qui, une fois la vitesse interpolee,
                              // va deplacer le noeud en paroi dans la bonne direction.
                              // C'est donc dans le cas 2 une methode tres indirecte d'imposer l'angle de contact.
                              // Mais qui marche? Voir fiche hysterisis...

                              const double costheta0 = tab_cos_theta(s0, 0);
                              const double costheta1 = tab_cos_theta(s0, 1);

                              // Normale unitaire a la face de bord (vers l'interieur)
                              double nf0[3] = {0., 0., 0.};
                              parcours.calculer_normale_face_bord(numface0,
                                                                  sommets_(s0,0), sommets_(s0,1),  sommets_(s0,2),
                                                                  nf0[0], nf0[1], nf0[2]);

                              // Calcul du produit scalaire entre la normale au bord et la normale unitaire a la facette
                              // On a besoin ici de celle orientee vers la vapeur pour etre coherent avec la definition de
                              // l'angle de contact (qui est pris du cote liquide). Donc il faut -normale dans le ps :
                              double ps = 0;
                              for (j = 0; j < dim; j++)
                                ps -= nf0[j]* normale(facette, j);
                              // Le produit scalaire entre les deux normales, c'est comme celui entre les 2 tangentes...
                              // On a donc ps=cos(theta)

                              // Correction de l'angle de contact (prise en compte de la tangente au cercle
                              // qui n'est pas rigoureusement donnee par la facette
                              if ((call> 1) && (correction_contact_courbure_coeff_!=0.))
                                //correction_costheta(store_courbure_sommets[s0], s0, s0s1, s0s2, ps);
                                correction_costheta(store_courbure_sommets[s0], s0, facette, ps);

                              double costheta=0.;
                              if ((ps-costheta0)*(ps-costheta1) <=0 )
                                {
                                  // L'angle reel est dans l'intervalle authorise. On l'utilise pour le calcul de la courbure.
                                  costheta = ps;
                                }
                              else
                                {
                                  // L'angle reel est en dehors de l'intervalle. On retient donc la borne de l'intervalle la
                                  // plus proche comme valeur pour utiliser dans le calcul de la courbure.
                                  if (std::fabs(ps-costheta0) <std::fabs(ps-costheta1))
                                    {
                                      // On est proche de la borne costheta0. C'est elle qu'on va utiliser.
                                      costheta = costheta0;
                                    }
                                  else
                                    {
                                      costheta = costheta1;
                                    }
                                }
                              {
                                // Pour remplacer une bonne part du bloc precedent :
                                const double costheta_bis = calculer_costheta_objectif(s0, facette, call, store_courbure_sommets[s0],
                                                                                       tab_cos_theta, drapeau_angle_in_range);
                                if (std::fabs(costheta-costheta_bis)>1e-8)
                                  {
                                    Cerr << "Oh oh... les 2 methodes different... " << finl;
                                    Process::exit();
                                  }
                              }

                              // On cherche quel type de contact on a (ponctuel ou lineique).
                              const long numface1 = sommet_face_bord_[s1];
                              const long numface2 = sommet_face_bord_[s2];

#ifdef DEBUG_HYSTERESIS_V2
                              const long inode = 6 ;
//                          if (facette == 11)
                              if ( (ii==0) && ((s0 == inode) || (s1 == inode) || (s2 == inode))
                                   && (!(numface1 >= 0 && numface2 >= 0)))  /* pour ne pas afficher les facettes avec 3 som bords*/
                                {
                                  Cerr << "TAG FACETTE "<< facette << " : " << s0 << " " << s1 << " " << s2  << finl;
                                  Cerr << som0 << finl;
                                  Cerr << som1 << finl;
                                  Cerr << som2 << finl;
                                  Cerr << som0 << finl;
                                }
#endif

                              if (numface1 < 0 && numface2 < 0)
                                {
                                  // Le contact est ponctuel au sommet s0

                                  // Construction du vecteur tb :
                                  ArrOfDouble v(3), nb(3), tb(3);
                                  // Pour cela, on a besoin d'un segment de la facette (v)
                                  for (j=0; j<dim; j++)
                                    v[j] = som2[j] - som1[j]; // s1s2
                                  // et de la normale au bord (nb):
                                  parcours.calculer_normale_face_bord(numface0,
                                                                      sommets_(s0,0), sommets_(s0,1),  sommets_(s0,2),
                                                                      nb[0], nb[1], nb[2]);
                                  // pour en deduire le vecteur du plan de la face de bord normale a la facette :
                                  produit_vectoriel(v,nb,tb); // tb non-unitaire
                                  normalize(tb);
                                  // D'apres la rotation, le vecteur ntheta se construit comme :
                                  // (attention, ici, theta est bien l'angle objectif souhaitee en accord avec la CL),
                                  //  Hors equilibre, c'est different de l'angle de contact reel.
                                  const double sintheta = sqrt(1.-costheta*costheta);
                                  for (j=0; j<dim; j++)
                                    ntheta[j] = sintheta*tb[j] + costheta*nb[j];

                                  // on cree un miroir de s1 et de s2 par rapport au plan parallel a ntheta et s2s1 et contenant s0.
                                  ArrOfDouble nplan(3); // la normale au plan.

                                  // Seul la direction de nplan compte. Pas son signe ou sa valeur.
                                  parcours.projeter_vecteur_sur_face(numface0,v[0], v[1],v[2]);
                                  produit_vectoriel(ntheta,v,nplan);

                                  ArrOfDouble reflexion_som1(dim),reflexion_som2(dim);
                                  miroir(som1,nplan/* normale_au_plan*/,som0/*pt_du_plan*/, reflexion_som1);
                                  miroir(som2,nplan/* normale_au_plan*/,som0/*pt_du_plan*/, reflexion_som2);

                                  // remise de points dans l'ordre pour que le miroir tourne dans le meme sens...
                                  // Il faut donc utiliser les sommets dans l'ordre inverse pour construire
                                  // la normale de la partie reflechie...
                                  ArrOfDouble s0s2p(3), s0s1p(3), s1ps2p(3);
                                  for (j=0; j<dim; j++)
                                    {
                                      s0s2p[j] = reflexion_som2[j]-som0[j];
                                      s0s1p[j] = reflexion_som1[j]-som0[j];
                                      s1ps2p[j] = reflexion_som2[j]-reflexion_som1[j];
                                    }

                                  // Calcul de la nouvelle normale... (non-unitaire)
                                  produit_vectoriel(s0s2p,s0s1p,nprime);

                                  // Le vecteur de base pour calcul de la surface est :
                                  vect_base=s1ps2p;
#ifdef DEBUG_HYSTERESIS_V2
                                  if ( (s0 == inode) || (s1 == inode) || (s2 == inode) )
                                    {
                                      Cerr << "TAG REFLEXION s1s2 FACETTE "<< facette << " : " << s0 << " " << s1 << " " << s2  << finl;
                                      Cerr << som0 << finl;
                                      Cerr << reflexion_som1 << finl;
                                      Cerr << reflexion_som2 << finl;
                                      Cerr << som0 << finl;
                                    }
#endif
                                }
                              else if (numface1 >= 0 && numface2 >= 0
                                       && !((s0==s1) &&  (s0==s2)))
                                {
                                  // les 3 sommets sont sur des bords...
                                  //   et la facette n'est pas degeneree (ie 3 sommets identiques...)
                                  Cerr << "Cas de 3 sommets de la facette "  << facette
                                       << " (" << s0 << " " << s1 << " " << s2 << " )"
                                       << " sur process "
                                       << Process::me() << "  pas prevu pour l'instant!" << finl;
                                  Cerr << "som0 " << som0 << finl;
                                  Cerr << "som1 " << som1 << finl;
                                  Cerr << "som2 " << som2 << finl;
                                  Cerr << "Do nothing instead of exiting... " << finl;
                                  nprime *=0.;
                                  vect_base *=0.;
                                  //Cerr << "Exiting... " << finl;
                                  //Process::exit();
                                }
                              else if (numface1 >= 0)
                                {
                                  // Le segment s0s1 est une ligne de contact:

                                  // Construction du vecteur tb :
                                  ArrOfDouble v(3), nb(3), tb(3);
                                  // Pour cela, on a besoin d'un segment de la facette (v)
                                  for (j=0; j<dim; j++)
                                    v[j] = som0[j] - som1[j]; // s1s0 // Pour des questions d'orientation..
                                  // et de la normale au bord (nb):
                                  parcours.calculer_normale_face_bord(numface0,
                                                                      sommets_(s0,0), sommets_(s0,1),  sommets_(s0,2),
                                                                      nb[0], nb[1], nb[2]);
                                  {
                                    ArrOfDouble nb1(3);
                                    parcours.calculer_normale_face_bord(numface1,
                                                                        sommets_(s1,0), sommets_(s1,1),  sommets_(s1,2),
                                                                        nb1[0], nb1[1], nb1[2]);
                                    if (std::fabs(nb[0]*nb1[0]+nb[1]*nb1[1]+nb[2]*nb1[2]-1.) > 0.05)
                                      {
                                        nb+=nb1;
                                        normalize(nb);
                                        Cerr << "Mismatching normals... " << finl;
                                        Cerr << "nb: " << nb[0] << " " << nb[1] << " " << nb[2] << finl;
                                        Cerr << "nb1: " << nb1[0] << " " << nb1[1] << " " << nb1[2] << finl;
                                        Cerr << "Le probleme sera resolu quand on aura conserve les facettes de coins" << finl;
                                        Cerr << "Skipping for the time being instead of exiting"  << finl;
                                        //Process::exit();
                                      }

                                  }

                                  // pour en deduire le vecteur du plan de la face de bord normale a la facette :
                                  produit_vectoriel(v,nb,tb); // tb non-unitaire
                                  normalize(tb);
                                  // D'apres la rotation, le vecteur ntheta se construit comme :
                                  // (attention, ici, theta est bien l'angle objectif souhaitee en accord avec la CL,
                                  //  Hors equilibre, c'est different de l'angle de contact reel.
                                  const double sintheta = sqrt(1.-costheta*costheta);
                                  for (j=0; j<dim; j++)
                                    ntheta[j] = sintheta*tb[j] + costheta*nb[j];

                                  // on cree un miroir de s2 par rapport au plan parallel a ntheta et s0s1 et contenant s0.
                                  ArrOfDouble nplan(3); // la normale au plan.
                                  produit_vectoriel(s0s1,ntheta,nplan);
                                  ArrOfDouble reflexion_som2(dim);
                                  miroir(som2,nplan/* normale_au_plan*/,som0/*pt_du_plan*/, reflexion_som2);

                                  // remise de points dans l'ordre pour que le miroir tourne dans le meme sens...
                                  // Il faut donc utiliser les sommets dans l'ordre inverse pour la partie reflechie...
                                  ArrOfDouble s0s2p(3), s2ps1(3);
                                  for (j=0; j<dim; j++)
                                    {
                                      s0s2p[j] = reflexion_som2[j]-som0[j];
                                      s2ps1[j] = som1[j]-reflexion_som2[j];
                                    }

                                  // Calcul de la nouvelle normale...(non-unitaire)
                                  produit_vectoriel(s0s2p,s0s1,nprime);

                                  // Le vecteur de base pour calcul de la surface est :
                                  vect_base=s2ps1;
#ifdef DEBUG_HYSTERESIS_V2
                                  if ( (s0 == inode) || (s1 == inode) || (s2 == inode) )
                                    {
                                      Cerr << "TAG REFLEXION s0s1 "<< facette << " : " << s0 << " " << s1 << " " << s2  << finl;
                                      Cerr << som0 << finl;
                                      Cerr << som1 << finl;
                                      Cerr << reflexion_som2 << finl;
                                      Cerr << som0 << finl;
                                    }
#endif
                                }
                              else if (numface2 >= 0)
                                {
                                  // Le segment s0s2 est une ligne de contact:

                                  // Construction du vecteur tb :
                                  ArrOfDouble v(3), nb(3), tb(3);
                                  // Pour cela, on a besoin d'un segment de la facette (v)
                                  for (j=0; j<dim; j++)
                                    v[j] = som2[j] - som0[j]; // s0s2 // Pour des questions d'orientation..
                                  // et de la normale au bord (nb):
                                  parcours.calculer_normale_face_bord(numface0,
                                                                      sommets_(s0,0), sommets_(s0,1),  sommets_(s0,2),
                                                                      nb[0], nb[1], nb[2]);
                                  {
                                    ArrOfDouble nb2(3);
                                    parcours.calculer_normale_face_bord(numface2,
                                                                        sommets_(s1,0), sommets_(s1,1),  sommets_(s1,2),
                                                                        nb2[0], nb2[1], nb2[2]);
                                    if (std::fabs(nb[0]*nb2[0]+nb[1]*nb2[1]+nb[2]*nb2[2]-1.) > 0.05)
                                      {
                                        nb+=nb2;
                                        normalize(nb);
                                        Cerr << "Mismatching normals... " << finl;
                                        Cerr << "nb: " << nb[0] << " " << nb[1] << " " << nb[2] << finl;
                                        Cerr << "nb2: " << nb2[0] << " " << nb2[1] << " " << nb2[2] << finl;
                                        Cerr << "Le probleme sera resolu quand on aura conserve les facettes de coins" << finl;
                                        Cerr << "Skipping for the time being instead of exiting"  << finl;
                                        //Process::exit();
                                      }
                                  }

                                  // pour en deduire le vecteur du plan de la face de bord normale a la facette :
                                  produit_vectoriel(v,nb,tb); // tb non-unitaire
                                  normalize(tb);
                                  // D'apres la rotation, le vecteur ntheta se construit comme :
                                  // (attention, ici, theta est bien l'angle objectif souhaitee en accord avec la CL,
                                  //  Hors equilibre, c'est different de l'angle de contact reel.
                                  const double sintheta = sqrt(1.-costheta*costheta);
                                  for (j=0; j<dim; j++)
                                    ntheta[j] = sintheta*tb[j] + costheta*nb[j];

                                  // on cree un miroir de s1 par rapport au plan parallel a ntheta et s0s2 et contenant s0.
                                  ArrOfDouble nplan(3); // la normale au plan.
                                  produit_vectoriel(s0s2,ntheta,nplan);
                                  ArrOfDouble reflexion_som1(dim);
                                  miroir(som1,nplan/* normale_au_plan*/,som0/*pt_du_plan*/, reflexion_som1);

                                  // remise de points dans l'ordre pour que le miroir tourne dans le meme sens...
                                  // Il faut donc utiliser les sommets dans l'ordre inverse pour la partie reflechie...
                                  ArrOfDouble s0s1p(3), s2s1p(3);
                                  for (j=0; j<dim; j++)
                                    {
                                      s0s1p[j] = reflexion_som1[j]-som0[j];
                                      s2s1p[j] = reflexion_som1[j]-som2[j];
                                    }

                                  // Calcul de la nouvelle normale...(non-unitaire)
                                  produit_vectoriel(s0s2,s0s1,nprime);

                                  // Le vecteur de base pour calcul de la surface est :
                                  vect_base=s2s1p;
#ifdef DEBUG_HYSTERESIS_V2
                                  if ( (s0 == inode) || (s1 == inode) || (s2 == inode) )
                                    //if (facette == 11)
                                    {
                                      Cerr << "TAG REFLEXION s0s2 "<< facette << " : " << s0 << " " << s1 << " " << s2  << finl;
                                      Cerr << som0 << finl;
                                      Cerr << reflexion_som1 << finl;
                                      Cerr << som2 << finl;
                                      Cerr << som0 << finl;
                                    }
#endif
                                }
                              else
                                {
                                  Cerr << "Cas impossible... WTF!" << finl;
                                  Process::exit();
                                }

                              // Rendre la normale unitaire :
                              normalize(nprime);

                              //Calculer contribution du miroir a d_surface +=
                              d_surface(s0,0) += (vect_base[1] * nprime[2] - vect_base[2] * nprime[1])*0.5;
                              d_surface(s0,1) += (vect_base[2] * nprime[0] - vect_base[0] * nprime[2])*0.5;
                              d_surface(s0,2) += (vect_base[0] * nprime[1] - vect_base[1] * nprime[0])*0.5;

                              //Calculer contribution du miroir a d_volume +=
                              for (j = 0;
                                   j < dim;
                                   j++)
                                d_volume(s0, j) += nprime[j] * surface_sur_dim;
                            }
                        }
#endif
#if defined(PATCH_HYSTERESIS_V3)
                      if ((methode_calcul_courbure_contact_line_ == IMPROVED)
                          || (methode_calcul_courbure_contact_line_ == HYSTERESIS))
                        {
                          // Traitement des lignes de contact: ajout d'une contribution
                          // cos(theta) * differentielle de la surface de bord mouillee par
                          // la phase 0
                          const long numface1 = sommet_face_bord_[s1];
                          const long numface2 = sommet_face_bord_[s2];
                          // Le segment s1s2 est une ligne de contact:
                          if (numface1 >= 0 && numface2 >= 0)
                            {
                              // On prend l'angle de contact au milieu du segment,
                              // a partir des valeurs cibles de l'angle de contact a chaque sommet.
                              const double costheta_s1 = calculer_costheta_objectif(s1,facette,call,store_courbure_sommets[s1],
                                                                                    tab_cos_theta, drapeau_angle_in_range);
                              const double costheta_s2 = calculer_costheta_objectif(s2,facette,call,store_courbure_sommets[s2],
                                                                                    tab_cos_theta, drapeau_angle_in_range);
                              //Cerr << "Som0 " << s0 << " cos(theta) en s1 " << costheta_s1 << " et s2 " << costheta_s2 << finl;
                              const double costheta = (costheta_s1+costheta_s2) * 0.5;
                              // Normale unitaire au bord
                              double nx, ny, nz;
                              parcours.calculer_normale_face_bord(numface1,
                                                                  sommets_(s1,0), sommets_(s1,1),  sommets_(s1,2),
                                                                  nx, ny, nz);
                              // produit vectoriel s1s2 vectoriel n
                              d_surface(s1,0) -= (s2s1[1] * nz - s2s1[2] * ny) * costheta * 0.5;
                              d_surface(s1,1) -= (s2s1[2] * nx - s2s1[0] * nz) * costheta * 0.5;
                              d_surface(s1,2) -= (s2s1[0] * ny - s2s1[1] * nx) * costheta * 0.5;

                              parcours.calculer_normale_face_bord(numface2,
                                                                  sommets_(s2,0), sommets_(s2,1),  sommets_(s2,2),
                                                                  nx, ny, nz);
                              // produit vectoriel s1s2 vectoriel n
                              d_surface(s2,0) -= (s2s1[1] * nz - s2s1[2] * ny) * costheta * 0.5;
                              d_surface(s2,1) -= (s2s1[2] * nx - s2s1[0] * nz) * costheta * 0.5;
                              d_surface(s2,2) -= (s2s1[0] * ny - s2s1[1] * nx) * costheta * 0.5;
                            }
                        }
#endif
                      if (methode_calcul_courbure_contact_line_ == STANDARD)
                        {
                          // Traitement des lignes de contact: ajout d'une contribution
                          // cos(theta) * differentielle de la surface de bord mouillee par
                          // la phase 0
                          const long numface1 = sommet_face_bord_[s1];
                          const long numface2 = sommet_face_bord_[s2];
                          // Le segment s1s2 est une ligne de contact:
                          if (numface1 >= 0 && numface2 >= 0)
                            {
                              // On prend l'angle de contact au milieu du segment:
                              const double costheta = (tab_cos_theta(s1, 0) + tab_cos_theta(s2, 0)) * 0.5;
                              // Normale unitaire au bord
                              double nx, ny, nz;
                              parcours.calculer_normale_face_bord(numface1,
                                                                  sommets_(s1,0), sommets_(s1,1),  sommets_(s1,2),
                                                                  nx, ny, nz);
                              // produit vectoriel s1s2 vectoriel n
                              d_surface(s1,0) -= (s2s1[1] * nz - s2s1[2] * ny) * costheta * 0.5;
                              d_surface(s1,1) -= (s2s1[2] * nx - s2s1[0] * nz) * costheta * 0.5;
                              d_surface(s1,2) -= (s2s1[0] * ny - s2s1[1] * nx) * costheta * 0.5;

                              parcours.calculer_normale_face_bord(numface2,
                                                                  sommets_(s2,0), sommets_(s2,1),  sommets_(s2,2),
                                                                  nx, ny, nz);
                              // produit vectoriel s1s2 vectoriel n
                              d_surface(s2,0) -= (s2s1[1] * nz - s2s1[2] * ny) * costheta * 0.5;
                              d_surface(s2,1) -= (s2s1[2] * nx - s2s1[0] * nz) * costheta * 0.5;
                              d_surface(s2,2) -= (s2s1[0] * ny - s2s1[1] * nx) * costheta * 0.5;
                            }
                        }
                    }
                }
            }
          else
            {
              // Case bidim_axi (calculation for a 1radian angle,
              // As it is a ratio surface/volume, the angle as no effect.
              // Volume differential:

              const long s1 = facettes_(facette, 0);
              const long s2 = facettes_(facette, 1);
              const double r1 = sommets_(s1, 0);
              const double y1 = sommets_(s1, 1);
              const double r2 = sommets_(s2, 0);
              const double y2 = sommets_(s2, 1);
              const double L2 = (r2-r1)*(r2-r1) + (y2-y1)*(y2-y1);
              const double L  = sqrt(L2);
              const double nx = normale(facette, 0);
              const double ny = normale(facette, 1);
              const double dv_dn1 = L * (r1 * un_tiers + r2 * un_sixieme);
              const double dv_dn2 = L * (r2 * un_tiers + r1 * un_sixieme);
              d_volume(s1, 0) += nx * dv_dn1;
              d_volume(s1, 1) += ny * dv_dn1;
              d_volume(s2, 0) += nx * dv_dn2;
              d_volume(s2, 1) += ny * dv_dn2;
              // Differentielle de surface
              d_surface(s1, 0) += 0.5 * L + 0.5 * (r1 + r2) * (-ny);
              d_surface(s1, 1) +=           0.5 * (r1 + r2) * (nx);
              d_surface(s2, 0) += 0.5 * L + 0.5 * (r1 + r2) * (ny);
              d_surface(s2, 1) +=           0.5 * (r1 + r2) * (-nx);

              // GB 09/01/2018. Correction to account for the surface differencial
              // from the boundary face wetted by phase 0
              long som[2];
              double r[2];
              som[0] = facettes_(facette, 0);
              som[1] = facettes_(facette, 1);
              r[0] = sommets_(s1, 0);
              r[1] = sommets_(s2, 0);
              for (long i2 = 0; i2 < 2; i2++)
                {
                  const long face = sommet_face_bord_[som[i2]];
                  if (face < 0) // pas une face de bord
                    continue;

                  double costheta = tab_cos_theta(som[i2], 0);
#if TCL_MODEL
                  // Modification of the value of costheta based on
                  // the effect of CL velocity. Stored in tcl.set_theta_app
                  // (but maybe unused there). Value used locally here afterward.
                  // If the TCL model is activated and we're not on a virtual node :
                  if ((sommet_elem_[som[i2]]>0) && flag_tcl)
                    {
                      const double t=temps_physique_;
                      long face_loc;
                      const Domaine_Cl_dis_base& zcl = equation_transport().get_probleme_base().equation(0).domaine_Cl_dis().valeur();
                      const Cond_lim_base& type_cl = zcl.condition_limite_de_la_face_reelle(face,face_loc);
                      const Nom& bc_name = type_cl.frontiere_dis().le_nom();
                      // For each BC, we check its type to see if it's a wall:
                      //   Cerr << "TCL boundary condition is " << type_cl << " for BC " << bc_name;
                      if ( sub_type(Dirichlet_paroi_fixe,type_cl)
                           || sub_type(Dirichlet_paroi_defilante,type_cl) )
                        {
                          Cerr << "  -> for this BC [" <<
                               bc_name <<"], we compute a specific contact angle at TCL." << finl;
                        }
                      else
                        {
                          Cerr << "  -> no modification of contact angle for this BC [" << bc_name <<"] at TCL." << finl;
                          // We're on the symmetry axis, or something else but not at the wall...
                          continue;
                        }
                      FTd_vecteur3 nface = {0., 0., 0.} ;
                      parcours.calculer_normale_face_bord(face, sommets_(som[i2],0), sommets_(som[i2],1), 0., nface[0], nface[1], nface[2]) ;
                      // The other som of the facette is som[1-i2] :
                      // Cerr << " sommet-1 x= " << sommets_(som[1-i2],0) << " y= " << sommets_(som[1-i2],1) << " time_sommet= " << t << finl;
                      // if(dt != 0.) double cl_v = (sommets_(som[i2],0) - x_cl_)/dt;
                      const long isom1 = som[1-i2];
                      // The second vertex of the segment should not be a contact line
                      // so we compute its tangential velocity:

                      const long nb_compo = vit.dimension(1);
                      double norm_vit_som1 = 0.;
                      double vn = 0.;
                      double v_cl = 0.;
                      double v_comp = 0.;
                      // Component of the velocity is calculated along the direction of wall as vw = v - (v.n)n where n is the wall face normal
                      // and v is the velocity vector of the node under consideration.
                      for (long k=0 ; k<nb_compo ; k++)
                        {
                          const double vk = (double) vit(isom1,k);
                          vn += vk*nface[k];
                          Cerr << "Estimated TCL velocity[som=" << isom1
                               << ", compo["<< k<<"]= " << vk << "m/s)" << " face-normal= " << nface[k] <<  finl;
                          Cerr << "Vn= " << vn << finl;
                          norm_vit_som1 += vk*vk;
                        }
                      nface[0] = vn*nface[0];
                      nface[1] = vn*nface[1];
                      nface[2] = vn*nface[2];
                      for (long k=0 ; k<nb_compo ; k++)
                        {
                          v_comp = vit(isom1,k) - nface[k];
                          v_cl += v_comp*v_comp;
                        }
                      v_cl = sqrt(v_cl);
                      Cerr << "v_cl= " << v_cl << " time_v_cl= " << t << finl;

                      // And we assume as a best guess that the contact line
                      // velocity should be approximately that of the first
                      // marker that is not a contact line. We use it
                      // directly to build the Capilary number Ca.
                      const double theta = acos(costheta);
                      Cerr << "theta=" << theta << " time_theta= " << t << finl;
                      const double bubble_center = 0.;
                      const double W = (sommets_(i2, 0) - bubble_center)/(2*2.71*2.71);
                      const double Ca = 2.8e-4*v_cl/5.89e-2;
                      Cerr << "Capillary_number = " << Ca << " time = " << t << finl;
                      double theta_app = pow((theta),3) - 9. * Ca * log(std::max(W, 1.e-20)/l_v);
                      theta_app = pow(std::max(theta_app, 0.),1./3.);
                      Cerr << "theta_after " << theta_app << " time_theta_after= " << t << finl;
                      // We store the apparent contact angle in costheta for
                      // later use in the calculation of the curvature
                      // (it is through this mean that we will consider
                      //  it and try to indirectly satisfy it).
                      costheta = cos(theta_app);
                      // unused?
                      //tcl.set_theta_app(theta_app);
                      Cerr << "[TCL-model] Contact_angle_micro= " << M_PI-theta << " apparent= " << theta_app
                           << " (velocity= " << norm_vit_som1 << " m/s)" << " time= " << t << " theta_app_degree= " << (theta_app/M_PI)*180 << finl;
                    }

                  if (refequation_transport_.non_nul() && (sommet_elem_[som[i2]]>0))
                    {
                      const Transport_Interfaces_FT_Disc& eq_interfaces = refequation_transport_.valeur();
                      const Probleme_base& pb = eq_interfaces.get_probleme_base();
                      Probleme_FT_Disc_gen& pb_ft = ref_cast_non_const(Probleme_FT_Disc_gen, pb);
                      Triple_Line_Model_FT_Disc& tcl = pb_ft.tcl();

                      long face_loc;
                      const Domaine_Cl_dis_base& zcl =
                        equation_transport ().get_probleme_base ().equation (
                          0).domaine_Cl_dis ().valeur ();
                      const Cond_lim_base& type_cl =
                        zcl.condition_limite_de_la_face_reelle (face,
                                                                face_loc);
                      bool is_wall = (sub_type(Dirichlet_paroi_fixe, type_cl) || sub_type(Dirichlet_paroi_defilante,type_cl));

                      if (tcl.is_read_via_file() && is_wall )
                        {
                          double theta_app = tcl.get_theta_app(face)/180.*M_PI;
                          costheta = cos(theta_app);
                          Cerr << "[TCL-model] Contact_angle apparent= " << (theta_app/M_PI)*180 << finl;
                        }
                    }
#endif
                  // Normale unitaire au bord
                  double nfx, nfy, nfz;
                  parcours.calculer_normale_face_bord(face,
                                                      sommets_(som[i2],0), sommets_(som[i2],1), 0.,
                                                      nfx, nfy, nfz);
                  // Ajout de la contribution de la surface:
                  double signe = (i2==0) ? -1. : 1.;
                  d_surface(som[i2], 0) +=  r[i2] * signe * nfy * costheta;
                  d_surface(som[i2], 1) +=  r[i2] * signe * nfx * costheta;
                }
            }
        }
    }

  // On a calcule la contribution de chaque facette reelle aux differents sommets.
  // Certaines contributions ont ete ajoutees a des sommets virtuels, il
  // faut recuperer ces contributions sur le sommet reel.
  desc_sommets().collecter_espace_virtuel(d_volume, MD_Vector_tools::EV_SOMME);
  desc_sommets().collecter_espace_virtuel(d_surface, MD_Vector_tools::EV_SOMME);

  // Calcul de la courbure :
  //   d_surface * d_volume / norme(d_volume)^2
  double dx[3] = { 0., 0., 0. };
#if (defined(PATCH_HYSTERESIS_V2) || defined(PATCH_HYSTERESIS_V3) )
  double norm_L2=0.;
  double norm_Linf=0.;
  long count=0;
  long isom_max_courb=-123456;
#endif

  for (i = 0; i < nsom; i++)
    {
      double c = 0.;
      if (!sommet_virtuel(i))
        {
#if (defined(PATCH_HYSTERESIS_V2) || defined(PATCH_HYSTERESIS_V3) )
          if ((methode_calcul_courbure_contact_line_ == STANDARD)
              || (methode_calcul_courbure_contact_line_ == IMPROVED)
              || (methode_calcul_courbure_contact_line_ == HYSTERESIS)  )
            {
#endif
              // Calcul du vecteur deplacement :
              // Si ce n'est pas une ligne de contact, c'est le vecteur normal,
              // sinon c'est la projection du vecteur normal sur le bord du domaine.
              dx[0] = d_volume(i, 0);
              dx[1] = d_volume(i, 1);
              if (dim == 3)
                dx[2] = d_volume(i, 2);

              const long face_bord = sommet_face_bord_[i];
              if (face_bord >= 0)
                parcours.projeter_vecteur_sur_face(face_bord, dx[0], dx[1], dx[2]);

              double ds_dx = 0.;
              double dv_dx = 0.;
              for (j = 0; j < dim; j++)
                {
                  const double ds = d_surface(i, j);
                  const double dv = d_volume(i, j);
                  ds_dx += ds * dx[j];
                  dv_dx += dv * dx[j];
                }
              if (dv_dx != 0.)
                {
                  c = - ds_dx / dv_dx;
                }
#if (defined(PATCH_HYSTERESIS_V2) || defined(PATCH_HYSTERESIS_V3) )
            }
          if (methode_calcul_courbure_contact_line_ == MIRROR)
            {

              double ds_dv = 0.;
              double dv_dv = 0.;
              for (j = 0; j < dim; j++)
                {
                  const double ds = d_surface(i, j);
                  const double dv = d_volume(i, j);
                  ds_dv += ds * dv;
                  dv_dv += dv * dv;
                }
              if (dv_dv != 0.)
                {
                  c = - ds_dv / dv_dv;
                }

            }
          // Quelle que soit la methode, on print si on a plus d'une passe :
          if ((call>1) && (sommet_face_bord_[i]>=0))
            {
              double delta=std::fabs(courbure_sommets[i] - c);
              norm_L2 += delta*delta;
              count++;
              if (delta>norm_Linf)
                {
                  norm_Linf = delta;
                  isom_max_courb = i;
                }
            }
#endif
        }
#if (defined(PATCH_HYSTERESIS_V2) || defined(PATCH_HYSTERESIS_V3) )
#else
      ArrOfDouble& store_courbure_sommets(courbure_sommets);
      ArrOfBit drapeau_angle_in_range(courbure_sommets.size_array());
      drapeau_angle_in_range = 1;
#endif
      courbure_sommets[i] = c;
      if (methode_calcul_courbure_contact_line_ == NONE)
        {
          if ((call>1) && (sommet_face_bord_[i]>=0))
            //if (call>1)
            courbure_sommets[i] = store_courbure_sommets[i];
        }
      if (methode_calcul_courbure_contact_line_ == WEIGHTED)
        {
          if ((call>1) && (sommet_face_bord_[i]>=0))
            {
              //if (call>1)
              const double c_in = store_courbure_sommets[i];
              const double c_bord = c;
              courbure_sommets[i] = (1-weight_CL_)*c_in+ weight_CL_*c_bord;
            }
        }
      if (methode_calcul_courbure_contact_line_ == HYSTERESIS)
        {
          if ((call>1) && (sommet_face_bord_[i]>=0))
            {
              // Au second passage, la courbure des noeuds de la ligne de contact est :
              const long contact_angle_inside_range = drapeau_angle_in_range.testsetbit(i);
              if (contact_angle_inside_range)
                {
                  // o Soit fournie par une valeur lissee issue de l'interieur du domaine
                  //       si le cos(theta) est dans la gamme autorisee.
                  //       Ainsi, la courbure sera assez reguliere au voisinage de la ligne de contact
                  //       et il n'y aura pas de force generee pour la faire bouger.
                  const double c_in = store_courbure_sommets[i];
                  courbure_sommets[i] = c_in;
                }
              else
                {
                  // o Soit calculee sur le bord avec la valeur limite de cos(theta) la plus proche du domaine authorise
                  //     si l'angle est sur un bord de la gamme.
                  //     Dans ce cas, une inomogeneite de courbure peut conduire a la creation d'un gradient de la
                  //     force de tension de surface et par suite a la mise en mouvement de la ligne de contact.
                  const double c_bord = c;
                  courbure_sommets[i] = c_bord;
                }
            }
        }
    }
#if (defined(PATCH_HYSTERESIS_V2) || defined(PATCH_HYSTERESIS_V3) )
  if ((call>1) && ((methode_calcul_courbure_contact_line_ == MIRROR)
                   || (methode_calcul_courbure_contact_line_ == IMPROVED)
                   || (methode_calcul_courbure_contact_line_ == HYSTERESIS) ))
    {
      norm_Linf = mp_max(norm_Linf);
      count=mp_sum(count);
      if (count>0)
        norm_L2=sqrt(mp_sum(norm_L2))/count;
      Cerr << "Evaluation de la convergence de la courbure. call= " << call << " norm_Linf= " << norm_Linf
           << " norm_L2= " << norm_L2 << " nsom_contact= " << count << finl;
      Journal() << "max_courb on " << isom_max_courb << " value  " << norm_Linf << finl;
    }
#endif
  desc_sommets().echange_espace_virtuel(courbure_sommets);
}

/*! @brief Prepare un tableau de donnees aux sommets ou aux facettes pour conserver les valeurs apres transport.
 *
 * Le transport modifie le descripteur et augmente
 *   le nombre de sommets ou de facettes. En general les tableaux qui contiennent
 *   des valeurs aux sommets ou aux facettes ne sont donc plus valables apres
 *   le transport. On peut les rendre valables comme suit
 *   (exemple avec le tableau courbure aux sommets) :
 *    preparer_tableau_avant_transport(courbure, desc_sommets());
 *    transporter(deplacement);
 *    update_tableau_apres_transport(courbure, desc_sommets());
 *   Attention, le tableau ne devient valable qu'apres avoir appele
 *    update_tableau_apres_transport.
 *
 */
void Maillage_FT_Disc::preparer_tableau_avant_transport(ArrOfDouble& tableau,
                                                        const Desc_Structure_FT& descripteur) const
{
  // Preparation de "collecter_espace_virtuel" : on annule la valeur pour tous
  // les elements virtuels:
  long voisin, i;
  const Descripteur_FT& esp_virt = descripteur.espace_virtuel();
  const ArrOfInt& liste_pe_voisins = esp_virt.pe_voisins();
  const long n_voisins = liste_pe_voisins.size_array();
  for (voisin = 0; voisin < n_voisins; voisin++)
    {
      const long pe_voisin = liste_pe_voisins[voisin];
      const ArrOfInt& elements = esp_virt.elements(pe_voisin);
      const long n = elements.size_array();
      for (i = 0; i < n; i++)
        {
          const long elem = elements[i];
          tableau[elem] = 0.;
        }
    }
}
void Maillage_FT_Disc::preparer_tableau_avant_transport(ArrOfInt& tableau,
                                                        const Desc_Structure_FT& descripteur) const
{
  // Preparation de "collecter_espace_virtuel" : on annule la valeur pour tous
  // les elements virtuels:
  long voisin, i;
  const Descripteur_FT& esp_virt = descripteur.espace_virtuel();
  const ArrOfInt& liste_pe_voisins = esp_virt.pe_voisins();
  const long n_voisins = liste_pe_voisins.size_array();
  for (voisin = 0; voisin < n_voisins; voisin++)
    {
      const long pe_voisin = liste_pe_voisins[voisin];
      const ArrOfInt& elements = esp_virt.elements(pe_voisin);
      const long n = elements.size_array();
      for (i = 0; i < n; i++)
        {
          const long elem = elements[i];
          tableau[elem] = 0;
        }
    }
}

/*! @brief Voir preparer_tableau_avant_transport
 *
 */
void Maillage_FT_Disc::preparer_tableau_avant_transport(DoubleTab& tableau,
                                                        const Desc_Structure_FT& descripteur) const
{
  const long dim = tableau.dimension(1);
  long voisin, i, j;
  const Descripteur_FT& esp_virt = descripteur.espace_virtuel();
  const ArrOfInt& liste_pe_voisins = esp_virt.pe_voisins();
  const long n_voisins = liste_pe_voisins.size_array();
  for (voisin = 0; voisin < n_voisins; voisin++)
    {
      const long pe_voisin = liste_pe_voisins[voisin];
      const ArrOfInt& elements = esp_virt.elements(pe_voisin);
      const long n = elements.size_array();
      for (i = 0; i < n; i++)
        {
          const long elem = elements[i];
          for (j = 0; j < dim; j++)
            tableau(elem,j) = 0.;
        }
    }
}

/*! @brief Voir preparer_tableau_avant_transport
 *
 */
void Maillage_FT_Disc::update_tableau_apres_transport(ArrOfDouble& tableau,
                                                      long new_size,
                                                      const Desc_Structure_FT& descripteur) const
{
  // Le transport des interfaces cree de nouveaux sommets/facettes en fin de tableau.
  // On doit donc augmenter la taille du tableau et determiner la valeur
  // de ces nouveaux elements.
  const long old_size = tableau.size_array();
  tableau.resize_array(new_size);
  long i;
  for (i = old_size; i < new_size; i++)
    tableau[i] = 0.;
  // Somme sur tous les processeurs de la valeur stockee pour chaque element.
  // La valeur est non nulle uniquement sur le processeur qui possedait
  // l'element avant le transport (voir preparer_tableau_avant_transport).
  // Le nouveau proprietaire recupere ainsi la valeur avant transport.
  descripteur.collecter_espace_virtuel(tableau, MD_Vector_tools::EV_SOMME);
  descripteur.echange_espace_virtuel(tableau);
}

void Maillage_FT_Disc::update_tableau_apres_transport(ArrOfInt& tableau,
                                                      long new_size,
                                                      const Desc_Structure_FT& descripteur) const
{
  // Le transport des interfaces cree de nouveaux sommets/facettes en fin de tableau.
  // On doit donc augmenter la taille du tableau et determiner la valeur
  // de ces nouveaux elements.
  const long old_size = tableau.size_array();
  tableau.resize_array(new_size);
  long i;
  for (i = old_size; i < new_size; i++)
    tableau[i] = 0;
  // Somme sur tous les processeurs de la valeur stockee pour chaque element.
  // La valeur est non nulle uniquement sur le processeur qui possedait
  // l'element avant le transport (voir preparer_tableau_avant_transport).
  // Le nouveau proprietaire recupere ainsi la valeur avant transport.
  descripteur.collecter_espace_virtuel(tableau, MD_Vector_tools::EV_SOMME);
  descripteur.echange_espace_virtuel(tableau);
}

/*! @brief Voir preparer_tableau_avant_transport
 *
 */
void Maillage_FT_Disc::update_tableau_apres_transport(DoubleTab& tableau,
                                                      long new_size,
                                                      const Desc_Structure_FT& descripteur) const
{
  const long old_size = tableau.dimension(0);
  const long dim = tableau.dimension(1);
  tableau.resize(new_size, dim);
  long i, j;
  for (i = old_size; i < new_size; i++)
    for (j = 0; j < dim; j++)
      tableau(i,j) = 0.;
  descripteur.collecter_espace_virtuel(tableau, MD_Vector_tools::EV_SOMME);
  descripteur.echange_espace_virtuel(tableau);
}

const Intersections_Elem_Facettes& Maillage_FT_Disc::intersections_elem_facettes() const
{
  assert(statut_ >= PARCOURU);
  return intersections_elem_facettes_;
}

Transport_Interfaces_FT_Disc& Maillage_FT_Disc::equation_transport()
{
  return refequation_transport_.valeur();
}

const Transport_Interfaces_FT_Disc& Maillage_FT_Disc::equation_transport() const
{
  return refequation_transport_.valeur();
}

const Equation_base& Maillage_FT_Disc::equation_associee() const
{
  return equation_transport();
}

long Maillage_FT_Disc::type_statut() const
{
  if (statut_==RESET)
    return 0;
  else if (statut_==MINIMAL)
    return 1;
  else if (statut_==PARCOURU)
    return 2;
  else if (statut_==COMPLET)
    return 3;
  else
    {
      Cerr<<"Le type de statut n est pas fixe"<<finl;
      exit();
    }
  return -1;
}

/*! @brief creation d'un tableau aux sommets du maillage Meme principe que Domaine::creer_tableau_sommets()
 *
 */
void Maillage_FT_Disc::creer_tableau_sommets(Array_base& x, Array_base::Resize_Options opt) const
{
  const MD_Vector& md = desc_sommets().get_md_vector();
  MD_Vector_tools::creer_tableau_distribue(md, x, opt);
}

/*! @brief creation d'un tableau aux sommets du maillage Meme principe que Domaine::creer_tableau_elements()
 *
 */
void Maillage_FT_Disc::creer_tableau_elements(Array_base& x, Array_base::Resize_Options opt) const
{
  const MD_Vector& md = desc_facettes().get_md_vector();
  MD_Vector_tools::creer_tableau_distribue(md, x, opt);
}

Schema_Comm_FT Maillage_FT_Disc::get_schema_comm_FT() const
{
  return schema_comm_domaine_;
}
void Maillage_FT_Disc::set_is_solid_particle(long is_solid_particle)
{
  is_solid_particle_=is_solid_particle;
}