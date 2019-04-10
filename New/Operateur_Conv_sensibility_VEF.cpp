/****************************************************************************
* Copyright (c) 2019, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Operateur_Conv_sensibility_VEF.cpp
// Directory : $BALTIK_COUPLAGE_ROOT/src/New
//
/////////////////////////////////////////////////////////////////////////////

#include <Operateur_Conv_sensibility_VEF.h>
#include <Op_Conv_VEF_base.h>
#include <Navier_Stokes_std_sensibility.h>
#include <DoubleTrav.h>
#include <Op_Conv_VEF_Face.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Debog.h>


Implemente_instanciable( Operateur_Conv_sensibility_VEF, "Op_Conv_sensibility_VEF_P1NC",Operateur_Conv_sensibility ) ;


Sortie&  Operateur_Conv_sensibility_VEF::printOn(Sortie& os) const
{
  return os;
}


Entree&  Operateur_Conv_sensibility_VEF::readOn(Entree& is)
{
  Cerr << " Operateur_Conv_sensibility_VEF ::readOn " << finl;
  op_conv.associer_eqn(equation());
  op_conv.associer_vitesse(la_vitesse.valeur());
  op_conv.lire(is);
  Cerr << "Operateur_Conv_sensibility_VEF : " << op_conv.que_suis_je() << finl;
  return is;
}

void  Operateur_Conv_sensibility_VEF::associer (const Zone_dis& zone_dis ,
                                                const Zone_Cl_dis& zone_cl_dis,
                                                const Champ_Inc& inco )
{
  Cerr << " Operateur_Conv_sensibility_VEF::associer" << finl;
  const Zone_VEF& zvef = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zclvef = ref_cast(Zone_Cl_VEF,zone_cl_dis.valeur());

  la_zone_vef = zvef;
  la_zcl_vef = zclvef;
  la_zone_vef.valeur().creer_tableau_faces(fluent);
  Operateur_Conv_sensibility::associer(zone_dis,zone_cl_dis,inco);
}

DoubleTab& Operateur_Conv_sensibility_VEF::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  Cerr << "Operateur_Conv_sensibility_VEF::ajouter" << finl;

  //Check convection scheme discretization type:
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());
  int convectionSchemeDiscrType; // amont=0, muscl=1, centre=2
  opConvVEFFace.get_type_op(convectionSchemeDiscrType);

  if(convectionSchemeDiscrType==0) // Convection scheme discr "amont".
    {
      const Navier_Stokes_std_sensibility& eq = ref_cast(Navier_Stokes_std_sensibility, equation());
      const DoubleTab& state = eq.get_state_field();
      ajouter_Lstate_sensibility_Amont(state, inco, resu);
      ajouter_Lsensibility_state_Amont(inco, state, resu);
      opConvVEFbase.modifier_flux(*this); // Multiplication by density in case of incompressible Navier Stokes
    }
  else
    {
      Cout << "Sensibility cannot use currently convection scheme " << op_conv.type() <<". Try available Sensibility amont convection scheme." << finl;
      exit();
    }
  resu.echange_espace_virtuel();
  Debog::verifier("resu dansOperateur_Conv_sensibility_VEF::ajouter", resu);
  return resu;
}

void Operateur_Conv_sensibility_VEF::ajouter_Lstate_sensibility_Amont(const DoubleTab& state_field, const DoubleTab& inco, DoubleTab& resu ) const
{
  Cout<<"ajouter_Lstate_sensibility "<<finl;
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, la_zone_vef.valeur());
  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& facenormales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const int nfa7 = zone_VEF.type_elem().nb_facette();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();
  int nfac = zone.nb_faces_elem();
  int nsom = zone.nb_som_elem();
  int nb_bord = zone_VEF.nb_front_Cl();
  const IntTab& les_elems=zone.les_elems();
  int option_appliquer_cl_dirichlet = 0 ;

  // Definition d'un tableau pour un traitement special des schemas pres des bords
  if (traitement_pres_bord_.size_array()!=nb_elem_tot)
    {
      traitement_pres_bord_.resize_array(nb_elem_tot);
      traitement_pres_bord_=0;
      ArrOfInt est_un_sommet_de_bord_(zone_VEF.nb_som_tot());
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              int size = zone_VEF.face_sommets().dimension(1);
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                for (int som=0; som<size; som++)
                  {
                    int face = le_bord.num_face(ind_face);
                    est_un_sommet_de_bord_(zone_VEF.face_sommets(face,som))=1;
                  }
            }
        }
      for (int elem=0; elem<nb_elem_tot; elem++)
        {
          if (rang_elem_non_std(elem)!=-1)
            traitement_pres_bord_(elem)=1;
          else
            {
              for (int n_som=0; n_som<nsom; n_som++)
                if (est_un_sommet_de_bord_(les_elems(elem,n_som)))
                  traitement_pres_bord_(elem)=1;
            }
        }
      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(zone_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_(num_face) = 1;
                }
            }
        }
    }

  // Pour le traitement de la convection on distingue les polyedres
  // standard qui ne "voient" pas les conditions aux limites et les
  // polyedres non standard qui ont au moins une face sur le bord.
  // Un polyedre standard a n facettes sur lesquelles on applique le
  // schema de convection.
  // Pour un polyedre non standard qui porte des conditions aux limites
  // de Dirichlet, une partie des facettes sont portees par les faces.
  // En bref pour un polyedre le traitement de la convection depend
  // du type (triangle, tetraedre ...) et du nombre de faces de Dirichlet.

  const Elem_VEF_base& type_elemvef= zone_VEF.type_elem().valeur();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  double psc, psc_inco;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20;
  int ncomp_ch_transporte=(inco.nb_dim() == 1?1:inco.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          nb_faces_perio+=le_bord.nb_faces();
        }
    }

  DoubleTab tab;
  if (ncomp_ch_transporte == 1)
    tab.resize(nb_faces_perio);
  else
    tab.resize(nb_faces_perio,ncomp_ch_transporte);

  nb_faces_perio=0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (ncomp_ch_transporte == 1)
                tab(nb_faces_perio) = resu(num_face);
              else
                for (int comp=0; comp<ncomp_ch_transporte; comp++)
                  tab(nb_faces_perio,comp) = resu(num_face,comp);
              nb_faces_perio++;
            }
        }
    }

  int comp0;
  int nb_faces_ = zone_VEF.nb_faces();
  ArrOfInt face(nfac);
  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble cc(dimension);
  DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab xsom(nsom,dimension);


  ArrOfDouble vs_inco(dimension);
  ArrOfDouble vc_inco(dimension);
  DoubleTab vsom_inco(nsom,dimension);

  // Dimensionnement du tableau des flux convectifs au bord du domaine de calcul
  DoubleTab& flux_b = ref_cast(DoubleTab,flux_bords_);
  int nb_faces_bord=zone_VEF.nb_faces_bord();
  flux_b.resize(nb_faces_bord,ncomp_ch_transporte);
  flux_b = 0.;

  const IntTab& KEL=type_elemvef.KEL();



  // Les polyedres non standard sont ranges en 2 groupes dans la Zone_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans la zone
  // boucle sur les polys

  for (poly=0; poly<nb_elem_tot; poly++)
    {
      int contrib = 0;
      // calcul des numeros des faces du polyedre
      for (face_adj=0; face_adj<nfac; face_adj++)
        {
          int face_ = elem_faces(poly,face_adj);
          face(face_adj)= face_;
          if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
        }
      //
      if (contrib)
        {
          for (j=0; j<dimension; j++)
            {
              vs(j) = state_field(face(0),j);
              vs_inco(j) = inco(face(0),j);
              for (i=1; i<nfac; i++)
                {
                  vs(j)+= state_field(face(i),j);
                  vs_inco(j)+= inco(face(i),j);
                }
            }
          // calcul de la vitesse aux sommets des polyedres
          if (istetra==1)
            {
              for (i=0; i<nsom; i++)
                for (j=0; j<dimension; j++)
                  {
                    vsom(i,j) = (vs(j) - dimension*state_field(face(i),j));
                    vsom_inco(i,j) = (vs_inco(j) - dimension*inco(face(i),j));
                  }
            }
          else
            {
              Cout << "Sensibility is currently working only with Tetra_VEF (3D) or Tri_VEF (2D)." << finl;
              exit();
            }

          // Determination du type de CL selon le rang
          rang = rang_elem_non_std(poly);
          int itypcl = (rang==-1 ? 0 : zone_Cl_VEF.type_elem_Cl(rang));

          // calcul de vc (a l'intersection des 3 facettes) vc vs vsom
          calcul_vc(face,vc,vs,vsom,state_field,itypcl);
          calcul_vc(face,vc_inco,vs_inco,vsom_inco,inco,itypcl);


          // Boucle sur les facettes du polyedre non standard:
          for (fa7=0; fa7<nfa7; fa7++)
            {
              num10 = face(KEL(0,fa7));
              num20 = face(KEL(1,fa7));

              // normales aux facettes
              if (rang==-1)
                for (i=0; i<dimension; i++)
                  cc[i] = facette_normales(poly, fa7, i);
              else
                for (i=0; i<dimension; i++)
                  cc[i] = normales_facettes_Cl(rang,fa7,i);

              // Calcul des vitesses en C,S,S2 les 3 extremites de la fa7 et M le centre de la fa7
              double psc_c=0,psc_s=0,psc_m,psc_s2=0;
              double psc_c_inco=0,psc_s_inco=0,psc_m_inco;
              if (dimension==2)
                {
                  for (i=0; i<2; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                      psc_c_inco+=vc_inco[i]*cc[i];
                      psc_s_inco+=vsom_inco(KEL(2,fa7),i)*cc[i];
                    }
                  psc_m=(psc_c+psc_s)/2.;
                  psc_m_inco=(psc_c_inco+psc_s_inco)/2.;
                }
              else
                {
                  for (i=0; i<3; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                      psc_s2+=vsom(KEL(3,fa7),i)*cc[i];
                    }
                  psc_m=(psc_c+psc_s+psc_s2)/3.;
                }
              // On applique les CL de Dirichlet si num1 ou num2 est une face avec CL de Dirichlet
              // auquel cas la fa7 coincide avec la face num1 ou num2 -> C est au centre de la face
              if (option_appliquer_cl_dirichlet)
                if (est_une_face_de_dirichlet_(num10) || est_une_face_de_dirichlet_(num20))
                  {
                    psc_m = psc_c;
                    psc_m_inco = psc_c_inco;
                  }
              // Determination de la face amont pour M
              int face_amont_m;
              if (option_appliquer_cl_dirichlet)  Cout<<"psc_m_inco =  "<<psc_m_inco<<finl;
              if (psc_m >= 0)
                {
                  face_amont_m = num10;
                }
              else
                {
                  face_amont_m = num20;
                }

              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                {
                  double flux;
                  double inco_m = (ncomp_ch_transporte==1?inco(face_amont_m):inco(face_amont_m,comp0));
                  flux = inco_m*psc_m;

                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num10) -= flux;
                      resu(num20) += flux;
                      if (num10<nb_faces_bord) flux_b(num10,0) += flux;
                      if (num20<nb_faces_bord) flux_b(num20,0) -= flux;
                    }
                  else
                    {
                      resu(num10,comp0) -= flux;
                      resu(num20,comp0) += flux;
                      if (num10<nb_faces_bord) flux_b(num10,comp0) += flux;
                      if (num20<nb_faces_bord) flux_b(num20,comp0) -= flux;
                    }

                }// boucle sur comp
            } // fin de la boucle sur les facettes
        }
    } // fin de la boucle
  resu.echange_espace_virtuel();
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0.;
              psc_inco = 0.;
              for (i=0; i<dimension; i++)
                {
                  psc += state_field(num_face,i)*facenormales(num_face,i);
                  psc_inco += inco(num_face,i)*facenormales(num_face,i);
                }
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= psc*inco(num_face);
                    flux_b(num_face,0) -= psc*inco(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*inco(num_face,i);
                      flux_b(num_face,i) -= psc*inco(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= psc_inco*la_sortie_libre.val_ext(num_face-num1);
                      flux_b(num_face,0) -= psc_inco*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) -= psc_inco*la_sortie_libre.val_ext(num_face-num1,i);
                        flux_b(num_face,i) -= psc_inco*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          ArrOfInt fait(le_bord.nb_faces());
          fait = 0;
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (fait[num_face-num1] == 0)
                {
                  voisine = la_cl_perio.face_associee(num_face-num1) + num1;

                  if (ncomp_ch_transporte == 1)
                    {
                      diff1 = resu(num_face)-tab(nb_faces_perio);
                      diff2 = resu(voisine)-tab(nb_faces_perio+voisine-num_face);
                      resu(voisine)  += diff1;
                      resu(num_face) += diff2;
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        resu(voisine,comp)  += diff1;
                        resu(num_face,comp) += diff2;
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }
}
void Operateur_Conv_sensibility_VEF::ajouter_Lsensibility_state_Amont(const DoubleTab& inco, const DoubleTab& state_field, DoubleTab& resu ) const
{
  Cout<<"ajouter_Lstate_sensibility "<<finl;
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, la_zone_vef.valeur());
  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& facenormales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const int nfa7 = zone_VEF.type_elem().nb_facette();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();
  int nfac = zone.nb_faces_elem();
  int nsom = zone.nb_som_elem();
  int nb_bord = zone_VEF.nb_front_Cl();
  const IntTab& les_elems=zone.les_elems();
  int option_appliquer_cl_dirichlet = 0 ;


  // Definition d'un tableau pour un traitement special des schemas pres des bords
  if (traitement_pres_bord_.size_array()!=nb_elem_tot)
    {
      traitement_pres_bord_.resize_array(nb_elem_tot);
      traitement_pres_bord_=0;
      ArrOfInt est_un_sommet_de_bord_(zone_VEF.nb_som_tot());
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              int size = zone_VEF.face_sommets().dimension(1);
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                for (int som=0; som<size; som++)
                  {
                    int face = le_bord.num_face(ind_face);
                    est_un_sommet_de_bord_(zone_VEF.face_sommets(face,som))=1;
                  }
            }
        }
      for (int elem=0; elem<nb_elem_tot; elem++)
        {
          if (rang_elem_non_std(elem)!=-1)
            traitement_pres_bord_(elem)=1;
          else
            {
              for (int n_som=0; n_som<nsom; n_som++)
                if (est_un_sommet_de_bord_(les_elems(elem,n_som)))
                  traitement_pres_bord_(elem)=1;
            }
        }

      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(zone_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_(num_face) = 1;
                }
            }
        }
    }

  // Pour le traitement de la convection on distingue les polyedres
  // standard qui ne "voient" pas les conditions aux limites et les
  // polyedres non standard qui ont au moins une face sur le bord.
  // Un polyedre standard a n facettes sur lesquelles on applique le
  // schema de convection.
  // Pour un polyedre non standard qui porte des conditions aux limites
  // de Dirichlet, une partie des facettes sont portees par les faces.
  // En bref pour un polyedre le traitement de la convection depend
  // du type (triangle, tetraedre ...) et du nombre de faces de Dirichlet.

  const Elem_VEF_base& type_elemvef= zone_VEF.type_elem().valeur();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20;
  int ncomp_ch_transporte=(state_field.nb_dim() == 1?1:state_field.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          nb_faces_perio+=le_bord.nb_faces();
        }
    }

  DoubleTab tab;
  if (ncomp_ch_transporte == 1)
    tab.resize(nb_faces_perio);
  else
    tab.resize(nb_faces_perio,ncomp_ch_transporte);

  nb_faces_perio=0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (ncomp_ch_transporte == 1)
                tab(nb_faces_perio) = resu(num_face);
              else
                for (int comp=0; comp<ncomp_ch_transporte; comp++)
                  tab(nb_faces_perio,comp) = resu(num_face,comp);
              nb_faces_perio++;
            }
        }
    }

  int comp0;
  int nb_faces_ = zone_VEF.nb_faces();
  ArrOfInt face(nfac);
  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble cc(dimension);
  DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab xsom(nsom,dimension);

  ArrOfDouble vs_state(dimension);
  ArrOfDouble vc_state(dimension);
  DoubleTab vsom_state(nsom,dimension);

  // Dimensionnement du tableau des flux convectifs au bord du domaine de calcul
  DoubleTab& flux_b = ref_cast(DoubleTab,flux_bords_);
  int nb_faces_bord=zone_VEF.nb_faces_bord();
  flux_b.resize(nb_faces_bord,ncomp_ch_transporte);
  flux_b = 0.;

  const IntTab& KEL=type_elemvef.KEL();

  // Les polyedres non standard sont ranges en 2 groupes dans la Zone_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans la zone
  // boucle sur les polys

  for (poly=0; poly<nb_elem_tot; poly++)
    {
      int contrib = 0;
      // calcul des numeros des faces du polyedre
      for (face_adj=0; face_adj<nfac; face_adj++)
        {
          int face_ = elem_faces(poly,face_adj);
          face(face_adj)= face_;
          if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
        }
      //
      if (contrib)
        {
          for (j=0; j<dimension; j++)
            {
              vs(j) = inco(face(0),j);
              vs_state(j) = state_field(face(0),j);
              for (i=1; i<nfac; i++)
                {
                  vs(j)+= inco(face(i),j);
                  vs_state(j) += state_field(face(i),j);
                }
            }
          // calcul de la vitesse aux sommets des polyedres
          if (istetra==1)
            {
              for (i=0; i<nsom; i++)
                for (j=0; j<dimension; j++)
                  {
                    vsom(i,j) = (vs(j) - dimension*inco(face(i),j));
                    vsom_state(i,j) = (vs_state(j) - dimension*state_field(face(i),j));
                  }
            }
          else
            {
              Cout << "Sensibility is currently working only with Tetra_VEF (3D) or Tri_VEF (2D)." << finl;
              exit();
            }

          // Determination du type de CL selon le rang
          rang = rang_elem_non_std(poly);
          int itypcl = (rang==-1 ? 0 : zone_Cl_VEF.type_elem_Cl(rang));

          // calcul de vc (a l'intersection des 3 facettes) vc vs vsom
          calcul_vc(face,vc,vs,vsom,inco,itypcl);
          calcul_vc(face,vc_state,vs_state,vsom_state,state_field,itypcl);

          // Boucle sur les facettes du polyedre non standard:
          for (fa7=0; fa7<nfa7; fa7++)
            {
              num10 = face(KEL(0,fa7));
              num20 = face(KEL(1,fa7));

              // normales aux facettes
              if (rang==-1)
                for (i=0; i<dimension; i++)
                  cc[i] = facette_normales(poly, fa7, i);
              else
                for (i=0; i<dimension; i++)
                  cc[i] = normales_facettes_Cl(rang,fa7,i);

              // Calcul des vitesses en C,S,S2 les 3 extremites de la fa7 et M le centre de la fa7
              double psc_c=0,psc_s=0,psc_m,psc_s2=0;
              double psc_c_state=0,psc_s_state=0,psc_m_state;
              if (dimension==2)
                {
                  for (i=0; i<2; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                      psc_c_state+=vc_state[i]*cc[i];
                      psc_s_state+=vsom_state(KEL(2,fa7),i)*cc[i];
                    }
                  psc_m=(psc_c+psc_s)/2.;
                  psc_m_state=(psc_c_state+psc_s_state)/2.;
                }
              else
                {
                  for (i=0; i<3; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                      psc_s2+=vsom(KEL(3,fa7),i)*cc[i];
                    }
                  psc_m=(psc_c+psc_s+psc_s2)/3.;
                }
              // On applique les CL de Dirichlet si num1 ou num2 est une face avec CL de Dirichlet
              // auquel cas la fa7 coincide avec la face num1 ou num2 -> C est au centre de la face
              //int appliquer_cl_dirichlet=0;
              if (option_appliquer_cl_dirichlet)
                if (est_une_face_de_dirichlet_(num10) || est_une_face_de_dirichlet_(num20))
                  {
                    //appliquer_cl_dirichlet = 1;
                    psc_m = psc_c;
                    psc_m_state = psc_c_state;
                  }

              // Determination de la face amont pour M
              int face_amont_m;
              if (option_appliquer_cl_dirichlet)  Cout<<"psc_m_state = "<<psc_m_state<<finl;
              if (psc_m >= 0)
                {
                  face_amont_m = num10;
                }
              else
                {
                  face_amont_m = num20;
                }

              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                {
                  double flux;
                  double inco_m = (ncomp_ch_transporte==1?state_field(face_amont_m):state_field(face_amont_m,comp0));

                  flux = inco_m*psc_m;

                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num10) -= flux;
                      resu(num20) += flux;
                      if (num10<nb_faces_bord) flux_b(num10,0) += flux;
                      if (num20<nb_faces_bord) flux_b(num20,0) -= flux;
                    }
                  else
                    {
                      resu(num10,comp0) -= flux;
                      resu(num20,comp0) += flux;
                      if (num10<nb_faces_bord) flux_b(num10,0) += flux;
                      if (num20<nb_faces_bord) flux_b(num20,0) -= flux;
                    }

                }// boucle sur comp
            } // fin de la boucle sur les facettes
        }
    } // fin de la boucle
  resu.echange_espace_virtuel();
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += inco(num_face,i)*facenormales(num_face,i);
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= psc*state_field(num_face);
                    flux_b(num_face,0) -= psc*state_field(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*state_field(num_face,i);
                      flux_b(num_face,i) -= psc*state_field(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= psc*la_sortie_libre.val_ext(num_face-num1);
                      flux_b(num_face) -= psc*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                        flux_b(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          ArrOfInt fait(le_bord.nb_faces());
          fait = 0;
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (fait[num_face-num1] == 0)
                {
                  voisine = la_cl_perio.face_associee(num_face-num1) + num1;

                  if (ncomp_ch_transporte == 1)
                    {
                      diff1 = resu(num_face)-tab(nb_faces_perio);
                      diff2 = resu(voisine)-tab(nb_faces_perio+voisine-num_face);
                      resu(voisine)  += diff1;
                      resu(num_face) += diff2;
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        resu(voisine,comp)  += diff1;
                        resu(num_face,comp) += diff2;
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }
}


void Operateur_Conv_sensibility_VEF::calcul_vc(const ArrOfInt& Face, ArrOfDouble& vc, const ArrOfDouble& vs, const DoubleTab& vsom,
                                               const DoubleTab& vitesse_face,int type_cl) const
{
// Operateur_Conv_sensibility_VEF::calcul_vc(...) is based on Tetra_VEF::calcul_vc(...) and on Tri_VEF::calcul_vc(...).
  if (dimension == 2) //Tri_VEF
    {

      switch(type_cl)
        {
        case 0: // le triangle n'a pas de Face de Dirichlet
          {
            vc[0] = vs[0]/3;
            vc[1] = vs[1]/3;
            break;
          }

        case 1: // le triangle a une Face de Dirichlet :la Face 2
          {
            vc[0]= vitesse_face(Face[2],0);
            vc[1]= vitesse_face(Face[2],1);
            break;
          }

        case 2: // le triangle a une Face de Dirichlet :la Face 1
          {
            vc[0]= vitesse_face(Face[1],0);
            vc[1]= vitesse_face(Face[1],1);
            break;
          }

        case 4: // le triangle a une Face de Dirichlet :la Face 0
          {
            vc[0]= vitesse_face(Face[0],0);
            vc[1]= vitesse_face(Face[0],1);
            break;
          }

        case 3: // le triangle a deux faces de Dirichlet :les faces 1 et 2
          {
            vc[0]= vsom(0,0);
            vc[1]= vsom(0,1);
            break;
          }

        case 5: // le triangle a deux faces de Dirichlet :les faces 0 et 2
          {
            vc[0]= vsom(1,0);
            vc[1]= vsom(1,1);
            break;
          }

        case 6: // le triangle a deux faces de Dirichlet :les faces 0 et 1
          {
            vc[0]= vsom(2,0);
            vc[1]= vsom(2,1);
            break;
          }

        } // switch end

    } // 2D end

  else if (dimension == 3) //Tetra_VEF
    {

      int comp;
      switch(type_cl)
        {

        case 0: // le tetraedre n'a pas de Face de Dirichlet
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.25*vs[comp];
            break;
          }

        case 1: // le tetraedre a une Face de Dirichlet : KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face(Face[3],comp);
            break;
          }

        case 2: // le tetraedre a une Face de Dirichlet : KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face(Face[2],comp);
            break;
          }

        case 4: // le tetraedre a une Face de Dirichlet : KEL1
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face(Face[1],comp);
            break;
          }

        case 8: // le tetraedre a une Face de Dirichlet : KEL0
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face(Face[0],comp);
            break;
          }

        case 3: // le tetraedre a deux faces de Dirichlet : KEL3 et KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(0,comp) + vsom(1,comp));
            break;
          }

        case 5: // le tetraedre a deux faces de Dirichlet : KEL3 et KEL1
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(0,comp) + vsom(2,comp));
            break;
          }

        case 6: // le tetraedre a deux faces de Dirichlet : KEL1 et KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(0,comp) + vsom(3,comp));
            break;
          }

        case 9: // le tetraedre a deux faces de Dirichlet : KEL0 et KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(1,comp) + vsom(2,comp));
            break;
          }

        case 10: // le tetraedre a deux faces de Dirichlet : KEL0 et KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(1,comp) + vsom(3,comp));
            break;
          }

        case 12: // le tetraedre a deux faces de Dirichlet : KEL0 et KEL1
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5*(vsom(2,comp) + vsom(3,comp));
            break;
          }

        case 7: // le tetraedre a trois faces de Dirichlet : KEL1, KEL2 et KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vsom(0,comp);
            break;
          }

        case 11: // le tetraedre a trois faces de Dirichlet : KEL0,KEL2 et KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vsom(1,comp);
            break;
          }

        case 13: // le tetraedre a trois faces de Dirichlet : KEL0, KEL1 et KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vsom(2,comp);
            break;
          }

        case 14: // le tetraedre a trois faces de Dirichlet : KEL0, KEL1 et KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vsom(3,comp);
            break;
          }

        } // switch end

    } // 3D end
}
void Operateur_Conv_sensibility_VEF::remplir_fluent(DoubleVect& tab_fluent) const
{

  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, la_zone_vef.valeur());
  const Champ_Inc_base& velocity=vitesse();
  const DoubleTab& vitesse_face=velocity.valeurs();
  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& face_normales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const int nfa7 = zone_VEF.type_elem().nb_facette();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();
  int nfac = zone.nb_faces_elem();
  int nsom = zone.nb_som_elem();
  const IntTab& sommet_elem = zone.les_elems();
  DoubleVect& fluent_ = ref_cast(DoubleVect, tab_fluent);


  const Elem_VEF_base& type_elemvef= zone_VEF.type_elem().valeur();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;


  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang ,itypcl;
  int num1,num2,num_som;

  ArrOfInt face(nfac);
  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  DoubleTab vsom(nsom,dimension);
  ArrOfDouble cc(dimension);

  // Dimensionnement du tableau des flux convectifs au bord du domaine de calcul
  const IntTab& KEL=type_elemvef.KEL();

  // On remet a zero le tableau qui sert pour
  // le calcul du pas de temps de stabilite
  fluent_ = 0;

  // Les polyedres non standard sont ranges en 2 groupes dans la Zone_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans la zone

  // boucle sur les polys
  for (poly=0; poly<nb_elem_tot; poly++)
    {
      rang = rang_elem_non_std(poly);
      // On cherche, pour un elem qui n'est pas de bord (rang==-1),
      // si un des sommets est sur un bord (tableau des sommets) (C MALOD 17/07/2007)

      if (rang==-1)
        itypcl=0;
      else
        itypcl=zone_Cl_VEF.type_elem_Cl(rang);

      // calcul des numeros des faces du polyedre
      for (face_adj=0; face_adj<nfac; face_adj++)
        face(face_adj)= elem_faces(poly,face_adj);

      for (j=0; j<dimension; j++)
        {
          vs(j) = vitesse_face(face(0),j);
          for (i=1; i<nfac; i++)
            vs(j)+= vitesse_face(face(i),j);
        }
      // calcul de la vitesse aux sommets des polyedres
      if (istetra==1)
        {
          for (i=0; i<nsom; i++)
            for (j=0; j<dimension; j++)
              vsom(i,j) = (vs(j) - dimension*vitesse_face(face(i),j));
        }
      else
        {
          // pour que cela soit valide avec les hexa (c'est + lent a calculer...)
          int ncomp;
          for (j=0; j<nsom; j++)
            {
              num_som = sommet_elem(poly,j);
              for (ncomp=0; ncomp<dimension; ncomp++)
                vsom(j,ncomp) =velocity.valeur_a_sommet_compo(num_som,poly,ncomp);
            }
        }

      // calcul de vc (a l'intersection des 3 facettes) vc vs vsom proportionnelles a la prosite
      calcul_vc(face,vc,vs,vsom,vitesse(),itypcl);
      // Boucle sur les facettes du polyedre non standard:
      for (fa7=0; fa7<nfa7; fa7++)
        {
          num1 = face(KEL(0,fa7));
          num2 = face(KEL(1,fa7));
          // normales aux facettes
          if (rang==-1)
            {
              for (i=0; i<dimension; i++)
                cc[i] = facette_normales(poly, fa7, i);
            }
          else
            {
              for (i=0; i<dimension; i++)
                cc[i] = normales_facettes_Cl(rang,fa7,i);
            }
          // On applique le schema de convection a chaque sommet de la facette

          double psc_c=0,psc_s=0,psc_m,psc_s2=0;
          if (dimension==2)
            {
              for (i=0; i<dimension; i++)
                {
                  psc_c+=vc[i]*cc[i];
                  psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                }
              psc_m=(psc_c+psc_s)/2.;
            }
          else
            {
              for (i=0; i<dimension; i++)
                {
                  psc_c+=vc[i]*cc[i];
                  psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                  psc_s2+=vsom(KEL(3,fa7),i)*cc[i];
                }
              psc_m=(psc_c+psc_s+psc_s2)/3.;
            }

          // int amont,dir;
          if (psc_m >= 0)
            {
              // amont = num1;
              fluent_(num2)  += psc_m;
              //dir=0;
            }
          else
            {
              //amont = num2;
              fluent_(num1)  -= psc_m;
              //dir=1;
            }

        } // fin de la boucle sur les facettes
    } // fin de la boucle

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  int nb_bord = zone_VEF.nb_front_Cl();
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1b = le_bord.num_premiere_face();
          int num2b = num1b + le_bord.nb_faces();
          for (num_face=num1b; num_face<num2b; num_face++)
            {
              psc = 0;
              for (i=0; i<dimension; i++)
                psc += vitesse_face(num_face,i)*face_normales(num_face,i);
              if (psc>0)
                ;
              else
                fluent_(num_face) -= psc;
            }
        }
    }
}
double Operateur_Conv_sensibility_VEF::calculer_dt_stab() const
{
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  const DoubleVect& volumes_entrelaces_Cl = zone_Cl_VEF.volumes_entrelaces_Cl();
  remplir_fluent(fluent);
  double dt_face,dt_stab =1.e30;

  // On traite les conditions aux limites
  // Si une face porte une condition de Dirichlet on n'en tient pas compte
  // dans le calcul de dt_stab
  for (int n_bord=0; n_bord<zone_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if ( (sub_type(Dirichlet,la_cl.valeur()))
           ||
           (sub_type(Dirichlet_homogene,la_cl.valeur()))
         )
        ;
      else
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              dt_face = volumes_entrelaces_Cl(num_face)/(fluent[num_face]+1.e-30);
              dt_stab = (dt_face < dt_stab) ? dt_face : dt_stab;
            }
        }
    }

  // On traite les faces internes non standard
  int ndeb = zone_VEF.premiere_face_int();
  int nfin = zone_VEF.premiere_face_std();

  for (int num_face=ndeb; num_face<nfin; num_face++)
    {
      dt_face = volumes_entrelaces_Cl(num_face)/(fluent[num_face]+DMINFLOAT);
      dt_stab =(dt_face < dt_stab) ? dt_face : dt_stab;
    }

  // On traite les faces internes standard
  ndeb = nfin;
  nfin = zone_VEF.nb_faces();
  for (int num_face=ndeb; num_face<nfin; num_face++)
    {
      dt_face = volumes_entrelaces(num_face)/(fluent[num_face]+DMINFLOAT);
      dt_stab =(dt_face < dt_stab) ? dt_face : dt_stab;
    }

  dt_stab = Process::mp_min(dt_stab);
  // astuce pour contourner le type const de la methode
  Operateur_Conv_sensibility_VEF& op = ref_cast_non_const(Operateur_Conv_sensibility_VEF,*this);
  op.fixer_dt_stab_conv(dt_stab);


  return dt_stab;
}
