/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Champ_Post_Operateur_Eqn.h>
#include <Probleme_base.h>
#include <Postraitement.h>
#include <Domaine_VF.h>
#include <Nom.h>
#include <Champ_Generique_refChamp.h>
#include <Discretisation_base.h>
#include <Equation_base.h>
#include <Operateur_base.h>
#include <Param.h>
#include <MD_Vector_composite.h>
#include <Transport_K_Eps_base.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Mod_turb_hyd_RANS.h>

Implemente_instanciable_sans_constructeur(Champ_Post_Operateur_Eqn,"Operateur_Eqn|Champ_Post_Operateur_Eqn",Champ_Generique_Operateur_base);

Sortie& Champ_Post_Operateur_Eqn::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}
Champ_Post_Operateur_Eqn::Champ_Post_Operateur_Eqn()
{
  numero_op_=-1;
  numero_source_=-1;
  sans_solveur_masse_=0;
  compo_=-1;
}

void Champ_Post_Operateur_Eqn::set_param(Param& param)
{
  Champ_Generique_Operateur_base::set_param(param);
  param.ajouter("numero_source",&numero_source_);
  param.ajouter("numero_op",&numero_op_);
  param.ajouter_flag("sans_solveur_masse",&sans_solveur_masse_);
  param.ajouter("compo",&compo_);
}

Entree& Champ_Post_Operateur_Eqn::readOn(Entree& s )
{
  Champ_Generique_Operateur_base::readOn(s);

  return s ;
}

void Champ_Post_Operateur_Eqn::verification_cas_compo() const
{
  // On applique compo a un vecteur
  const Nature_du_champ& nature_ch=ref_eq_.valeur().inconnue().valeur().nature_du_champ();
  if ((nature_ch != vectoriel) && (compo_ != -1 ))
    {
      Cerr<<"Error in Champ_Post_Operateur_Eqn::verification_cas_compo()"<<finl;
      Cerr<<"It isn't possible to get a component from a non vectoriel field " <<finl;
      exit();
    }

  // Verification de compo
  const int nb_compo= ref_eq_.valeur().inconnue().valeur().nb_comp();
  if ((compo_<-1)||(compo_>nb_compo-1))
    {
      Cerr<<"Error in Champ_Post_Operateur_Eqn::verification_cas_compo()"<<finl;
      Cerr<<"compo="<<compo_<<" is not allowed."<<" You must give a value between 0 and "<< nb_compo-1<<finl;
      exit();
    }
  // Verifier qu'on n'est pas en VDF
  const Domaine_dis_base& domaine_dis = ref_eq_.valeur().domaine_dis().valeur();
  if ((domaine_dis.que_suis_je().debute_par("Domaine_VDF")) && (compo_ != -1 ))
    {
      Cerr<<"Error in Champ_Post_Operateur_Eqn::verification_cas_compo()"<<finl;
      Cerr<<"The option compo is not available in case of VDF discretization"<<finl;
      exit();
    }
}

void Champ_Post_Operateur_Eqn::completer(const Postraitement_base& post)
{

  Champ_Gen_de_Champs_Gen::completer(post);

  const Probleme_base& Pb = ref_cast(Postraitement,post).probleme();
  int numero_eq_=-1;
  bool iskeps = false;
  if (sub_type(Champ_Generique_refChamp,get_source(0)))
    {

      Champ espace_stockage;
      const Champ_base& mon_champ = get_source(0).get_champ(espace_stockage);
      if (sub_type(Champ_Inc_base,mon_champ))
        {
          const Champ_Inc_base& mon_champ_inc = ref_cast(Champ_Inc_base,mon_champ);

          {

            //On recupere l equation alors qu elle n est pas encore associee au Champ_Inc
            //On parcours les equatiosn du probleme et on identifie celle qui correspond au champ inconnu

            int nb_eq = Pb.nombre_d_equations();
            int i=0;

            while (i<nb_eq)
              {
                const Equation_base& eq_test = Pb.equation(i);
                if ((eq_test.inconnue().le_nom() == mon_champ_inc.le_nom()))
                  {
                    numero_eq_=i;
                    break;
                  }
                else if (mon_champ_inc.le_nom() == "K_Eps")
                  {
                    const RefObjU& modele_turbulence = eq_test.get_modele(TURBULENCE);
                    if (sub_type(Modele_turbulence_hyd_K_Eps, modele_turbulence.valeur()) || sub_type(Modele_turbulence_hyd_K_Eps_Realisable, modele_turbulence.valeur()) )
                      {
                        const Mod_turb_hyd_RANS& le_mod_RANS = ref_cast(Mod_turb_hyd_RANS, eq_test.get_modele(TURBULENCE).valeur());
                        const Transport_K_Eps_base& transportkeps = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
                        if ((transportkeps.inconnue().le_nom() == mon_champ_inc.le_nom()))
                          {
                            numero_eq_=i;
                            iskeps = true;
                            break;
                          }
                      }
                  }
                i++;
              }
          }
        }
    }
  if (numero_eq_==-1)
    {
      Cerr<<"Champ_Post_Operateur_Eqn can be applied only to equation unknown."<<finl;
      exit();
    }

  if (!iskeps)
    ref_eq_=Pb.equation(numero_eq_);
  else
    {
      const Mod_turb_hyd_RANS& le_mod_RANS = ref_cast(Mod_turb_hyd_RANS, Pb.equation(numero_eq_).get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      ref_eq_= ref_cast(Equation_base,eqn);
    }
  int ok=0;
  const Equation_base& eqn=ref_eq_.valeur();
  const MD_Vector& mdf = eqn.inconnue().valeurs().get_md_vector(),
                   md = sub_type(MD_Vector_composite, mdf.valeur()) ? ref_cast(MD_Vector_composite, mdf.valeur()).get_desc_part(0) : mdf;
  const Domaine_VF& zvf= ref_cast( Domaine_VF,ref_eq_.valeur().domaine_dis().valeur());
  if (md== zvf.face_sommets().get_md_vector())
    {
      localisation_inco_=FACE;
      ok=1;
    }
  if (md== zvf.domaine().les_elems().get_md_vector())
    {
      localisation_inco_=ELEMENT;
      ok=1;
    }
  if (md == zvf.domaine().les_sommets().get_md_vector())
    {
      ok=1;
      localisation_inco_=NODE;
    }
  if (ok==0)
    {
      Cerr<<"Error in "<<que_suis_je()<<" unknown localisation"<<finl;
      exit();
    }
  verification_cas_compo();
}

const Champ_base& Champ_Post_Operateur_Eqn::get_champ_without_evaluation(Champ& espace_stockage) const
{
  espace_stockage = ref_eq_.valeur().inconnue();
  return espace_stockage;
}

const Champ_base& Champ_Post_Operateur_Eqn::get_champ_compo_without_evaluation(Champ& espace_stockage) const
{

  Champ_Fonc espace_stockage_fonc;
  //Champ source_espace_stockage;

  double temps=0.;
  Nom directive;

  switch (localisation_inco_)
    {
    case ELEMENT:
      directive="CHAMP_ELEM";
      break;
    case NODE:
      directive="CHAMP_SOMMETS";
      break;
    case FACE:
      directive="CHAMP_FACE";
      break;
    default:
      Cerr<<"error in Champ_Post_Operateur_Eqn::get_champ"<<finl;
      exit();

    }
  int nb_comp = 1;
  ref_eq_.valeur().discretisation().discretiser_champ(directive,ref_eq_->domaine_dis().valeur(),"oooo","unit", nb_comp,temps,espace_stockage_fonc);
  espace_stockage = espace_stockage_fonc;
  espace_stockage.valeur().fixer_nature_du_champ(scalaire);

  return espace_stockage;
}

const Champ_base& Champ_Post_Operateur_Eqn::get_champ(Champ& espace_stockage) const
{
  // On commence par construire le champ vectoriel complet
  Champ espace_stockage_complet;
  espace_stockage_complet = get_champ_without_evaluation(espace_stockage_complet);
  DoubleTab& es = (espace_stockage_complet.valeurs());

  //if (ref_eq_->schema_temps().temps_courant()!=0)
  {
    if (numero_op_!=-1)
      {
        // certains calculer  sont faux !!!! il faudrait tous les recoder en res =0 ajouter();
        es=0;
        Operateur().ajouter(ref_eq_->operateur(numero_op_).mon_inconnue().valeurs(),es);
      }
    else
      ref_eq_->sources()(numero_source_).calculer(es);
    if (!sans_solveur_masse_)
      ref_eq_->solv_masse().valeur().appliquer_impl(es); //On divise par le volume
    // Hack: car Masse_PolyMAC_Face::appliquer_impl ne divise pas par le volume (matrice de masse)....
    if (ref_eq_->solv_masse().valeur().que_suis_je()=="Masse_PolyMAC_Face")
      {
        //Cerr << "Volumic source terms on faces with PolyMAC can't be post-processed yet." << finl;
        Cerr << "Warning, source terms on faces with PolyMAC are post-processed as S*dV not as volumic source terms S." << finl;
        Cerr << "Cause Masse_PolyMAC_Face::appliquer_impl do not divide per volume." << finl;
        //Process::exit();
      }
  }
  es.echange_espace_virtuel();
  // apres je recupere juste la composante si elle est demandee
  if (compo_>-1)
    {
      // on prepare l'espace de stockage pour une composante
      Champ espace_stockage_compo;
      espace_stockage_compo = get_champ_compo_without_evaluation(espace_stockage_compo);
      DoubleTab& es_compo = (espace_stockage_compo.valeurs());
      int nb_pos = es.dimension(0);
      for (int i=0; i<nb_pos; i++)
        {
          es_compo(i) = es(i,compo_);
        }
      espace_stockage = espace_stockage_compo;
    }
  else
    espace_stockage = espace_stockage_complet;

  return espace_stockage.valeur();
}

const Noms Champ_Post_Operateur_Eqn::get_property(const Motcle& query) const
{
//Creation des composantes serait a faire de maniere dynamique (Eqn_...)

  Motcles motcles(2);
  motcles[0] = "composantes";
  motcles[1] = "unites";

  int rang = motcles.search(query);

  switch(rang)
    {
    case 0:
      {
        if (compo_==-1)
          {
            int nb_comp= ref_eq_.valeur().inconnue().valeur().nb_comp();
            Noms compo(nb_comp);
            for (int i=0; i<nb_comp; i++)
              {
                Nom nume(i);
                compo[i] = nom_post_+nume;
              }
            return compo;
          }
        else
          {
            Noms compo(1);
            compo[0] = nom_post_;
            return compo;
          }
      }

    case 1:
      {
        if (compo_==-1)
          {
            int nb_comp= ref_eq_.valeur().inconnue().valeur().nb_comp();
            Noms unites(nb_comp);
            //Noms source_unites = get_source(0).get_property("unites");
            for (int i=0; i<nb_comp; i++)
              {
                unites[i] = "unit";
              }
            return unites;
          }
        else
          {
            // J'utilise un vecteur car la methode renvoie Noms
            Noms unites(1);
            unites[0] = "unit";
            return unites;
          }
      }

    }
  return Champ_Gen_de_Champs_Gen::get_property(query);
}

Entity Champ_Post_Operateur_Eqn::get_localisation(const int index) const
{

  return localisation_inco_;

}
//Nomme le champ en tant que source par defaut
//"Eqn_" + nom_champ_source
void Champ_Post_Operateur_Eqn::nommer_source()
{
  if (nom_post_=="??")
    {
      Nom nom_post_source;
      //const Noms nom = get_source(0).get_property("nom");
      nom_post_source = "Eqn_s";
      nom_post_source += Nom(numero_source_);
      nom_post_source += "_o" ;
      nom_post_source += Nom(numero_op_);
      if (compo_!=-1)
        {
          Nom nume(compo_);
          nom_post_source += nume;
        }
      nommer(nom_post_source);
    }

}
const Operateur_base& Champ_Post_Operateur_Eqn::Operateur() const
{
  return ref_eq_->operateur(numero_op_).l_op_base();
}

Operateur_base& Champ_Post_Operateur_Eqn::Operateur()
{
  return ref_eq_->operateur(numero_op_).l_op_base();
}
