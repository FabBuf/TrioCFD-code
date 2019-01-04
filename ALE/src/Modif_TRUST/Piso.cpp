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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Piso.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Schemas_Temps
// Version:     /main/29
//
//////////////////////////////////////////////////////////////////////////////

#include <Piso.h>
#include <Zone_VF.h>
#include <Navier_Stokes_std.h>
#include <EChaine.h>
#include <Debog.h>
#include <Matrice_Bloc.h>
#include <Assembleur_base.h>
#include <Statistiques.h>
#include <Schema_Temps_base.h>
#include <DoubleTrav.h>
#include <Fluide_Quasi_Compressible.h>
#include <Dirichlet.h>
#include <Probleme_base.h>

#include <Op_Conv_ALE_VEF.h>
#include <Domaine_ALE.h>

Implemente_instanciable_sans_constructeur(Piso,"Piso",Simpler);

Implemente_instanciable_sans_constructeur(Implicite,"Implicite",Piso);

Piso::Piso()
{
  nb_corrections_max_ = 21;
  avancement_crank_ = 0;
}

Implicite::Implicite()
{
  avancement_crank_ = 1;
}

Sortie& Piso::printOn(Sortie& os ) const
{
  return Simpler::printOn(os);
}

Entree& Piso::readOn(Entree& is )
{
  Simpler::readOn(is);
  return is;
}

Sortie& Implicite::printOn(Sortie& os ) const
{
  return Piso::printOn(os);
}

Entree& Implicite::readOn(Entree& is )
{
  Piso::readOn(is);
  return is;
}


Entree& Piso::lire(const Motcle& motlu,Entree& is)
{
  Motcles les_mots(1);
  {
    les_mots[0] = "nb_corrections_max";
  }

  int rang = les_mots.search(motlu);
  switch(rang)
    {
    case 0:
      {
        is >> nb_corrections_max_;
        if (nb_corrections_max_ < 2)
          {
            Cerr<<"There must be at least two corrections steps for the PISO algorithm."<<finl;
            exit();
          }
        break;
      }

    default :
      {
        Cerr << "Keyword : " << motlu << " is not understood in " << que_suis_je() << finl;
        exit();
      }
    }
  return is;
}



void test_imposer_cond_lim(Equation_base& eqn,DoubleTab& current2,const char * msg,int flag)
{
  return;
  DoubleTab& present = eqn.inconnue().futur();
  DoubleTab sauv(present);
  const Schema_Temps_base& sch = eqn.probleme().schema_temps();
  eqn.zone_Cl_dis()->imposer_cond_lim(eqn.inconnue(),sch.temps_courant()+sch.pas_de_temps());
  present -= sauv;
  // BM, je remplace max_abs par mp_pax_abs: du coup la methode doit etre appelee simultanement par tous les procs.
  double ecart_max=mp_max_abs_vect(present);
  Cout<<msg <<" "<<ecart_max<<finl;

  if ((ecart_max>1e-10))
    abort();
  present = sauv;
}

//Entree Un ; Pn
//Sortie Un+1 = U***_k ; Pn+1 = P**_k
//n designe une etape temporelle

void Piso::iterer_NS(Equation_base& eqn,DoubleTab& current,DoubleTab& pression,
                     double dt,Matrice_Morse& matrice,double seuil_resol,DoubleTrav& secmem,int nb_ite,int& converge)
{
  Parametre_implicite& param_eqn = get_and_set_parametre_implicite(eqn);
  SolveurSys& le_solveur_ = param_eqn.solveur();
  converge = 1;
  if (nb_ite>1) return;

  Navier_Stokes_std& eqnNS = ref_cast(Navier_Stokes_std,eqn);
  DoubleTrav gradP(current);
  DoubleTrav correction_en_pression(pression);
  DoubleTrav resu(current);
  int is_QC = eqn.probleme().is_QC();

  double vitesse_norme,vitesse_norme_old ;
  double pression_norme,pression_norme_old ;
  vitesse_norme_old = 1.e10 ;
  pression_norme_old = 1.e10 ;

  // int deux_entrees = 0;
  //if (current.nb_dim()==2) deux_entrees = 1;
  Operateur_Grad& gradient = eqnNS.operateur_gradient();
  Operateur_Div& divergence = eqnNS.operateur_divergence();

  //Renewing ALE Jacobians
  if (eqn.probleme().domaine().que_suis_je()=="Domaine_ALE")
    {
      int TimeStepNr=eqn.probleme().schema_temps().nb_pas_dt();
      Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, eqn.probleme().domaine());

      DoubleTab New_ALEjacobian_Old=dom_ale.getNewJacobian(); //New  value for ALEjacobian_old
      DoubleTab New_ALEjacobian_New(New_ALEjacobian_Old);
      Op_Conv_ALE_VEF& opALEforJacob=ref_cast(Op_Conv_ALE_VEF, eqnNS.get_terme_convectif().valeur());
      opALEforJacob.calculateALEjacobian(New_ALEjacobian_New); //New  value for ALEjacobian_new
      dom_ale.update_ALEjacobians(New_ALEjacobian_Old, New_ALEjacobian_New, TimeStepNr); // Update new values of ALEjacobian_old and ALEjacobian_new saved in Domaine_ALE
    }
  //End of renewing ALE Jacobians

  //int nb_comp = 1;
  //int nb_dim = current.nb_dim();

  //  if(nb_dim==2)    nb_comp = current.dimension(1);

  //Construction de matrice et resu
  //matrice = A[Un] = M/delta_t + CONV +DIFF
  //resu =  A[Un]Un -(A[Un]Un-Ss) + Sv -BtPn
  gradient.calculer(pression,gradP);
  resu -= gradP;

  //Adding ALE convection term start
  if (eqn.probleme().domaine().que_suis_je()=="Domaine_ALE")
    {
      Cerr << "Adding ALE contribution..." << finl;
      Op_Conv_ALE& opale=ref_cast(Op_Conv_ALE, eqnNS.get_terme_convectif().valeur());
      DoubleTrav ALE(resu); // copie de la structure, initialise a zero
      opale.ajouterALE(current, ALE);
      ALE.echange_espace_virtuel();
      //solveur_masse.appliquer(ALE); do not need to divide by mass
      //ALE.echange_espace_virtuel();
      resu+=ALE;                        //resu + ALE convection tem
      resu.echange_espace_virtuel();
      Debog::verifier("Piso::iterer_NS resu after adding ALE",resu);
    }
  //Adding ALE convection term end

  eqnNS.assembler_avec_inertie(matrice,current,resu);
  le_solveur_.valeur().reinit();

  int nb_faces = resu.size()/dimension;

  // Adding Jacobians to "matrice" start. Jn+1[A[Un]]
  if (eqn.probleme().domaine().que_suis_je()=="Domaine_ALE")
    {
      Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, eqn.probleme().domaine());
      DoubleTab ALEjacobian_New=dom_ale.getNewJacobian();

      int MatriceNbLines=matrice.nb_lignes();
      int dim=0;
      int num_faceJacobian=0;
      for (int num_face=0; num_face<MatriceNbLines; num_face++)
        {
          if(num_face % nb_faces == 0 && num_face != 0)
            {
              dim++;
              num_faceJacobian-=nb_faces;
            }
          matrice(num_face,num_face)*=ALEjacobian_New(num_faceJacobian,dim);
          num_faceJacobian++;
        }
    }
  // Adding Jacobians to "matrice" end


  //Adding ALE Jacobians to "resu" start to obtain Jn+1[-gradP + ALE convective term + Sv + Ss]+Jn[(M/dt)Un].
  if (eqn.probleme().domaine().que_suis_je()=="Domaine_ALE")
    {
      Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, eqn.probleme().domaine());
      DoubleTab ALEjacobian_New=dom_ale.getNewJacobian();
      DoubleTab ALEjacobian_Old=dom_ale.getOldJacobian();

      // Obtaining (M/dt)Un start.
      DoubleTrav deltaJMdtUn(resu); // copie de la structure, initialise a zero. To hold (M/dt)Un.
      double timestep=eqn.probleme().schema_temps().pas_de_temps();
      int pen=0;
      eqn.solv_masse().ajouter_masse(timestep,deltaJMdtUn,eqn.inconnue().passe(),pen);
      // Obtaining (M/dt)Un end.

      for (int num_face=0; num_face<nb_faces; num_face++) // Obtaining Jn+1[-gradP + ALE convective term + (M/dt)Un + Sv + Ss]+(Jn-Jn+1)[(M/dt)Un]=Jn+1[-gradP + ALE convective term + Sv + Ss]+Jn[(M/dt)Un].
        {
          for (int dim=0; dim<dimension; dim++)
            {
              resu(num_face,dim)=ALEjacobian_New(num_face,dim)*resu(num_face,dim)+(ALEjacobian_Old(num_face,dim)-ALEjacobian_New(num_face,dim))*deltaJMdtUn(num_face,dim);
            }
        }
      resu.echange_espace_virtuel();
    }
  // Adding Jacobians to "resu" end

  //Construction de matrice_en_pression_2 = BD-1Bt[Un]
  //Assemblage reeffectue seulement pour algorithme Piso (avancement_crank_==0)
  Matrice& matrice_en_pression_2 = eqnNS.matrice_pression();
  SolveurSys& solveur_pression_ = eqnNS.solveur_pression();
  if (avancement_crank_==0)
    {
      assembler_matrice_pression_implicite(eqnNS,matrice,matrice_en_pression_2);
      solveur_pression_.valeur().reinit();
    }

  //Etape predicteur
  //Resolution du systeme A[Un]U* = -BtPn + Sv + Ss
  //current = U*
  le_solveur_.resoudre_systeme(matrice,resu,current);

  test_imposer_cond_lim(eqn,current,"apres resolution ",0);
  current.echange_espace_virtuel();
  Debog::verifier("Piso::iterer_NS current apres solveur",current);

  //Calcul de secmem = BU* (en incompressible) BU* -drho/dt (en quasi-compressible)
  if (is_QC)
    {
      if (with_d_rho_dt_)
        {
          Fluide_Quasi_Compressible& fluide_QC = ref_cast(Fluide_Quasi_Compressible,eqn.milieu());
          fluide_QC.secmembre_divU_Z(secmem);
          secmem *= -1;
        }
      else secmem = 0;
      divergence.ajouter(current,secmem);
    }
  else
    divergence.calculer(current,secmem);
  secmem *= -1;
  secmem.echange_espace_virtuel();
  Debog::verifier("Piso::iterer_NS secmem",secmem);
  // GF il ne faut pas modifier le scd membre le terme en du/dt au bord a deja ete pris en compte dans la resolution precedente
  //  eqnNS.assembleur_pression().valeur().modifier_secmem(secmem);

  //Etape de correction 1
  Cout << "Solving mass equation :" << finl;
  //Description du cas implicite
  //Resolution du systeme (BD-1Bt)P' = Bu* (D-1 = M-1 pour le cas implicite)
  //correction_en_pression = P' pour Piso et correction_en_pression = delta_t*P' pour implicite
  solveur_pression_.resoudre_systeme(matrice_en_pression_2.valeur(),
                                     secmem,correction_en_pression);
  correction_en_pression.echange_espace_virtuel();
  Debog::verifier("Piso::iterer_NS arpes correction_pression",correction_en_pression);

  if (avancement_crank_==1)
    {
      //Calcul de Bt(delta_t*delta_P)
      gradient.valeur().multvect(correction_en_pression,gradP);
      eqn.solv_masse().appliquer(gradP);

      //Calcul de Un+1 = U* -delta_t*delta_P
      current -= gradP;
      current.echange_espace_virtuel();
      divergence.calculer(current,secmem);

      //Calcul de Pn+1 = Pn + (delta_t*delta_P)/delat_t
      correction_en_pression /= dt;
      pression += correction_en_pression;
      eqnNS.assembleur_pression().valeur().modifier_solution(pression);
      pression.echange_espace_virtuel();
      if (is_QC)
        {
          // on redivise par rho_np_1 avant de sortir
          diviser_par_rho_np1_face(eqn,current);
        }
      return;
    }

  // calcul de la correction en vitesse premiere etape (DU' =-Bt P)

  //Calcul de P* = Pn + P'
  pression += correction_en_pression;
  eqnNS.assembleur_pression().valeur().modifier_solution(pression);

  //Resolution du systeme D[Un]U' = -BtP'
  DoubleTrav correction_en_vitesse(current);
  calculer_correction_en_vitesse(correction_en_pression,gradP,correction_en_vitesse,matrice,gradient);

  //Calcul de U** = U* + U'
  current += correction_en_vitesse;
  test_imposer_cond_lim(eqn,current,"apres premiere correction ",0);
  Debog::verifier("Piso::iterer_NS arpes cor pression",pression);
  Debog::verifier("Piso::iterer_NS arpes cor vitesse",current);


  //Etape correcteur 2
  for (int compt=0; compt<nb_corrections_max_-1; compt++)
    {
      correction_en_vitesse.echange_espace_virtuel();

      //Resolution du systeme D resu = EU' + (resu2=0) pour stocker resu = D-1EU'
      DoubleTrav resu2(resu);
      int status = inverser_par_diagonale(matrice,resu2,correction_en_vitesse,resu);

      if (status!=0) exit();
      // calcul de P''  BD-1Bt P''= -div(D-1EU')

      resu.echange_espace_virtuel();
      //Calcul de B(D-1EU')
      divergence.calculer(resu,secmem);
      secmem *= -1;
      secmem.echange_espace_virtuel();

      //Resolution du systeme (BD-1Bt)P'' = (BD-1E)U'
      //correction_en_pression = P''
      correction_en_pression = 0;
      solveur_pression_.resoudre_systeme(matrice_en_pression_2.valeur(),
                                         secmem,correction_en_pression);

#ifdef DEBUG
      // Pour verifier que le calcul du gradient ne modifie pas la pression
      DoubleTrav correction_en_pression_mod(pression);
      correction_en_pression_mod = correction_en_pression;
#endif
      //Resolution du systeme D[Un]U'' = -BtP''
      //correction_en_vitesse = U''
      calculer_correction_en_vitesse(correction_en_pression,gradP,correction_en_vitesse,matrice,gradient);

#ifdef DEBUG
      correction_en_pression_mod -= correction_en_pression;
      assert(max_abs(correction_en_pression_mod)==0);
#endif

      //Calcul de U'' = U'' + D-1EU'
      correction_en_vitesse += resu;
      // ajout des increments

      vitesse_norme = mp_norme_vect(correction_en_vitesse);
      pression_norme = mp_norme_vect(correction_en_pression);

      if ( (vitesse_norme>vitesse_norme_old) || (pression_norme>pression_norme_old) )
        {
          Cout <<"PISO : "<< compt+1 <<" corrections to perform the projection."<< finl;
          if (is_QC)
            {
              // on redivise par rho_np_1 avant de sortir
              diviser_par_rho_np1_face(eqn,current);
            }
          return ;
        }

      vitesse_norme_old = vitesse_norme;
      pression_norme_old = pression_norme;

      //Calcul de P** = P* + P''
      pression += correction_en_pression;
      eqnNS.assembleur_pression().valeur().modifier_solution(pression);

      //Calcul de U*** = U** + U''
      current += correction_en_vitesse;
      test_imposer_cond_lim(eqn,current,"apres correction (int)__LINE__",0);

      Debog::verifier("Piso::iterer_NS apres correction pression",pression);
      Debog::verifier("Piso::iterer_NS apres correction vitesse",current);
    }
  if (is_QC)
    {
      diviser_par_rho_np1_face(eqn,current);
      //ref_cast(Navier_Stokes_QC,eqn).impr_impl(eqnNS, Cout);
    }
  current.echange_espace_virtuel();
  // divergence.calculer(current, secmem); Cerr<<" ici DIVU  "<<mp_max_abs_vect(secmem)<<finl;;
  Cout <<"PISO : "<<nb_corrections_max_<<" corrections to perform the projection."<< finl;
}
