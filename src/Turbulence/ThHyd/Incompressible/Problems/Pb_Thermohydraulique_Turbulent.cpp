/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        Pb_Thermohydraulique_Turbulent.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Problems
//
//////////////////////////////////////////////////////////////////////////////

#include <Pb_Thermohydraulique_Turbulent.h>
#include <Fluide_Incompressible.h>
#include <Les_mod_turb.h>
#include <Verif_Cl.h>
#include <Verif_Cl_Turb.h>
#include <Mod_turb_hyd_RANS.h>

Implemente_instanciable(Pb_Thermohydraulique_Turbulent,"Pb_Thermohydraulique_Turbulent",Pb_Fluide_base);


/*! @brief Simple appel a: Pb_Fluide_base::printOn(Sortie&) Ecrit le probleme sur un flot de sortie.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Pb_Thermohydraulique_Turbulent::printOn(Sortie& os) const
{
  return Pb_Fluide_base::printOn(os);
}


/*! @brief Simple appel a: Pb_Fluide_base::readOn(Entree&) Lit le probleme a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Pb_Thermohydraulique_Turbulent::readOn(Entree& is)
{
  return Pb_Fluide_base::readOn(is);
}

/*! @brief Renvoie le nombre d'equation, Renvoie 2 car il y a 2 equations a un probleme de
 *
 *     thermo-hydraulique turbulent:
 *      - l'equation de Navier Stokes
 *      - l'equation de la thermique de type Convection_Diffusion_Temperature_Turbulent
 *
 * @return (int) le nombre d'equation
 */
int Pb_Thermohydraulique_Turbulent::nombre_d_equations() const
{
  return 2;
}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_Turbulent si i=0 Renvoie l'equation de la thermique de type
 *
 *     Convection_Diffusion_Temperature_Turbulent si i=1
 *     (version const)
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation correspondante a l'index
 */
const Equation_base& Pb_Thermohydraulique_Turbulent::equation(int i) const
{
  if ( !( i==0 || i==1 ) )
    {
      Cerr << "\nError in Pb_Thermohydraulique_Turbulent::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_thermique;

}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_Turbulent si i=0 Renvoie l'equation de la thermique de type
 *
 *     Convection_Diffusion_Temperature_Turbulent si i=1
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation correspondante a l'index
 */
Equation_base& Pb_Thermohydraulique_Turbulent::equation(int i)
{
  if ( !( i==0 || i==1 ) )
    {
      Cerr << "\nError in Pb_Thermohydraulique_Turbulent::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_thermique;
}



/*! @brief Associe le milieu au probleme Le milieu doit etre de type fluide incompressible
 *
 * @param (Milieu_base& mil) le milieu physique a associer au probleme
 * @throws mauvais type de milieu physique
 */
void Pb_Thermohydraulique_Turbulent::associer_milieu_base(const Milieu_base& mil)
{
  if sub_type(Fluide_Incompressible,mil)
    {
      const Fluide_base& mi = ref_cast(Fluide_base,mil);
      eq_hydraulique.associer_milieu_base(mi);
      eq_thermique.associer_milieu_base(mi);
    }
  else
    {
      Cerr << "Un milieu de type " << mil.que_suis_je() << " ne peut etre associe a " << finl;
      Cerr << "un probleme de type Pb_Thermohydraulique_Turbulent " << finl;
      exit();
    }
}


/*! @brief Teste la compatibilite des equations de la thermique et de l'hydraulique.
 *
 * Les tests se font sur les conditions
 *     aux limites discretisees de chaque equation et sur les
 *     modeles de turbulences respectifs des equations
 *     de l'hydraulique et de la thermique (qui doivent etre de la meme famille).
 *     Appel la fonction de librairie hors classe:
 *       tester_compatibilite_hydr_thermique(const Domaine_Cl_dis&,const Domaine_Cl_dis&)
 *
 * @return (int) renvoie toujours 1
 * @throws modeles de turbulence de famille differente pour
 * l'hydraulique et la thermique
 */
int Pb_Thermohydraulique_Turbulent::verifier()
{
  const Domaine_Cl_dis& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis& domaine_Cl_th = eq_thermique.domaine_Cl_dis();

  // Verification de la compatibilite des conditions aux limites:
  tester_compatibilite_hydr_thermique(domaine_Cl_hydr,domaine_Cl_th);
  if ( sub_type(Mod_turb_hyd_RANS, eq_hydraulique.get_modele(TURBULENCE).valeur() ))
    {
      const Mod_turb_hyd_RANS& le_mod_RANS = ref_cast(Mod_turb_hyd_RANS, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Domaine_Cl_dis& domaine_Cl_turb = eqn.domaine_Cl_dis();
      tester_compatibilite_hydr_turb(domaine_Cl_hydr, domaine_Cl_turb);
    }
  /*
    // Verification de la compatibilite des modeles de turbulence:
    const Mod_turb_hyd& le_mod_turb_hyd = eq_hydraulique.modele_turbulence();
    const Modele_turbulence_scal_base& le_mod_turb_th = ref_cast(Modele_turbulence_scal_base,eq_thermique.get_modele(TURBULENCE).valeur());

    if (!sub_type(Modele_turbulence_hyd_K_Eps_Bas_Reynolds,le_mod_turb_hyd.valeur()))
      {
        if ((!sub_type(Modele_turbulence_scal_Prandtl,le_mod_turb_th))
            && (le_mod_turb_th.que_suis_je()!="Modele_turbulence_scal_sous_maille_dyn_VDF")
            && (le_mod_turb_th.que_suis_je()!="Modele_turbulence_scal_Dispersion_Thermique_InterAss"))
          {
            Cerr << "Les modeles de turbulence ne sont pas de la meme famille" << finl;
            Cerr << "pour l'hydraulique et la thermique" << finl;
            exit();
          }
      }
    else
      {
        if  ( (!sub_type(Modele_turbulence_scal_Fluctuation_Temperature,le_mod_turb_th)) &&
              (!sub_type(Modele_turbulence_scal_Prandtl,le_mod_turb_th)) &&
              (!sub_type(Modele_turbulence_scal_Fluctuation_Temperature_W,le_mod_turb_th)))
          {
            Cerr << "Les modeles de turbulence ne sont pas de la meme famille" << finl;
            Cerr << "pour l'hydraulique et la thermique" << finl;
            exit();
          }
      }
  */
  return 1;
}

