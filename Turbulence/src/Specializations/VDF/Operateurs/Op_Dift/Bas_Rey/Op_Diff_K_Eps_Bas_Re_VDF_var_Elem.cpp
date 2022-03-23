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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Op_Diff_K_Eps_Bas_Re_VDF_var_Elem.cpp
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Eps_Bas_Re_VDF_var_Elem.h>
#include <Champ_P0_VDF.h>
#include <Modele_turbulence_hyd_K_Eps_2_Couches.h>
#include <Modele_turbulence_hyd_K_Eps_Bas_Reynolds.h>

Implemente_instanciable_sans_constructeur(Op_Diff_K_Eps_Bas_Re_VDF_var_Elem,"Op_Diff_K_Eps_Bas_Re_VDF_var_P0_VDF",Op_Diff_K_Eps_Bas_Re_VDF_base);
implemente_It_VDF_Elem(Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem)

Sortie& Op_Diff_K_Eps_Bas_Re_VDF_var_Elem::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Op_Diff_K_Eps_Bas_Re_VDF_var_Elem::readOn(Entree& s ) { return s ; }

void Op_Diff_K_Eps_Bas_Re_VDF_var_Elem::associer_diffusivite_turbulente()
{
  assert(mon_equation.non_nul());
  if(sub_type(Transport_K_KEps,mon_equation.valeur()))
    {
      const Transport_K_KEps& eqn_transport = ref_cast(Transport_K_KEps,mon_equation.valeur());
      const Modele_turbulence_hyd_K_Eps_2_Couches& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_2_Couches,eqn_transport.modele_turbulence());
      const Champ_Fonc& diff_turb = mod_turb.viscosite_turbulente();
      Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem& eval_diff = static_cast<Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem&> (iter->evaluateur());
      eval_diff.associer_diff_turb(diff_turb);
    }
  else if(sub_type(Transport_K_Eps_Bas_Reynolds,mon_equation.valeur()))
    {
      const Transport_K_Eps_Bas_Reynolds& eqn_transport = ref_cast(Transport_K_Eps_Bas_Reynolds,mon_equation.valeur());
      const Modele_turbulence_hyd_K_Eps_Bas_Reynolds& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bas_Reynolds,eqn_transport.modele_turbulence());
      const Champ_Fonc& diff_turb = mod_turb.viscosite_turbulente();
      Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem& eval_diff = static_cast<Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem&> (iter->evaluateur());
      eval_diff.associer_diff_turb(diff_turb);
    }
  else if(sub_type(Transport_K_Eps,mon_equation.valeur()))
    {
      const Transport_K_Eps& eqn_transport = ref_cast(Transport_K_Eps,mon_equation.valeur());
      const Modele_turbulence_hyd_K_Eps& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps,eqn_transport.modele_turbulence());
      const Champ_Fonc& diff_turb = mod_turb.viscosite_turbulente();
      Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem& eval_diffc = static_cast<Eval_Diff_K_Eps_Bas_Re_VDF_const_Elem&> (iter->evaluateur());
      eval_diffc.associer_diff_turb(diff_turb);
    }
}

double Op_Diff_K_Eps_Bas_Re_VDF_var_Elem::calculer_dt_stab() const
{
  double dt_stab, coef;
  const Zone_VDF& zone_VDF = iter->zone();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const DoubleVect& alpha = diffusivite().valeurs(), &alpha_t = diffusivite_turbulente().valeurs();

  // Calcul du pas de temps de stabilite :
  //  - La diffusivite est non constante donc: dt_stab = Min (1/(2*(diff_lam(i)+diff_turb(i))*coeff(elem))
  //     avec :
  //            coeff =  1/(dx*dx) + 1/(dy*dy) + 1/(dz*dz)
  //            i decrivant l'ensemble des elements du maillage

  coef= -1.e10;
  double alpha_local,h_x,h_y,h_z;

  if (dimension == 2)
    {
      int numfa[4];
      for (int elem=0; elem<zone_VDF.nb_elem(); elem++)
        {
          for (int i=0; i<4; i++) numfa[i] = elem_faces(elem,i);

          h_x = zone_VDF.dist_face(numfa[0],numfa[2],0);
          h_y = zone_VDF.dist_face(numfa[1],numfa[3],1);
          alpha_local = (alpha(elem)+alpha_t(elem))*(1/(h_x*h_x) + 1/(h_y*h_y));
          coef = std::max(coef,alpha_local);
        }
    }
  else if (dimension == 3)
    {
      int numfa[6];
      for (int elem=0; elem<zone_VDF.nb_elem(); elem++)
        {
          for (int i=0; i<6; i++) numfa[i] = elem_faces(elem,i);

          h_x = zone_VDF.dist_face(numfa[0],numfa[3],0);
          h_y = zone_VDF.dist_face(numfa[1],numfa[4],1);
          h_z = zone_VDF.dist_face(numfa[2],numfa[5],2);
          alpha_local = (alpha(elem)+alpha_t(elem))*(1/(h_x*h_x) + 1/(h_y*h_y) + 1/(h_z*h_z));
          coef = std::max(coef,alpha_local);
        }
    }

  coef = Process::mp_max(coef);
  dt_stab = 1/(2*(coef+DMINFLOAT));

  return dt_stab;
}

Op_Diff_K_Eps_Bas_Re_VDF_var_Elem::Op_Diff_K_Eps_Bas_Re_VDF_var_Elem() : Op_Diff_K_Eps_Bas_Re_VDF_base(It_VDF_Elem(Eval_Diff_K_Eps_Bas_Re_VDF_var_Elem)()) { }
