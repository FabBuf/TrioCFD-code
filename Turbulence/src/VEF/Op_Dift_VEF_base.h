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
// File:        Op_Dift_VEF_base.h
// Directory:   $TRUST_ROOT/src/VEF/Turbulence
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////



#ifndef Op_Dift_VEF_base_included
#define Op_Dift_VEF_base_included

#include <Op_Diff_VEF_base.h>
#include <Op_Diff_Turbulent_base.h>
#include <Ref_Mod_turb_hyd_base.h>
#include <Ref_Modele_turbulence_scal_base.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION class  Op_Dift_VEF_base
//
//////////////////////////////////////////////////////////////////////////////

class Op_Dift_VEF_base : public Op_Diff_VEF_base, public Op_Diff_Turbulent_base
{


  Declare_base(Op_Dift_VEF_base);

public:

  virtual void calculer_borne_locale(DoubleVect& ,double ,double ) const;
  void associer_modele_turbulence(const Mod_turb_hyd_base& );
  void associer_modele_turbulence_temp(const Modele_turbulence_scal_base& );
  void mettre_a_jour(double temps);
  void associer(const Zone_dis& , const Zone_Cl_dis& , const Champ_Inc& );
  void completer();

protected:
  REF(Mod_turb_hyd_base) le_modele_turbulence; // A deplacer dans Op_Diff_turb ?
  REF(Modele_turbulence_scal_base) le_modele_turb_temp;  // A deplacer dans Op_Diff_turb ?
  DoubleTab tau_tan_;
  DoubleTab k;
  int indic_lp_neg_;
  int indic_bas_Re_;
  int indic_Pr;  // mod turb temp de Prandtl Oui = 1
};

#endif
