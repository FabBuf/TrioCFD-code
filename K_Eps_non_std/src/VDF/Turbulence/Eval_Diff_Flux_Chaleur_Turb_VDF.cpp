/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Eval_Diff_Flux_Chaleur_Turb_VDF.cpp
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/5
//
//////////////////////////////////////////////////////////////////////////////

#include <Eval_Diff_Flux_Chaleur_Turb_VDF.h>
#include <Champ_Fonc.h>
#include <Zone_VDF.h>


void Eval_Diff_Flux_Chaleur_Turb_VDF::associer_diff_turb(const Champ_Fonc& diffu)
{
  diffusivite_turbulente_ = diffu;
}

void Eval_Diff_Flux_Chaleur_Turb_VDF::mettre_a_jour( )
{
  dv_diffusivite_turbulente.ref(diffusivite_turbulente_->valeurs());
}

// Description:
// renvoie la distance entre les faces fac1 et fac2 dans la direction k
double Eval_Diff_Flux_Chaleur_Turb_VDF::dist_face(int fac1, int fac2, int k) const
{
  return la_zone->dist_face(fac1, fac2, k);
}

double Eval_Diff_Flux_Chaleur_Turb_VDF::dist_face_period(int fac1, int fac2, int k) const
{
  return la_zone->dist_face_period(fac1,fac2,k);
}
