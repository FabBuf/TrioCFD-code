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
// File:        Eval_Dift_VDF_var_Elem_Axi.cpp
// Directory:   $TRUST_ROOT/src/VDF/Axi/Turbulence
// Version:     /main/4
//
//////////////////////////////////////////////////////////////////////////////

#include <Eval_Dift_VDF_var_Elem_Axi.h>
#include <Paroi_std_scal_hyd_VDF.h>
#include <Turbulence_paroi_scal.h>

void Eval_Dift_VDF_var_Elem_Axi::associer_loipar(const Turbulence_paroi_scal& loi_paroi)
{
  loipar = ref_cast(Paroi_std_scal_hyd_VDF,loi_paroi.valeur());
}


void Eval_Dift_VDF_var_Elem_Axi::mettre_a_jour( )
{
  int s=loipar->tab_equivalent_distance_size();
  equivalent_distance.dimensionner(s);
  for(int i=0; i<s; i++)
    {
      equivalent_distance[i].ref(loipar->tab_equivalent_distance(i));
    }

}
