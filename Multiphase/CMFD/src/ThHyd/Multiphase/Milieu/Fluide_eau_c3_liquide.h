/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Fluide_eau_c3_liquide.h
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Milieu
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Fluide_eau_c3_liquide_included
#define Fluide_eau_c3_liquide_included

#include <Fluide_reel_base.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Fluide_eau_c3_liquide
//    Cette classe represente un milieu reel
//    dont les lois viennent de "Lois_Na"
//////////////////////////////////////////////////////////////////////////////
class Fluide_eau_c3_liquide: public Fluide_reel_base
{
  Declare_instanciable(Fluide_eau_c3_liquide);

  /* bornes (p, T_l) de H2OPROP.H */
  std::map<std::string, std::array<double, 2>> unknown_range() const override
  {
    return { { "pression"   , { 0.01e5, 260e5 } },
      { "temperature", {      0, 400   } }
    };
  }

protected :
  // densite
  double     rho_(const double T, const double P) const override;
  double  dP_rho_(const double T, const double P) const override;
  double  dT_rho_(const double T, const double P) const override;
  // enthalpie
  double       h_(const double T, const double P) const override;
  double    dP_h_(const double T, const double P) const override;
  double    dT_h_(const double T, const double P) const override;
  // lois champs "faibles" -> pas de derivees
  double      cp_(const double T, const double P) const override;
  double    beta_(const double T, const double P) const override;
  double      mu_(const double T, const double P) const override;
  double  lambda_(const double T, const double P) const override;
};

#endif
