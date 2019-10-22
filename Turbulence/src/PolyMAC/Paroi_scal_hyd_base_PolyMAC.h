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
// File:        Paroi_scal_hyd_base_PolyMAC.h
// Directory:   $TRUST_ROOT/src/PolyMAC/Turbulence
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Paroi_scal_hyd_base_PolyMAC_included
#define Paroi_scal_hyd_base_PolyMAC_included

#include <Turbulence_paroi_scal_base.h>
#include <Zone_PolyMAC.h>
#include <Ref_Zone_PolyMAC.h>
#include <Ref_Zone_Cl_PolyMAC.h>

//.DESCRIPTION classe Paroi_scal_hyd_base_PolyMAC
//
//
//

//.SECTION  voir aussi
//  Paroi_std_hyd_PolyMAC

class Paroi_scal_hyd_base_PolyMAC : public Turbulence_paroi_scal_base
{

  Declare_base(Paroi_scal_hyd_base_PolyMAC);

public:
  void associer(const Zone_dis& ,const Zone_Cl_dis& );
  virtual int init_lois_paroi();
  virtual void imprimer_nusselt(Sortie&) const;
  virtual DoubleVect& equivalent_distance_name(DoubleVect& d_equiv, const Nom& nom_bord) const;

protected:

  REF(Zone_PolyMAC) la_zone_PolyMAC;
  REF(Zone_Cl_PolyMAC) la_zone_Cl_PolyMAC;

  mutable int nb_impr_;        // Compteur d'impression

};

#endif

