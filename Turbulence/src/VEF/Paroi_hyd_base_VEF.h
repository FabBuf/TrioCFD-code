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
// File:        Paroi_hyd_base_VEF.h
// Directory:   $TRUST_ROOT/src/VEF/Turbulence
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////



#ifndef Paroi_hyd_base_VEF_included
#define Paroi_hyd_base_VEF_included

#include <Turbulence_paroi_base.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>


//.DESCRIPTION
//
// CLASS: Paroi_hyd_base_VEF
// Classe de base des lois de paroi hydraulique en VEF

//.SECTION  voir aussi
// Turbulence_paroi_base

class Paroi_hyd_base_VEF : public Turbulence_paroi_base
{

  Declare_base(Paroi_hyd_base_VEF);

public:

  void associer(const Zone_dis& ,const Zone_Cl_dis& );
  void init_lois_paroi_();
  DoubleTab& corriger_derivee_impl(DoubleTab& d) const;
  inline const ArrOfInt& face_keps_imposee() const
  {
    return face_keps_imposee_ ;
  };
  void imprimer_premiere_ligne_ustar(int boundaries_, const LIST(Nom)& boundaries_list, const Nom& nom_fichier_) const;
  void imprimer_ustar_mean_only(Sortie&, int, const LIST(Nom)&, const Nom& ) const;

protected:

  REF(Zone_VEF) la_zone_VEF;
  REF(Zone_Cl_VEF) la_zone_Cl_VEF;
  IntVect face_keps_imposee_; // avec descripteur parallele
  int flag_face_keps_imposee_;
};

#endif
