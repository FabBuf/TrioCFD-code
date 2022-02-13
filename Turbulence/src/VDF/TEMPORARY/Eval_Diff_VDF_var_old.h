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
// File:        Eval_Diff_VDF_var_old.h
// Directory:   $TRUST_ROOT/src/VDF/Operateurs/Evaluateurs
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Eval_Diff_VDF_var_old_included
#define Eval_Diff_VDF_var_old_included

#include <Eval_Diff_VDF_old.h>
#include <Ref_Champ_base.h>
#include <Champ_base.h>
#include <Champ_Uniforme.h>
#include <Zone_VDF.h>

//
// .DESCRIPTION class Eval_Diff_VDF_var_old
//
// Cette classe represente un evaluateur de flux diffusif
// avec un coefficient de diffusivite qui n'est pas constant en espace.

//.SECTION voir aussi Evaluateur_VDF

class Eval_Diff_VDF_var_old : public Eval_Diff_VDF_old
{

public:

  inline Eval_Diff_VDF_var_old();
  inline void associer(const Champ_base& ) override;
  inline void mettre_a_jour( ) override;
  inline const Champ_base& diffusivite() const;
  inline const Champ_base& get_diffusivite() const override;

protected:

  REF(Champ_base) diffusivite_;
  DoubleVect dv_diffusivite;
  inline double dist_face(int, int, int) const;
  inline double dist_face_period(int, int, int) const;
};

//
//  Fonctions inline de la classe Eval_Diff_VDF_var_old
//
// Description:
// constructeur par defaut
inline Eval_Diff_VDF_var_old::Eval_Diff_VDF_var_old()
{
}

// Description:
// associe le champ de diffusivite
inline void Eval_Diff_VDF_var_old::associer(const Champ_base& diffu)
{
  diffusivite_ = diffu;
  dv_diffusivite.ref(diffu.valeurs());
}

// Description:
// mise a jour de DoubleVect diffusivite
inline void Eval_Diff_VDF_var_old::mettre_a_jour( )
{
  (diffusivite_->valeurs().echange_espace_virtuel());
  dv_diffusivite.ref(diffusivite_->valeurs());
}

// Description:
// renvoie la distance entre les faces fac1 et fac2 dans la direction k
inline double Eval_Diff_VDF_var_old::dist_face(int fac1, int fac2, int k) const
{
  return la_zone->dist_face(fac1, fac2, k);
}

// Description:
// renvoie la distance entre les faces fac1 et fac2 de part et d'autre
// d'une condition aux limites de periodicite dans la direction k
inline double Eval_Diff_VDF_var_old::dist_face_period(int fac1, int fac2, int k) const
{
  return la_zone->dist_face_period(fac1, fac2, k);
}

inline const Champ_base& Eval_Diff_VDF_var_old::diffusivite() const
{
  assert(diffusivite_.non_nul());
  return diffusivite_.valeur();
}

inline const Champ_base&  Eval_Diff_VDF_var_old::get_diffusivite() const
{
  return diffusivite();
}

#endif /* Eval_Diff_VDF_var_old_included */
