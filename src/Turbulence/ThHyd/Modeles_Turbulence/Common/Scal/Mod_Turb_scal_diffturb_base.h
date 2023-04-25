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
// File:        Mod_Turb_scal_diffturb_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/Common/Scal
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Mod_Turb_scal_diffturb_base_included
#define Mod_Turb_scal_diffturb_base_included

#include <Modele_turbulence_scal_base.h>
#include <Probleme_base.h>
#include <Milieu_base.h>
#include <TRUST_Ref.h>

class Champ_Fonc;

/*! @brief Classe Mod_Turb_scal_diffuturb_base Cette classe represente la classe de base pour le modele de calcul suivant
 *
 *     pour la diffusivite turbulente:
 *                   alpha_turb = visco_turb / coeff_turb;
 *
 * @sa Modele_turbulence_scal_base
 */
class Mod_Turb_scal_diffturb_base : public Modele_turbulence_scal_base
{

  Declare_base(Mod_Turb_scal_diffturb_base);

public:

  //virtual int comprend_mot(const Motcle& ) const;
  void associer_viscosite_turbulente(const Champ_Fonc& );
  void completer() override;
  int reprendre(Entree& ) override;

protected:

  REF(Champ_Fonc) la_viscosite_turbulente;
};

#endif
