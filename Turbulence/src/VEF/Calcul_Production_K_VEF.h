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
// File:        Calcul_Production_K_VEF.h
// Directory:   $TRUST_ROOT/src/VEF/Turbulence
// Version:     /main/2
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Calcul_Production_K_VEF_included
#define Calcul_Production_K_VEF_included

#include <Champ_Don.h>

/////////////////////////////////////////////////////////////////////////
//
//  CLASS Calcul_Production_K_VEF
//
// Classe qui porte les fonctions de calcul des termes de production
// de l'energie cinetique turbulente et de destruction de cette energie
// Cette classe ne derive d'aucune autre car nous voulons l'utiliser
// pour faire de l'heritage multiple.
//
/////////////////////////////////////////////////////////////////////////
class DoubleTab;
class Zone_VEF;
class Zone_Cl_VEF;
class DoubleVect;

class Calcul_Production_K_VEF
{
protected:
  Calcul_Production_K_VEF() {}

  DoubleTab&  calculer_terme_production_K(const Zone_VEF& ,const Zone_Cl_VEF& ,
                                          DoubleTab& ,const DoubleTab& ,
                                          const DoubleTab& ,const DoubleTab& ) const;

  DoubleTab&  calculer_terme_production_K_EASM(const Zone_VEF& ,const Zone_Cl_VEF& ,
                                               DoubleTab& ,const DoubleTab& ,
                                               const DoubleTab& , const DoubleTab&,
                                               const DoubleTab&) const;

  DoubleTab& calcul_tenseur_face(DoubleTab&, const DoubleTab&,
                                 const Zone_VEF&, const Zone_Cl_VEF&) const;

  DoubleTab& calculer_terme_destruction_K_gen(const Zone_VEF& , const Zone_Cl_VEF& , DoubleTab& ,
                                              const DoubleTab& , const DoubleTab& ,
                                              const Champ_Don& ,const DoubleVect& ,int ) const;

  void mettre_a_jour(double temps)
  {
    ;
  };

};

#endif
