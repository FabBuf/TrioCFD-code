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
// File:        MuLambda_TBLE.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////

#include <MuLambda_TBLE.h>


Implemente_deriv(MuLambda_TBLE_base);
Implemente_instanciable(MuLambda_TBLE,"MuLambda_TBLE",DERIV(MuLambda_TBLE_base));

// Description:
//
// Precondition:
// Parametre: Sortie& os
//    Signification: le flot de sortie de sauvegarde
//    Valeurs par defaut:
//    Contraintes:
//    Acces: sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord: le flot de sortie est modifie
// Postcondition: la methode ne modifie pas l'objet

Sortie& MuLambda_TBLE::printOn(Sortie& os) const
{
  return DERIV(MuLambda_TBLE_base)::printOn(os);
}


// Description:
//
// Precondition:
// Parametre: Entree& is
//    Signification: le flot d'entree pour la lecture d'une equation
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord: le flot d'entree est modifie
// Postcondition: l'equation est construite avec les elements lus.
Entree& MuLambda_TBLE::readOn( Entree& is)
{
  return DERIV(MuLambda_TBLE_base)::readOn(is);
}


