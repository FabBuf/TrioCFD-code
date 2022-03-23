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
// File:        Mod_echelle_LRM.h
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Mod_echelle_LRM_included
#define Mod_echelle_LRM_included

#include <Deriv.h>
#include <Mod_echelle_LRM_base.h>

Declare_deriv(Mod_echelle_LRM_base);

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Mod_echelle_LRM
//    Classe generique calculant les echelle de vitesse et de longueur pour le
//        modele bas-Reynolds  dans les equations de
//    couche limite simplifiees necessaires
//    a l'utilisation des lois de parois de type TBLE_LRM
//
// .SECTION voir aussi
//    Mod_echelle_LRM_base
//
//////////////////////////////////////////////////////////////////////////////

class Mod_echelle_LRM : public DERIV(Mod_echelle_LRM_base)

{
  Declare_instanciable(Mod_echelle_LRM);

public :

};


#endif
