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
// File:        Pb_Phase_field.h
// Directory:   $TRUST_ROOT/src/Phase_field
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Pb_Phase_field_included
#define Pb_Phase_field_included

#include <Pb_Hydraulique.h>
#include <Convection_Diffusion_Phase_field.h>
#include <Navier_Stokes_phase_field.h>



//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Pb_Phase_field
//    Cette classe represente un probleme d'hydraulique avec transport
//    d'un ou plusieurs constituants:
//       - Equations de Navier_Stokes en regime laminaire
//         pour un fluide incompressible
//       - Equations de convection-diffusion en regime laminaire
//         En fait si on transporte plusieurs constituants on utilisera une
//         seule equation de convection-diffusion avec une inconnue vectorielle.
//         En general, on couple les 2 equations par l'intermediaire du terme
//         source des forces de volume de Navier_Stokes dans lequel on prend
//         en compte de petites variations de la masse volumique en fonction
//         du ou des constituants
// .SECTION voir aussi
//     Pb_qdm_fluide
//////////////////////////////////////////////////////////////////////////////
class Pb_Phase_field : public Pb_qdm_fluide
{

  Declare_instanciable(Pb_Phase_field);

public:

  int nombre_d_equations() const;
  const Equation_base& equation(int) const ;
  Equation_base& equation(int);
  void associer_milieu_base(const Milieu_base& );
  int verifier();

protected:

  Navier_Stokes_phase_field eq_hydraulique;
  Convection_Diffusion_Phase_field eq_concentration;

};



#endif
