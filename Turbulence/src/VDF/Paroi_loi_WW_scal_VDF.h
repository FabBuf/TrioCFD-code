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
// File:        Paroi_loi_WW_scal_VDF.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Paroi_loi_WW_scal_VDF_included
#define Paroi_loi_WW_scal_VDF_included

#include <Paroi_scal_hyd_base_VDF.h>
#include <Table.h>

//.DESCRIPTION classe Paroi_loi_WW_scal_VDF
//
//
//

//.SECTION  voir aussi
//  Paroi_loi_WW_scal_VDF

class Paroi_loi_WW_scal_VDF : public Paroi_scal_hyd_base_VDF
{

  Declare_instanciable_sans_constructeur(Paroi_loi_WW_scal_VDF);

public:


  int calculer_scal(Champ_Fonc_base& );
  int init_lois_paroi();


private:

  static double Fthpar(double y_plus,double Pr,double Beta);
  DoubleVect tab_u_star;

};



#endif

