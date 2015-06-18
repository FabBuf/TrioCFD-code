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
// File:        Source_DWF_VDF.h
// Directory:   $TRUST_ROOT/src/Zoom/VDF/Turbulence
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_DWF_VDF_included
#define Source_DWF_VDF_included



#include <Source_base.h>
#include <Zone_VDF.h>

class Probleme_base;
class Zone_dis;
class Zone_Cl_dis;

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_DWF_VDF
// Cette classe modelise les termes sources pour l'equation fine des
// lois de paroi DWF.
//
//////////////////////////////////////////////////////////////////////////////

class Source_DWF_VDF : public Source_base
{

  Declare_instanciable(Source_DWF_VDF);

public:

  inline void associer_zones(const Zone_dis&, const Zone_Cl_dis& );

  inline DoubleTab& ajouter(DoubleTab& ) const ;
  inline DoubleTab& calculer(DoubleTab& ) const ;
  void completer();
  inline void associer_pb(const Probleme_base& ) { };

  inline  DoubleTab& getValeurs()
  {
    return les_valeurs ;
  }

  void mettre_a_jour(double temps)
  {
    ;
  }

private:
  DoubleTab les_valeurs;


};

inline DoubleTab& Source_DWF_VDF::ajouter(DoubleTab& resu) const
{
  resu += les_valeurs;
  return resu;
}

inline DoubleTab& Source_DWF_VDF::calculer(DoubleTab& resu) const
{
  resu = 0;
  resu = les_valeurs;
  return resu;
}

void Source_DWF_VDF::associer_zones(const Zone_dis& zone_dis,
                                    const Zone_Cl_dis& zone_cl_dis)
{
  const Zone_VDF& zvdf = ref_cast(Zone_VDF,zone_dis.valeur());

  int nb_face = zvdf.nb_faces();
  les_valeurs.resize(nb_face);

}

#endif
