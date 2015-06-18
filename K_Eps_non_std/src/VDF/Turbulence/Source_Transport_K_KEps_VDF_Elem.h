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
// File:        Source_Transport_K_KEps_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Transport_K_KEps_VDF_Elem_included
#define Source_Transport_K_KEps_VDF_Elem_included

#include <Source_Transport_K_Eps_VDF_Elem.h>
#include <Ref_Zone_Cl_VDF.h>
#include <Ref_Transport_K_KEps.h>

class Probleme_base;
class Champ_Don_base;
class DoubleVect;
class DoubleTab;
class Zone_dis;
class Zone_Cl_dis;
class Zone_Cl_VDF;
class Champ_Face;



//.DESCRIPTION class Source_Transport_K_KEps_VDF_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) avec le modele a deux couches et dans le cas
// ou les equations de Navier-Stokes
// ne sont pas couplees a la thermique ou a l'equation de convection-diffusion
// d'une concentration.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_KEps_VDF_Elem : public Source_base,
  public Calcul_Production_K_VDF
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_KEps_VDF_Elem);

public:

  inline Source_Transport_K_KEps_VDF_Elem(double cte1 = C1_DEFAULT,
                                          double cte2 = C2_DEFAULT );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double temps)
  {
    Calcul_Production_K_VDF::mettre_a_jour(temps);
  }

protected:

  double C1;
  double C2;
  REF(Zone_VDF) la_zone_VDF;
  REF(Zone_Cl_VDF) la_zone_Cl_VDF;
  REF(Equation_base) eq_hydraulique;
  REF(Transport_K_KEps)  mon_eq_transport_K_Eps;

  virtual void associer_pb(const Probleme_base& pb);
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );

};


//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Eps_anisotherme_VDF_Elem
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) pour le modele a deux couches et dans le cas
// ou les equations de Navier_Stokes
// sont couplees a l'equation de la thermique
// On suppose que le coefficient de variation de la masse volumique
// du fluide en fonction de ce scalaire est un coefficient uniforme.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_KEps_anisotherme_VDF_Elem :
  public Source_Transport_K_KEps_VDF_Elem
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_KEps_anisotherme_VDF_Elem);

public:

  inline Source_Transport_K_KEps_anisotherme_VDF_Elem(double cte1 = C1_DEFAULT,
                                                      double cte2 = C2_DEFAULT,
                                                      double cte3 = C3_DEFAULT);
  virtual void associer_pb(const Probleme_base& );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;

protected:

  double C3;
  REF(Convection_Diffusion_Temperature) eq_thermique;
  REF(Champ_Don) beta_t;
  REF(Champ_Don_base) gravite;

};



//////////////////////////////////////////////////////////////////////////////
//
//   Fonctions inline de la classe Source_Transport_K_KEps_VDF_Elem
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_KEps_VDF_Elem::
Source_Transport_K_KEps_VDF_Elem(double cte1,double cte2)

  : C1(cte1), C2(cte2) {}


//////////////////////////////////////////////////////////////////////////////
//
//                     Fonctions inline de la classe
//
//             Source_Transport_K_KEps_anisotherme_VDF_Elem
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_KEps_anisotherme_VDF_Elem::
Source_Transport_K_KEps_anisotherme_VDF_Elem(double cte1,double cte2,double cte3)

  : Source_Transport_K_KEps_VDF_Elem(cte1,cte2) , C3(cte3) {}


#endif
