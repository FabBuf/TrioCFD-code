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
// File:        Source_Transport_K_Eps_VEF_Face.h
// Directory:   $TRUST_ROOT/src/VEF/Turbulence
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Transport_K_Eps_VEF_Face_included
#define Source_Transport_K_Eps_VEF_Face_included

#define C1_DEFAULT 1.44   // Valeurs par defaut des constantes qui interviennent
#define C2_DEFAULT 1.92   // dans le calcul des termes sources des equations
#define C3_DEFAULT 1.0    // de transport de K et Eps source: Chabard de N3S

#include <Source_base.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Champ_Don.h>
#include <Ref_Champ_Don_base.h>
#include <Ref_Convection_Diffusion_Temperature.h>
#include <Ref_Convection_Diffusion_Concentration.h>
#include <Ref_Equation_base.h>
#include <Ref_Transport_K_Eps.h>
#include <Calcul_Production_K_VEF.h>

class Probleme_base;
class Champ_Don_base;
class DoubleVect;
class DoubleTab;
class Zone_dis;
class Zone_Cl_dis;
class Zone_Cl_VEF;
class Champ_Face;
class Champ_base;
//
//  Remarque:  On n'a pas integre les classes qui suivent dans la hierarchie
//             Terme_Source_VEF_base car il etait plus simple de les traiter
//             hors hierarchie



//.DESCRIPTION class Source_Transport_K_Eps_VEF_Face
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de Navier-Stokes
// ne sont pas couplees a la thermique ou a l'equation de convection-diffusion
// d'une concentration.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Eps_VEF_Face : public Source_base,
  public Calcul_Production_K_VEF
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_VEF_Face);

public:

  inline Source_Transport_K_Eps_VEF_Face(double cte1 = C1_DEFAULT,
                                         double cte2 = C2_DEFAULT );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const ;
  void mettre_a_jour(double temps)
  {
    Calcul_Production_K_VEF::mettre_a_jour(temps);
  }

protected:

  double C1;
  double C2;
  REF(Zone_VEF) la_zone_VEF;
  REF(Equation_base) eq_hydraulique;
  REF(Transport_K_Eps)  mon_eq_transport_K_Eps;

  virtual void associer_pb(const Probleme_base& pb);
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
};

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Eps_anisotherme_VEF_Face
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de Navier_Stokes
// sont couplees a l'equation de la thermique
// On suppose que le coefficient de variation de la masse volumique
// du fluide en fonction de ce scalaire est un coefficient uniforme.
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Eps_anisotherme_VEF_Face :
  public Source_Transport_K_Eps_VEF_Face
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_anisotherme_VEF_Face);

public:

  inline Source_Transport_K_Eps_anisotherme_VEF_Face(double cte1 = C1_DEFAULT,
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
// CLASS: Source_Transport_K_Eps_aniso_concen_VEF_Face
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de
// Navier_Stokes sont couplees a l'equation de convection diffusion
// d'une concentration
// Le champ beta_c est uniforme
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Eps_aniso_concen_VEF_Face :
  public Source_Transport_K_Eps_VEF_Face
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_aniso_concen_VEF_Face);

public:

  inline Source_Transport_K_Eps_aniso_concen_VEF_Face(double cte1 = C1_DEFAULT,
                                                      double cte2 = C2_DEFAULT,
                                                      double cte3 = C3_DEFAULT);
  virtual void associer_pb(const Probleme_base& );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;

protected:

  double C3;
  REF(Convection_Diffusion_Concentration) eq_concentration;
  REF(Champ_Don) beta_c;
  REF(Champ_Don_base) gravite;

};

//////////////////////////////////////////////////////////////////////////////
//
// CLASS: Source_Transport_K_Eps_aniso_therm_concen_VEF_Face
//
// Cette classe represente le terme source qui figure dans l'equation
// de transport du couple (k,eps) dans le cas ou les equations de
// Navier_Stokes sont couplees a l'equation de convection diffusion
// d'une concentration et a l'equation de la thermique
// Les champs beta_t et beta_c sont uniformes
//
//////////////////////////////////////////////////////////////////////////////

class Source_Transport_K_Eps_aniso_therm_concen_VEF_Face :
  public Source_Transport_K_Eps_VEF_Face
{

  Declare_instanciable_sans_constructeur(Source_Transport_K_Eps_aniso_therm_concen_VEF_Face);

public:

  inline Source_Transport_K_Eps_aniso_therm_concen_VEF_Face(double cte1 = C1_DEFAULT,
                                                            double cte2 = C2_DEFAULT,
                                                            double cte3 = C3_DEFAULT);
  virtual void associer_pb(const Probleme_base& );
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;

protected:

  double C3;
  REF(Convection_Diffusion_Temperature) eq_thermique;
  REF(Convection_Diffusion_Concentration) eq_concentration;
  REF(Champ_Don) beta_t;
  REF(Champ_Don) beta_c;
  REF(Champ_Don_base) gravite;

};

//////////////////////////////////////////////////////////////////////////////
//
//   Fonctions inline de la classe Source_Transport_K_Eps_VEF_Face
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_VEF_Face::
Source_Transport_K_Eps_VEF_Face(double cte1,double cte2)

  : C1(cte1), C2(cte2) {}

//////////////////////////////////////////////////////////////////////////////
//
//                     Fonctions inline de la classe
//
//             Source_Transport_K_Eps_anisotherme_VEF_Face
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_anisotherme_VEF_Face::
Source_Transport_K_Eps_anisotherme_VEF_Face(double cte1,double cte2,double cte3)

  : Source_Transport_K_Eps_VEF_Face(cte1,cte2) , C3(cte3) {}

//////////////////////////////////////////////////////////////////////////////
//
//                        Fonctions inline de la classe
//
//                Source_Transport_K_Eps_aniso_concen_VEF_Face
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_aniso_concen_VEF_Face::
Source_Transport_K_Eps_aniso_concen_VEF_Face(double cte1,double cte2,double cte3)

  : Source_Transport_K_Eps_VEF_Face(cte1,cte2) , C3(cte3) {}


//////////////////////////////////////////////////////////////////////////////
//
//                        Fonctions inline de la classe
//
//                Source_Transport_K_Eps_aniso_therm_concen_VEF_Face
//
//////////////////////////////////////////////////////////////////////////////

inline Source_Transport_K_Eps_aniso_therm_concen_VEF_Face::
Source_Transport_K_Eps_aniso_therm_concen_VEF_Face(double cte1,double cte2,double cte3)

  : Source_Transport_K_Eps_VEF_Face(cte1,cte2) , C3(cte3) {}

#endif
