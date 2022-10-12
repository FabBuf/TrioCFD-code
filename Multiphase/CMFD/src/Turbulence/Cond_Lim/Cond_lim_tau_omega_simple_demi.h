/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Cond_lim_tau_omega_simple_demi.h
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Cond_lim_tau_omega_simple_demi_included
#define Cond_lim_tau_omega_simple_demi_included

#include <TRUSTTab.h>
#include <Echange_global_impose.h>
#include <Cond_lim_base.h>
#include <Ref_Correlation.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    Classe Cond_lim_tau_omega_simple_demi
// .SECTION voir aussi
//////////////////////////////////////////////////////////////////////////////
class Cond_lim_tau_omega_simple_demi : public Echange_global_impose
{

  Declare_instanciable(Cond_lim_tau_omega_simple_demi);

public :
  int compatible_avec_eqn(const Equation_base&) const override;
  virtual int initialiser(double temps) override;
  virtual int avancer(double temps) override {return 1;}; // Avancer ne fait rien car le champ est modifie dans mettre_a_jour
  void mettre_a_jour(double tps) override;
  double calc_tau(double y, double u_tau, double visc);
  double calc_omega(double y, double u_tau, double visc);
  virtual void liste_faces_loi_paroi(IntTab&) override;
  virtual void completer() override;


  void associer_fr_dis_base(const Frontiere_dis_base& fr) override {la_frontiere_dis=fr;};
  void associer_zone_cl_dis_base(const Zone_Cl_dis_base& zcl)  override { ma_zone_cl_dis=zcl;};

  // fonctions de cond_lim_base qui necessitent le champ_front qu'on met a zero car on fait abstraction du champ_front
  void fixer_nb_valeurs_temporelles(int nb_cases) override {};
  inline Frontiere_dis_base& frontiere_dis() override {return la_frontiere_dis;};
  inline const Frontiere_dis_base& frontiere_dis() const override {return la_frontiere_dis;};
  void changer_temps_futur(double temps,int i) override {};
  void set_temps_defaut(double temps) override {};
  void calculer_coeffs_echange(double temps) override {};
  void verifie_ch_init_nb_comp() const override {};

  Champ_front& T_ext() override {Process::exit(que_suis_je() + " : You shouldn't go through T_ext ! ") ; return Echange_impose_base::T_ext();};
  const Champ_front& T_ext() const override {Process::exit(que_suis_je() + " : You shouldn't go through T_ext ! ") ; return Echange_impose_base::T_ext();};
  inline virtual Champ_front& h_imp() override {Process::exit(que_suis_je() + " : You shouldn't go through h_imp ! ") ; return Echange_impose_base::h_imp();};
  inline virtual const Champ_front& h_imp() const override {Process::exit(que_suis_je() + " : You shouldn't go through h_imp ! ") ; return Echange_impose_base::h_imp();};
  double h_imp(int num) const override ;
  double h_imp(int num,int k) const override;
  double h_imp_grad(int num) const override ;
  double h_imp_grad(int num,int k) const override;
  double T_ext(int num) const override;
  double T_ext(int num,int k) const override;

protected :
  void me_calculer();

  REF(Correlation) correlation_loi_paroi_;
  REF(Frontiere_dis_base) la_frontiere_dis;

  double mon_temps = -1.e8;

  DoubleTab h_;
  DoubleTab h_grad_;
  DoubleTab d_;
  double von_karman_ = 0.41 ;
  double beta_omega = 0.075;
  double beta_k = 0.09;
  double is_tau_=-1 ; // 0 : omega ; 1 : tau
};

#endif
