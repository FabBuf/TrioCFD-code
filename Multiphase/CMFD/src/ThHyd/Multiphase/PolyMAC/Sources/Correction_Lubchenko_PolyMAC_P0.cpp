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
// File:        Correction_Lubchenko_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Correction_Lubchenko_PolyMAC_P0.h>
#include <Pb_Multiphase.h>
#include <Portance_interfaciale_PolyMAC_P0.h>
#include <Portance_interfaciale_base.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Dispersion_bulles_PolyMAC_P0.h>
#include <Dispersion_bulles_base.h>
#include <Champ_Face_PolyMAC_P0.h>
#include <Milieu_composite.h>
#include <Viscosite_turbulente_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <math.h>

Implemente_instanciable(Correction_Lubchenko_PolyMAC_P0, "Correction_Lubchenko_Face_PolyMAC_P0", Source_base);

Sortie& Correction_Lubchenko_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Correction_Lubchenko_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_lift", &beta_lift_);
  param.ajouter("beta_disp", &beta_disp_);
  param.lire_avec_accolades_depuis(is);

  //identification des phases
  Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : NULL;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;

  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  pbm->creer_champ("distance_paroi_globale"); // Besoin de distance a la paroi

  return is;
}

void Correction_Lubchenko_PolyMAC_P0::completer() // We must wait for all readOn's to be sure that the bubble dispersion and lift correlations are created
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());

  for (int i = 0 ; i < equation().sources().size() ; i++)
    if sub_type(Portance_interfaciale_PolyMAC_P0, equation().sources()(i).valeur()) correlation_lift_ = ref_cast(Portance_interfaciale_PolyMAC_P0, equation().sources()(i).valeur()).correlation();

  for (int i = 0 ; i < equation().sources().size() ; i++)
    if sub_type(Dispersion_bulles_PolyMAC_P0, equation().sources()(i).valeur()) correlation_dispersion_ = ref_cast(Dispersion_bulles_PolyMAC_P0, equation().sources()(i).valeur()).correlation();

  if ( (!correlation_lift_.non_nul()) || (!correlation_dispersion_.non_nul()) ) Process::exit("Correction_Lubchenko_PolyMAC_P0::completer() : a dispersion_bulles and a portance_interfaciale force must be defined !");

  if sub_type(Op_Diff_Turbulent_PolyMAC_P0_Face, equation().operateur(0).l_op_base()) is_turb = 1;

  if (!pbm.has_champ("diametre_bulles")) Process::exit("Correction_Lubchenko_PolyMAC_P0::completer() : a bubble diameter must be defined !");
}


void Correction_Lubchenko_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const // The necessary dimensionner_bloc is taken care of in the dispersion_bulles_PolyMAC_P0 and Portance_interfaciale_PolyMAC_P0 functions
{
}

void Correction_Lubchenko_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  ajouter_blocs_disp(matrices, secmem, semi_impl);
  ajouter_blocs_lift(matrices, secmem, semi_impl);
}

void Correction_Lubchenko_PolyMAC_P0::ajouter_blocs_disp(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl ) const
{
  const Champ_Face_PolyMAC_P0& ch = ref_cast(Champ_Face_PolyMAC_P0, equation().inconnue().valeur());
  const Zone_PolyMAC_P0& zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const IntTab& f_e = zone.face_voisins(), &fcl = ch.fcl(), &e_f = zone.elem_faces();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &pf = equation().milieu().porosite_face(), &ve = zone.volumes(), &vf = zone.volumes_entrelaces(), &fs = zone.face_surfaces();
  const DoubleTab& vf_dir = zone.volumes_entrelaces_dir(), &n_f = zone.face_normales();
  const DoubleTab& pvit = ch.passe(),
                   &alpha = ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().passe(),
                    &press = ref_cast(Pb_Multiphase, equation().probleme()).eq_qdm.pression().passe(),
                     &temp  = ref_cast(Pb_Multiphase, equation().probleme()).eq_energie.inconnue().passe(),
                      &rho   = equation().milieu().masse_volumique().passe(),
                       &mu    = ref_cast(Fluide_base, equation().milieu()).viscosite_dynamique().passe(),
                        &y_elem = zone.y_elem(),
                         &y_faces = zone.y_faces(),
                          &n_y_elem = zone.normale_paroi_elem(),
                           &n_y_faces = zone.normale_paroi_faces(),
                            &d_bulles = equation().probleme().get_champ("diametre_bulles").valeurs(),
                             *k_turb = (equation().probleme().has_champ("k")) ? &equation().probleme().get_champ("k").passe() : NULL ;
  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());

  int N = pvit.line_size() , Np = press.line_size(), Nk = (k_turb) ? (*k_turb).dimension(1) : 1, D = dimension,
      nf_tot = zone.nb_faces_tot(), nf = zone.nb_faces(), ne_tot = zone.nb_elem_tot(),
      cR = (rho.dimension_tot(0) == 1), cM = (mu.dimension_tot(0) == 1);
  DoubleTrav a_l(N), p_l(N), T_l(N), rho_l(N), mu_l(N), sigma_l(N,N), dv(N, N), nut_l(N), k_l(Nk), d_b_l(N), coeff(N, N, 2); //arguments pour coeff

  DoubleTrav nut(zone.nb_elem_tot(), N); //viscosite turbulente
  if (is_turb) ref_cast(Viscosite_turbulente_base, ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, equation().operateur(0).l_op_base()).correlation().valeur()).eddy_viscosity(nut); //remplissage par la correlation

  const Dispersion_bulles_base& correlation_db = ref_cast(Dispersion_bulles_base, correlation_dispersion_->valeur());

  // There is no need to calculate the gradient of alpha here

  /* faces */
  for (int f = 0; f < nf; f++)
    if (fcl(f, 0) < 2)
      {
        a_l = 0 ;
        p_l = 0 ;
        T_l = 0 ;
        rho_l = 0;
        mu_l = 0 ;
        d_b_l = 0;
        nut_l = 0 ;
        k_l = 0 ;
        sigma_l=0;
        dv = 0;
        int e;
        for (int c = 0; c < 2 && (e = f_e(f, c)) >= 0; c++)
          {
            for (int n = 0; n < N; n++)
              {
                a_l(n)   += vf_dir(f, c)/vf(f) * alpha(e, n);
                p_l(n)   += vf_dir(f, c)/vf(f) * press(e, n * (Np > 1));
                T_l(n)   += vf_dir(f, c)/vf(f) * temp(e, n);
                rho_l(n) += vf_dir(f, c)/vf(f) * rho(!cR * e, n);
                mu_l(n)  += vf_dir(f, c)/vf(f) * mu(!cM * e, n);
                nut_l(n) += is_turb    ? vf_dir(f, c)/vf(f) * nut(e,n) : 0;
                d_b_l(n) += vf_dir(f, c)/vf(f) * d_bulles(e,n) ;
                for (int k = 0; k < N; k++)
                  if (milc.has_interface(n,k))
                    {
                      Interface_base& sat = milc.get_interface(n, k);
                      sigma_l(n,k) += vf_dir(f, c)/vf(f) * sat.sigma(temp(e,n),press(e,n * (Np > 1)));
                    }
                for (int k = 0; k < N; k++)
                  dv(k, n) += vf_dir(f, c)/vf(f) * ch.v_norm(pvit, pvit, e, f, k, n, nullptr, nullptr);
              }
            for (int n = 0; n <Nk; n++) k_l(n)   += (k_turb)   ? vf_dir(f, c)/vf(f) * (*k_turb)(e,0) : 0;
          }

        correlation_db.coefficient(a_l, p_l, T_l, rho_l, mu_l, sigma_l, nut_l, k_l, d_b_l, dv, coeff);

        double sum_alphag_wall = 0 ;
        for (int k = 0; k<N ; k++)
          if (k!=n_l) sum_alphag_wall += (y_faces(f)<d_b_l(k)/2.) ? a_l(k) * (d_b_l(k)-2*y_faces(f))/(d_b_l(k)-y_faces(f)) :0 ;

        for (int k = 0; k < N; k++)
          if (k != n_l)
            if (y_faces(f)<d_b_l(k)/2.)
              {
                double fac = 0 ;
                for (int d = 0 ; d<D ; d++) fac += n_y_faces(f, d) * n_f(f, d)/fs(f);

                fac *= beta_disp_*pf(f) * vf(f) ;
                secmem(f, k)   += fac * coeff(k, n_l, 0) * 1/y_faces(f) * a_l(k) * (d_b_l(k)-2*y_faces(f))/(d_b_l(k)-y_faces(f));
                secmem(f, k)   += fac * coeff(n_l, k, 0) * 1/y_faces(f) * sum_alphag_wall;
                secmem(f, n_l) -= fac * coeff(k, n_l, 0) * 1/y_faces(f) * a_l(k) * (d_b_l(k)-2*y_faces(f))/(d_b_l(k)-y_faces(f));
                secmem(f, n_l) -= fac * coeff(n_l, k, 0) * 1/y_faces(f) * sum_alphag_wall;
              }
      }

  /* elements */
  for (int e = 0; e < ne_tot; e++)
    {
      /* arguments de coeff */
      for (int n = 0; n < N; n++)
        {
          a_l(n)   = alpha(e, n);
          p_l(n)   = press(e, n * (Np > 1));
          T_l(n)   =  temp(e, n);
          rho_l(n) =   rho(!cR * e, n);
          mu_l(n)  =    mu(!cM * e, n);
          nut_l(n) = is_turb    ? nut(e,n) : 0;
          d_b_l(n) = d_bulles(e,n) ;
          for (int k = 0; k < N; k++)
            if (milc.has_interface(n,k))
              {
                Interface_base& sat = milc.get_interface(n, k);
                sigma_l(n,k) = sat.sigma(temp(e,n),press(e,n * (Np > 1)));
              }

          for (int k = 0; k < N; k++)
            dv(k, n) = ch.v_norm(pvit, pvit, e, -1, k, n, nullptr, nullptr);
        }

      for (int n = 0; n <Nk; n++) k_l(n)   = (k_turb)   ? (*k_turb)(e,0) : 0;

      correlation_db.coefficient(a_l, p_l, T_l, rho_l, mu_l, sigma_l, nut_l, k_l, d_b_l, dv, coeff);

      double sum_alphag_wall = 0 ;
      for (int k = 0; k<N ; k++)
        if (k!=n_l) sum_alphag_wall += (y_elem(e)<d_b_l(k)/2.) ? a_l(k) *(d_b_l(k)-2*y_elem(e))/(d_b_l(k)-y_elem(e)) :0 ;
      for (int d = 0, i = nf_tot + D * e; d < D; d++, i++)
        for (int k = 0; k < N; k++)
          if (k != n_l)
            if (y_elem(e)<d_b_l(k)/2)
              {
                double fac = beta_disp_*pe(e) * ve(e);

                secmem(i, k)   += fac * coeff(k, n_l, 0) * 1/y_elem(e) * a_l(k) * (d_b_l(k)-2*y_elem(e))/(d_b_l(k)-y_elem(e)) * n_y_elem(e, d);
                secmem(i, k)   += fac * coeff(n_l, k, 0) * 1/y_elem(e) * sum_alphag_wall                                      * n_y_elem(e, d);
                secmem(i, n_l) -= fac * coeff(k, n_l, 0) * 1/y_elem(e) * a_l(k) * (d_b_l(k)-2*y_elem(e))/(d_b_l(k)-y_elem(e)) * n_y_elem(e, d);
                secmem(i, n_l) -= fac * coeff(n_l, k, 0) * 1/y_elem(e) * sum_alphag_wall                                      * n_y_elem(e, d);
              }
    }
}


void Correction_Lubchenko_PolyMAC_P0::ajouter_blocs_lift(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl ) const
{
  const Champ_Face_PolyMAC_P0& ch = ref_cast(Champ_Face_PolyMAC_P0, equation().inconnue().valeur());
  const Zone_PolyMAC_P0& zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const IntTab& f_e = zone.face_voisins(), &fcl = ch.fcl(), &e_f = zone.elem_faces();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &pf = equation().milieu().porosite_face(), &ve = zone.volumes(), &vf = zone.volumes_entrelaces(), &fs = zone.face_surfaces();
  const DoubleTab& vf_dir = zone.volumes_entrelaces_dir(), &n_f = zone.face_normales();
  const DoubleTab& pvit = ch.passe(),
                   &alpha = ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().passe(),
                    &press = ref_cast(Pb_Multiphase, equation().probleme()).eq_qdm.pression().passe(),
                     &temp  = ref_cast(Pb_Multiphase, equation().probleme()).eq_energie.inconnue().passe(),
                      &rho   = equation().milieu().masse_volumique().passe(),
                       &mu    = ref_cast(Fluide_base, equation().milieu()).viscosite_dynamique().passe(),
                        &vort  = equation().probleme().get_champ("vorticite").valeurs(),
                         &y_elem = zone.y_elem(),
                          &y_faces = zone.y_faces(),
                           &n_y_elem = zone.normale_paroi_elem(),
                            &n_y_faces = zone.normale_paroi_faces(),
                             &d_bulles = equation().probleme().get_champ("diametre_bulles").valeurs(),
                              &grad_v = equation().probleme().get_champ("gradient_vitesse").valeurs(),
                               * k_turb = (equation().probleme().has_champ("k")) ? &equation().probleme().get_champ("k").passe() : nullptr ;

  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());

  int N = pvit.line_size() , Np = press.line_size(), Nk = (k_turb) ? (*k_turb).dimension(1) : 1, D = dimension,
      nf_tot = zone.nb_faces_tot(), nf = zone.nb_faces(), ne_tot = zone.nb_elem_tot(),
      cR = (rho.dimension_tot(0) == 1), cM = (mu.dimension_tot(0) == 1);

  DoubleTrav a_l(N), p_l(N), T_l(N), rho_l(N), mu_l(N), sigma_l(N,N), k_l(Nk), d_b_l(N), dv(N, N), ddv_c(4), coeff(N, N), //arguments pour coeff
             vr_l(N,D), scal_ur(N), scal_u(N), pvit_l(N, D), vort_l( D==2 ? 1 :D), grad_l(D,D), scal_grad(D); // Requis pour corrections vort et u_l-u-g

  const Portance_interfaciale_base& correlation_pi = ref_cast(Portance_interfaciale_base, correlation_lift_->valeur());

  double vl_norm ;

  /* elements */
  int f;
  for (int e = 0; e < zone.nb_elem_tot(); e++)
    {
      /* arguments de coeff */
      for (int n = 0; n < N; n++)
        {
          a_l(n)   = alpha(e, n);
          p_l(n)   = press(e, n * (Np > 1));
          T_l(n)   = temp(e, n);
          rho_l(n) = rho(!cR * e, n);
          mu_l(n)  = mu(!cM * e, n);
          d_b_l(n) = d_bulles(e,n) ;
          for (int k = 0; k < N; k++)
            if(milc.has_interface(n, k))
              {
                Interface_base& sat = milc.get_interface(n, k);
                sigma_l(n,k) = sat.sigma(temp(e,n), press(e,n * (Np > 1)));
              }

          for (int k = 0; k < N; k++)
            dv(k, n) = ch.v_norm(pvit, pvit, e, -1, k, n, nullptr, nullptr);
        }

      for (int n = 0; n <Nk; n++)   k_l(n)   = (k_turb)   ? (*k_turb)(e,0) : 0;


      correlation_pi.coefficient(a_l, p_l, T_l, rho_l, mu_l, sigma_l, k_l, d_b_l, dv, e, coeff);

      double fac_e = beta_lift_*pe(e) * ve(e);
      int i = nf_tot + D * e;

      // Experimentation sur la portance
      vl_norm = 0;
      scal_ur = 0 ;
      for (int d = 0 ; d < D ; d++) vl_norm += pvit(i+d, n_l)*pvit(i+d, n_l);
      vl_norm = std::sqrt(vl_norm);
      if (vl_norm > 1.e-6)
        {
          for (int k = 0; k < N; k++)
            for (int d = 0 ; d < D ; d++) scal_ur(k) += pvit(i+d, n_l)/vl_norm * (pvit(i+d, k) -pvit(i+d, n_l));
          for (int k = 0; k < N; k++)
            for (int d = 0 ; d < D ; d++) vr_l(k, d)  = pvit(i+d, n_l)/vl_norm * scal_ur(k) ;
        }
      else for (int k=0 ; k<N ; k++)
          for (int d=0 ; d<D ; d++) vr_l(k, d) = pvit(i+d, k) -pvit(i+d, n_l) ;

      if (D==2)
        {
          for (int k = 0; k < N; k++)
            if (k!= n_l) // gas phase
              {
                // Damping of the lift force close to the wall;
                if      (y_elem(e) < .5*d_bulles(e,k)) fac_e *= -1 ; // suppresses lift
                else if (y_elem(e) >    d_bulles(e,k)) fac_e *=  0 ; // no effect
                else                                   fac_e *= (3*std::pow(2*y_elem(e)/d_bulles(e,k)-1, 2) - 2*std::pow(2*y_elem(e)/d_bulles(e,k)-1, 3)) - 1; // partial damping

                secmem(i, n_l) += fac_e * coeff(n_l, k) * vr_l(k, 1) * vort(e, 0) ;
                secmem(i,  k ) -= fac_e * coeff(n_l, k) * vr_l(k, 1) * vort(e, 0) ;
                secmem(i+1,n_l)-= fac_e * coeff(n_l, k) * vr_l(k, 0) * vort(e, 0) ;
                secmem(i+1, k )+= fac_e * coeff(n_l, k) * vr_l(k, 0) * vort(e, 0) ;
              } // 100% explicit
          for (int b = 0; b < e_f.dimension(1) && (f = e_f(e, b)) >= 0; b++)
            if (f<zone.nb_faces())
              if (fcl(f, 0) < 2)
                for (int k = 0; k < N; k++)
                  if (k!= n_l) // gas phase
                    {
                      int c = (e == f_e(f, 0)) ? 0 : 1 ;
                      double fac_f = beta_lift_*pf(f) * vf_dir(f, c);  // Coherence with portance_interfaciale that calculates the correlation at the element

                      if   (y_elem(e) < .5*d_bulles(e,k)) fac_f *= -1 ; // suppresses lift
                      else if (y_elem(e) > d_bulles(e,k)) fac_f *=  0 ; // no effect
                      else                                fac_f *= (3*std::pow(2*y_elem(e)/d_bulles(e,k)-1, 2) - 2*std::pow(2*y_elem(e)/d_bulles(e,k)-1, 3)) - 1; // partial damping

                      secmem(f, n_l) += fac_f * n_f(f, 0)/fs(f) * coeff(n_l, k) * vr_l(k, 1) * vort(e, 0) ;
                      secmem(f,  k ) -= fac_f * n_f(f, 0)/fs(f) * coeff(n_l, k) * vr_l(k, 1) * vort(e, 0) ;
                      secmem(f, n_l) -= fac_f * n_f(f, 1)/fs(f) * coeff(n_l, k) * vr_l(k, 0) * vort(e, 0) ;
                      secmem(f,  k ) += fac_f * n_f(f, 1)/fs(f) * coeff(n_l, k) * vr_l(k, 0) * vort(e, 0) ;
                    } // 100% explicit

        }

      if (D==3)
        {
          for (int k = 0; k < N; k++)
            if (k!= n_l) // gas phase
              {
                // Damping of the lift force close to the wall;
                if (y_elem(e) < .5*d_bulles(e,k)) fac_e *= -1 ; // suppresses lift
                else if (y_elem(e) >    d_bulles(e,k)) fac_e *=  0 ; // no effect
                else                              fac_e *= (3*std::pow(2*y_elem(e)/d_bulles(e,k)-1, 2) - 2*std::pow(2*y_elem(e)/d_bulles(e,k)-1, 3)) - 1; // partial damping

                secmem(i, n_l) += fac_e * coeff(n_l, k) * (vr_l(k, 1) * vort(e, n_l*D+ 2) - vr_l(k, 2) * vort(e, n_l*D+ 1)) ;
                secmem(i,  k ) -= fac_e * coeff(n_l, k) * (vr_l(k, 1) * vort(e, n_l*D+ 2) - vr_l(k, 2) * vort(e, n_l*D+ 1)) ;
                secmem(i+1,n_l)+= fac_e * coeff(n_l, k) * (vr_l(k, 2) * vort(e, n_l*D+ 0) - vr_l(k, 0) * vort(e, n_l*D+ 2)) ;
                secmem(i+1, k )-= fac_e * coeff(n_l, k) * (vr_l(k, 2) * vort(e, n_l*D+ 0) - vr_l(k, 0) * vort(e, n_l*D+ 2)) ;
                secmem(i+2,n_l)+= fac_e * coeff(n_l, k) * (vr_l(k, 0) * vort(e, n_l*D+ 1) - vr_l(k, 1) * vort(e, n_l*D+ 0)) ;
                secmem(i+2, k )-= fac_e * coeff(n_l, k) * (vr_l(k, 0) * vort(e, n_l*D+ 1) - vr_l(k, 1) * vort(e, n_l*D+ 0)) ;
              } // 100% explicit

          /*          for (int b = 0; b < e_f.dimension(1) && (f = e_f(e, b)) >= 0; b++)
                      if (f<zone.nb_faces())
                        if (fcl(f, 0) < 2)
                          for (int k = 0; k < N; k++)
                            if (k!= n_l) // gas phase
                              {
                                int c = (e == f_e(f, 0)) ? 0 : 1 ;
                                double fac_f = beta_*pf(f) * vf_dir(f, c);

                                if (y_elem(e) < .5*d_bulles(e,k)) fac_f *= -1 ; // suppresses lift
                                else if (y_elem(e) >    d_bulles(e,k)) fac_f *=  0 ; // no effect
                                else                              fac_f *= (3*std::pow(2*y_elem(e)/d_bulles(e,k)-1, 2) - 2*std::pow(2*y_elem(e)/d_bulles(e,k)-1, 3)) - 1; // partial damping

                                secmem(f, n_l) += fac_f * n_f(f, 0)/fs(f) * coeff(n_l, k) * (vr_l(k, 1) * vort(e, n_l*D+ 2) - vr_l(k, 2) * vort(e, n_l*D+ 1)) ;
                                secmem(f,  k ) -= fac_f * n_f(f, 0)/fs(f) * coeff(n_l, k) * (vr_l(k, 1) * vort(e, n_l*D+ 2) - vr_l(k, 2) * vort(e, n_l*D+ 1)) ;
                                secmem(f, n_l) += fac_f * n_f(f, 1)/fs(f) * coeff(n_l, k) * (vr_l(k, 2) * vort(e, n_l*D+ 0) - vr_l(k, 0) * vort(e, n_l*D+ 2)) ;
                                secmem(f,  k ) -= fac_f * n_f(f, 1)/fs(f) * coeff(n_l, k) * (vr_l(k, 2) * vort(e, n_l*D+ 0) - vr_l(k, 0) * vort(e, n_l*D+ 2)) ;
                                secmem(f, n_l) += fac_f * n_f(f, 2)/fs(f) * coeff(n_l, k) * (vr_l(k, 0) * vort(e, n_l*D+ 1) - vr_l(k, 1) * vort(e, n_l*D+ 0)) ;
                                secmem(f,  k ) -= fac_f * n_f(f, 2)/fs(f) * coeff(n_l, k) * (vr_l(k, 0) * vort(e, n_l*D+ 1) - vr_l(k, 1) * vort(e, n_l*D+ 0)) ;
                              } // 100% explicit */

        }

    }

  int c, e, n, k, d, d2;
  double fac_f ;

  if (D==3)
    for (f = 0 ; f<zone.nb_faces() ; f++)
      if (fcl(f, 0) < 2)
        {
          a_l = 0;
          p_l = 0;
          T_l = 0;
          rho_l = 0;
          mu_l = 0;
          d_b_l = 0 ;
          sigma_l = 0;
          dv = 0 ;
          k_l = 0;
          for (c = 0; c < 2 && (e = f_e(f, c)) >= 0; c++)
            {
              for (n = 0; n < N; n++)
                {
                  a_l(n)   += vf_dir(f, c)/vf(f) * alpha(e, n);
                  p_l(n)   += vf_dir(f, c)/vf(f) * press(e, n * (Np > 1));
                  T_l(n)   += vf_dir(f, c)/vf(f) * temp(e, n);
                  rho_l(n) += vf_dir(f, c)/vf(f) * rho(!cR * e, n);
                  mu_l(n)  += vf_dir(f, c)/vf(f) * mu(!cM * e, n);
                  d_b_l(n) += vf_dir(f, c)/vf(f) * d_bulles(e,n);
                  for (k = 0; k < N; k++)
                    if(milc.has_interface(n, k))
                      {
                        Interface_base& sat = milc.get_interface(n, k);
                        sigma_l(n,k) += vf_dir(f, c)/vf(f) * sat.sigma(temp(e,n), press(e,n * (Np > 1)));
                      }
                  for (k = 0; k < N; k++)
                    dv(k, n) += vf_dir(f, c)/vf(f) * ch.v_norm(pvit, pvit, e, f, k, n, nullptr, nullptr);
                }
              for (n = 0; n < Nk; n++)  k_l(n)   += (k_turb)   ? vf_dir(f, c)/vf(f) * (*k_turb)(e,0) : 0;
            }

          correlation_pi.coefficient(a_l, p_l, T_l, rho_l, mu_l, sigma_l, k_l, d_b_l, dv, e, coeff);

          grad_l = 0; // we fill grad_l so that grad_l(d, d2) = du_d/dx_d2 by averaging between both elements
          for (d = 0 ; d<D ; d++)
            for (d2 = 0 ; d2<D ; d2++)
              for (c=0 ; c<2  && (e = f_e(f, c)) >= 0; c++)
                grad_l(d, d2) += vf_dir(f, c)/vf(f)*grad_v(nf_tot + D*e + d2 , n_l * D + d) ;
          //We replace the n_l components by the one calculated without interpolation to elements
          scal_grad = 0 ; // scal_grad(d) = grad(u_d).n_f
          for (d = 0 ; d<D ; d++)
            for (d2 = 0 ; d2<D ; d2++)
              scal_grad(d) += grad_l(d, d2)*n_f(f, d2)/fs(f);
          for (d = 0 ; d<D ; d++)
            for (d2 = 0 ; d2<D ; d2++)
              grad_l(d, d2) += (grad_v(f ,n_l*D+d) - scal_grad(d)) * n_f(f, d2)/fs(f);
          // We calculate the local vorticity using this local gradient
          vort_l(0) = grad_l(2, 1) - grad_l(1, 2); // dUz/dy - dUy/dz
          vort_l(1) = grad_l(0, 2) - grad_l(2, 0); // dUx/dz - dUz/dx
          vort_l(2) = grad_l(1, 0) - grad_l(0, 1); // dUy/dx - dUx/dy

          // We also need to calculate relative velocity at the face
          pvit_l = 0 ;
          for (d = 0 ; d<D ; d++)
            for (k = 0 ; k<N ; k++)
              for (c=0 ; c<2 && (e = f_e(f, c)) >= 0; c++)
                pvit_l(k, d) += vf_dir(f, c)/vf(f)*pvit(zone.nb_faces_tot()+D*e+d, k) ;
          scal_u = 0;
          for (k = 0 ; k<N ; k++)
            for (d = 0 ; d<D ; d++)
              scal_u(k) += pvit_l(k, d)*n_f(f, d)/fs(f);
          for (k = 0 ; k<N ; k++)
            for (d = 0 ; d<D ; d++)
              pvit_l(k, d) += (pvit(f, k) - scal_u(k)) * n_f(f, d)/fs(f) ; // Corect velocity at the face
          vl_norm = 0;
          scal_ur = 0;
          for (d = 0 ; d < D ; d++) vl_norm += pvit_l(n_l, d)*pvit_l(n_l, d);
          vl_norm = std::sqrt(vl_norm);
          if (vl_norm > 1.e-6)
            {
              for (k = 0; k < N; k++)
                for (d = 0 ; d < D ; d++) scal_ur(k) += pvit_l(n_l, d)/vl_norm * (pvit_l(k, d) -pvit_l(n_l, d));
              for (k = 0; k < N; k++)
                for (d = 0 ; d < D ; d++) vr_l(k, d)  = pvit_l(n_l, d)/vl_norm * scal_ur(k) ;
            }
          else for (k=0 ; k<N ; k++)
              for (d=0 ; d<D ; d++) vr_l(k, d) = pvit_l(k, d)-pvit_l(n_l, d) ;


          // Use local vairables for the calculation of secmem

          for (k = 0; k < N; k++)
            if (k!= n_l) // gas phase
              {
                fac_f = beta_lift_*pf(f) * vf(f);

                if   (y_faces(f) < .5*d_b_l(k)) fac_f *= -1 ; // suppresses lift
                else if (y_faces(f) > d_b_l(k)) fac_f *=  0 ; // no effect
                else                            fac_f *= (3*std::pow(2*y_faces(f)/d_b_l(k)-1, 2) - 2*std::pow(2*y_faces(f)/d_b_l(k)-1, 3)) - 1; // partial damping

                secmem(f, n_l) += fac_f * n_f(f, 0)/fs(f) * coeff(n_l, k) * (vr_l(k, 1) * vort_l(2) - vr_l(k, 2) * vort_l(1)) ;
                secmem(f,  k ) -= fac_f * n_f(f, 0)/fs(f) * coeff(n_l, k) * (vr_l(k, 1) * vort_l(2) - vr_l(k, 2) * vort_l(1)) ;
                secmem(f, n_l) += fac_f * n_f(f, 1)/fs(f) * coeff(n_l, k) * (vr_l(k, 2) * vort_l(0) - vr_l(k, 0) * vort_l(2)) ;
                secmem(f,  k ) -= fac_f * n_f(f, 1)/fs(f) * coeff(n_l, k) * (vr_l(k, 2) * vort_l(0) - vr_l(k, 0) * vort_l(2)) ;
                secmem(f, n_l) += fac_f * n_f(f, 2)/fs(f) * coeff(n_l, k) * (vr_l(k, 0) * vort_l(1) - vr_l(k, 1) * vort_l(0)) ;
                secmem(f,  k ) -= fac_f * n_f(f, 2)/fs(f) * coeff(n_l, k) * (vr_l(k, 0) * vort_l(1) - vr_l(k, 1) * vort_l(0)) ;
              } // 100% explicit
        }


}
