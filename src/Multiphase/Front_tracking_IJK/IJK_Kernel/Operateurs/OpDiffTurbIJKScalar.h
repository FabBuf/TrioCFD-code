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

#ifndef OpDiffTurbIJKScalar_H
#define OpDiffTurbIJKScalar_H

#include <Operateur_IJK_base.h>

class OpDiffIJKScalarGeneric_double : public Operateur_IJK_elem_base_double
{
public:
  OpDiffIJKScalarGeneric_double();
  void initialize(const IJK_Splitting& splitting);

  inline void compute_flux_x(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::X>(resu,k_layer);
  }
  inline void compute_flux_y(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Y>(resu,k_layer);
  }
  inline void compute_flux_z(IJK_Field_local_double& resu, const int k_layer) override
  {
    compute_flux_<DIRECTION::Z>(resu,k_layer);
  }

protected:
  template <DIRECTION _DIR_>
  inline void compute_flux_(IJK_Field_local_double& resu, const int k_layer);

  const IJK_Field_local_double& get_lambda_vectorial(DIRECTION _DIR_);

  const IJK_Field_local_double& get_structural_model(DIRECTION _DIR_);


  Operateur_IJK_data_channel channel_data_;
  bool perio_k_;

  // Pointers to input data (set by calculer, used by compute_flux_...)
  const IJK_Field_local_double *input_field_;

  const IJK_Field_local_double *lambda_;

  const IJK_Field_local_double *lambda_vector_x_;
  const IJK_Field_local_double *lambda_vector_y_;
  const IJK_Field_local_double *lambda_vector_z_;

  const IJK_Field_local_double *structural_model_x_;
  const IJK_Field_local_double *structural_model_y_;
  const IJK_Field_local_double *structural_model_z_;

  const IJK_Field_local_double *boundary_flux_kmin_;
  const IJK_Field_local_double *boundary_flux_kmax_;

  bool is_anisotropic_;
  bool is_vectorial_;
  bool is_structural_;
};


class OpDiffIJKScalar_double : public OpDiffIJKScalarGeneric_double
{
public:
  OpDiffIJKScalar_double() : OpDiffIJKScalarGeneric_double() {}
  void calculer(const IJK_Field_double& field,
                const IJK_Field_double& molecular_alpha,
                IJK_Field_double& result,
                const IJK_Field_local_double& boundary_flux_kmin,
                const IJK_Field_local_double& boundary_flux_kmax);
  void ajouter(const IJK_Field_double& field,
               const IJK_Field_double& molecular_alpha,
               IJK_Field_double& result,
               const IJK_Field_local_double& boundary_flux_kmin,
               const IJK_Field_local_double& boundary_flux_kmax);
};

class OpDiffAnisotropicIJKScalar_double : public OpDiffIJKScalarGeneric_double
{
public:
  OpDiffAnisotropicIJKScalar_double() : OpDiffIJKScalarGeneric_double() { is_anisotropic_ = true; }
  void calculer(const IJK_Field_double& field,
                const IJK_Field_double& molecular_alpha,
                IJK_Field_double& result,
                const IJK_Field_local_double& boundary_flux_kmin,
                const IJK_Field_local_double& boundary_flux_kmax);
  void ajouter(const IJK_Field_double& field,
               const IJK_Field_double& molecular_alpha,
               IJK_Field_double& result,
               const IJK_Field_local_double& boundary_flux_kmin,
               const IJK_Field_local_double& boundary_flux_kmax);
};

class OpDiffVectorialIJKScalar_double : public OpDiffIJKScalarGeneric_double
{
public:
  OpDiffVectorialIJKScalar_double() : OpDiffIJKScalarGeneric_double() { is_vectorial_ = true; }
  void calculer(const IJK_Field_double& field,
                const IJK_Field_double& molecular_alpha_vector_x,
                const IJK_Field_double& molecular_alpha_vector_y,
                const IJK_Field_double& molecular_alpha_vector_z,
                IJK_Field_double& result,
                const IJK_Field_local_double& boundary_flux_kmin,
                const IJK_Field_local_double& boundary_flux_kmax);
  void ajouter(const IJK_Field_double& field,
               const IJK_Field_double& molecular_alpha_vector_x,
               const IJK_Field_double& molecular_alpha_vector_y,
               const IJK_Field_double& molecular_alpha_vector_z,
               IJK_Field_double& result,
               const IJK_Field_local_double& boundary_flux_kmin,
               const IJK_Field_local_double& boundary_flux_kmax);
};

class OpDiffVectorialAnisotropicIJKScalar_double : public OpDiffIJKScalarGeneric_double
{
public:
  OpDiffVectorialAnisotropicIJKScalar_double() : OpDiffIJKScalarGeneric_double() { is_vectorial_ = true, is_anisotropic_ = true; }
  void calculer(const IJK_Field_double& field,
                const IJK_Field_double& molecular_alpha_vector_x,
                const IJK_Field_double& molecular_alpha_vector_y,
                const IJK_Field_double& molecular_alpha_vector_z,
                IJK_Field_double& result,
                const IJK_Field_local_double& boundary_flux_kmin,
                const IJK_Field_local_double& boundary_flux_kmax);
  void ajouter(const IJK_Field_double& field,
               const IJK_Field_double& molecular_alpha_vector_x,
               const IJK_Field_double& molecular_alpha_vector_y,
               const IJK_Field_double& molecular_alpha_vector_z,
               IJK_Field_double& result,
               const IJK_Field_local_double& boundary_flux_kmin,
               const IJK_Field_local_double& boundary_flux_kmax);
};

class OpDiffIJKScalarStructuralOnly_double : public OpDiffIJKScalarGeneric_double
{
public:
  OpDiffIJKScalarStructuralOnly_double() : OpDiffIJKScalarGeneric_double() { is_structural_ = true; }
  void calculer(const IJK_Field_double& field,
                const IJK_Field_double& structural_model_x,
                const IJK_Field_double& structural_model_y,
                const IJK_Field_double& structural_model_z,
                IJK_Field_double& result,
                const IJK_Field_local_double& boundary_flux_kmin,
                const IJK_Field_local_double& boundary_flux_kmax);
  void ajouter(const IJK_Field_double& field,
               const IJK_Field_double& structural_model_x,
               const IJK_Field_double& structural_model_y,
               const IJK_Field_double& structural_model_z,
               IJK_Field_double& result,
               const IJK_Field_local_double& boundary_flux_kmin,
               const IJK_Field_local_double& boundary_flux_kmax);

};

#include <OpDiffTurbIJKScalar.tpp>


#endif
