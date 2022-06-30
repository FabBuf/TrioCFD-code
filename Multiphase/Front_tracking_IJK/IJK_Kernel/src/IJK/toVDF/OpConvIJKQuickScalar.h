//TRUST_NO_INDENT
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : OpConvIJKQuickScalar.h
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////
//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file OpConvIJKQuickScalar.h.P
//
#ifndef OpConvIJKQuickScalar_include
#define OpConvIJKQuickScalar_include

#include <IJK_Splitting.h>
#include <Operateur_IJK_base.h>
#include <Operateur_IJK_data_channel.h>


class OpConvIJKQuickScalar_double : public Operateur_IJK_elem_base_double
{
 public:
  OpConvIJKQuickScalar_double() { stored_curv_fram_layer_z_ = -1000; }
  void initialize(const IJK_Splitting & splitting);
  void calculer(const IJK_Field_double & field,
		const IJK_Field_double & vx, const IJK_Field_double & vy, const IJK_Field_double & vz,
		IJK_Field_double & result);
  void ajouter(const IJK_Field_double & field,
	       const IJK_Field_double & vx, const IJK_Field_double & vy, const IJK_Field_double & vz,
	       IJK_Field_double & result);
 protected:
  
  void compute_curv_fram_x(int k_layer);
  void compute_flux_x(IJK_Field_local_double & resu, const int k_layer) override;
  void compute_curv_fram_y(int k_layer);
  void compute_flux_y(IJK_Field_local_double & resu, const int k_layer) override;
  void compute_curv_fram_z(int k_layer);
  void compute_flux_z(IJK_Field_local_double & resu, const int k_layer) override;
  Operateur_IJK_data_channel channel_data_;

  // Pointers to input data (set by calculer, used by compute_flux_...)
  const IJK_Field_local_double *input_field_;
  const IJK_Field_local_double *input_velocity_x_;
  const IJK_Field_local_double *input_velocity_y_;
  const IJK_Field_local_double *input_velocity_z_;
  bool perio_k_ ;

  // Temporary array to store curvature and fram coefficients
  // for the current computed flux.
  // layer k=0 and k=1 are used for "curv", k=2 and k=3 are used for "fram".
  // layer k=0 and k=2 store the previous values computed in direction "k" (which is used 2 times)
  IJK_Field_local_double tmp_curv_fram_;
  int stored_curv_fram_layer_z_; // which (local) layer is currently stored in layer 0 of the tmp array ?
  
};

class OpConvCentre2IJKScalar_double : public Operateur_IJK_elem_base_double
{
 public:
  OpConvCentre2IJKScalar_double() { stored_curv_fram_layer_z_ = -1000; }
  void initialize(const IJK_Splitting & splitting);
  void calculer(const IJK_Field_double & field,
		const IJK_Field_double & vx, const IJK_Field_double & vy, const IJK_Field_double & vz,
		IJK_Field_double & result);
  void ajouter(const IJK_Field_double & field,
	       const IJK_Field_double & vx, const IJK_Field_double & vy, const IJK_Field_double & vz,
	       IJK_Field_double & result);
 protected:
  
  void compute_curv_fram_x(int k_layer);
  void compute_flux_x(IJK_Field_local_double & resu, const int k_layer) override;
  void compute_curv_fram_y(int k_layer);
  void compute_flux_y(IJK_Field_local_double & resu, const int k_layer) override;
  void compute_curv_fram_z(int k_layer);
  void compute_flux_z(IJK_Field_local_double & resu, const int k_layer) override;
  Operateur_IJK_data_channel channel_data_;

  // Pointers to input data (set by calculer, used by compute_flux_...)
  const IJK_Field_local_double *input_field_;
  const IJK_Field_local_double *input_velocity_x_;
  const IJK_Field_local_double *input_velocity_y_;
  const IJK_Field_local_double *input_velocity_z_;
  bool perio_k_ ;

  // Temporary array to store curvature and fram coefficients
  // for the current computed flux.
  // layer k=0 and k=1 are used for "curv", k=2 and k=3 are used for "fram".
  // layer k=0 and k=2 store the previous values computed in direction "k" (which is used 2 times)
  IJK_Field_local_double tmp_curv_fram_;
  int stored_curv_fram_layer_z_; // which (local) layer is currently stored in layer 0 of the tmp array ?
  
};


#endif
