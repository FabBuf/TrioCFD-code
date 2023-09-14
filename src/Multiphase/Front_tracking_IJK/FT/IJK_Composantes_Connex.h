/****************************************************************************
* Copyright (c) 2023, CEA
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
// File      : IJK_Composantes_Connex.h
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef IJK_Composantes_Connex_included
#define IJK_Composantes_Connex_included

#include <Objet_U.h>
#include <IJK_Field.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class IJK_Composantes_Connex
//
// <Description of class IJK_Composantes_Connex>
//
/////////////////////////////////////////////////////////////////////////////
class IJK_FT_double;
class IJK_Interfaces;
class IJK_Composantes_Connex : public Objet_U
{

  Declare_instanciable( IJK_Composantes_Connex ) ;

public :
  int initialize(const IJK_Splitting& splitting, const IJK_Interfaces& interfaces);
  void associer(const IJK_FT_double& ijk_ft);
  void compute_bounding_box_fill_compo_connex();

  const IJK_Field_double& get_eulerian_compo_connex_ft() const
  {
    return eulerian_compo_connex_ft_;
  }
  const IJK_Field_double& get_eulerian_compo_connex() const
  {
    return eulerian_compo_connex_ns_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_ghost_ft() const
  {
    return eulerian_compo_connex_ghost_ft_;
  }
  const IJK_Field_double& get_eulerian_compo_connex_ghost() const
  {
    return eulerian_compo_connex_ghost_ns_;
  }

  const DoubleTab& get_bounding_box() const
  {
    return bounding_box_;
  }
  const DoubleTab& get_bubbles_barycentre() const
  {
    return bubbles_barycentre_;
  }

protected :
  REF(IJK_FT_double) ref_ijk_ft_;
  const IJK_Interfaces * interfaces_ = nullptr;

  IJK_Field_double eulerian_compo_connex_ft_;
  IJK_Field_double eulerian_compo_connex_ns_;
  IJK_Field_double eulerian_compo_connex_ghost_ft_;
  IJK_Field_double eulerian_compo_connex_ghost_ns_;

  DoubleTab bounding_box_;
  DoubleTab bubbles_barycentre_;

  bool is_updated_ = false;
};

#endif /* IJK_Composantes_Connex_included */
