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
// File      : Corrige_flux_FT_temperature_subresolution.cpp
// Directory : $TRIOCFD_ROOT/src/Multiphase/Front_tracking_IJK/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <Corrige_flux_FT_temperature_subresolution.h>
#include <IJK_FT.h>

Implemente_instanciable_sans_constructeur( Corrige_flux_FT_temperature_subresolution, "Corrige_flux_FT_temperature_subresolution", Corrige_flux_FT_base ) ;

Corrige_flux_FT_temperature_subresolution::Corrige_flux_FT_temperature_subresolution()
{
  thermal_subproblems_ = nullptr;
  has_checked_consistency_ = false;
}

void Corrige_flux_FT_temperature_subresolution::associate_thermal_problems(const IJK_One_Dimensional_Subproblems& thermal_subproblems)
{
  thermal_subproblems_ = &thermal_subproblems;
}

Sortie& Corrige_flux_FT_temperature_subresolution::printOn( Sortie& os ) const
{
  Corrige_flux_FT_base::printOn( os );
  return os;
}

Entree& Corrige_flux_FT_temperature_subresolution::readOn( Entree& is )
{
  Corrige_flux_FT_base::readOn( is );
  return is;
}

void Corrige_flux_FT_temperature_subresolution::initialize_with_subproblems(const IJK_Splitting& splitting,
                                                                            const IJK_Field_double& field,
                                                                            const IJK_Interfaces& interfaces,
                                                                            const IJK_FT_double& ijk_ft,
                                                                            Intersection_Interface_ijk_face& intersection_ijk_face,
                                                                            Intersection_Interface_ijk_cell& intersection_ijk_cell,
                                                                            const IJK_One_Dimensional_Subproblems& thermal_subproblems)
{
  Corrige_flux_FT_base::initialize_with_subproblems(splitting, field, interfaces, ijk_ft, intersection_ijk_face, intersection_ijk_cell, thermal_subproblems);
  associate_thermal_problems(thermal_subproblems);
}

void Corrige_flux_FT_temperature_subresolution::update_intersections()
{

  // On commence par calculer les temperatures aux faces mouillées
  intersection_ijk_cell_->update_interpolations_cell_centres_on_interface();
  /*
   * TODO update with face cell centres positions
   */
  if (!convection_negligible_ || !diffusion_negligible_)
    intersection_ijk_cell_->update_interpolations_cell_faces_on_interface();

  Cerr << "The intersections have been updated" << finl;
}

void Corrige_flux_FT_temperature_subresolution::update()
{
//  ArrOfDouble temp_vap, temp_liqu;
//  temp_vap.set_smart_resize(1);
//  temp_liqu.set_smart_resize(1);
//  // On commence par calculer les temperatures aux faces mouillées
//  intersection_ijk_cell_->update_interpolations_cell_centres_on_interface();
//  /*
//   * TODO update with face cell centres positions
//   */
//  if (!convection_negligible_ && !diffusion_negligible_)
//    intersection_ijk_cell_->update_interpolations_cell_faces_on_interface();

  associate_indices_and_check_subproblems_consistency();
}

void Corrige_flux_FT_temperature_subresolution::associate_indices_and_check_subproblems_consistency()
{
  if (!has_checked_consistency_)
    {
      const int nb_diph = intersection_ijk_cell_->get_nb_diph();
      const int nb_subproblems = thermal_subproblems_->size();
      has_checked_consistency_ = (nb_diph==nb_subproblems);
      assert(has_checked_consistency_);
      ijk_intersections_subproblems_indices_.reset();
      ijk_intersections_subproblems_indices_.resize(nb_subproblems);
      int index_i_problem = 0;
      int index_j_problem = 0;
      int index_k_problem = 0;
      int problem_index;
      for (problem_index=0; problem_index<nb_subproblems; problem_index++)
        {
          thermal_subproblems_->get_subproblem_ijk_indices(index_i_problem, index_j_problem, index_k_problem, problem_index);
          const int ijk_intersections_index = (*intersection_ijk_cell_)(index_i_problem, index_j_problem, index_k_problem);
          if (ijk_intersections_index > -1)
            {
              ijk_intersections_subproblems_indices_[problem_index] = ijk_intersections_index;
              has_checked_consistency_ = has_checked_consistency_ && true;
            }
          else
            {
              Cerr << "Inconsistency between intersection_ijk_cell and the thermal_subproblem (index " << problem_index << ")" << finl;
              has_checked_consistency_ = has_checked_consistency_ && false;
            }
        }
      assert(has_checked_consistency_);
    }
  else
    Cerr << "Inconsistency has already be checked" << finl;
}

void Corrige_flux_FT_temperature_subresolution::clean()
{
  has_checked_consistency_=false;
  intersection_ijk_cell_->set_pas_a_jour();
}


void Corrige_flux_FT_temperature_subresolution::compute_temperature_cell_centre(IJK_Field_double& temperature) const
{
  /*
   * For each subproblem fill the right interfacial_cell
   */
  const DoubleTab dist_interf = intersection_ijk_cell_->dist_interf();
  const double min_temperature = thermal_subproblems_->get_min_temperature_domain_ends();
  const double max_temperature = thermal_subproblems_->get_max_temperature_domain_ends();
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      const double dist = dist_interf(intersection_ijk_cell_index, 0);

      double temperature_ghost = 0.;
      temperature_ghost = thermal_subproblems_->get_temperature_profile_at_point(i, dist);

      const int ijk_indices_i = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 0);
      const int ijk_indices_j = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 1);
      const int ijk_indices_k = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 2);

      const IJK_Field_double& indicator = ref_ijk_ft_->itfce().I();
      const double indic = indicator(ijk_indices_i, ijk_indices_j, ijk_indices_k);
      if (temperature_ghost < min_temperature && indic > 0.5)
        Cerr << "Ghost temperature: " << temperature_ghost << " is lower than the minimum temperature:" << min_temperature << finl;
      if (temperature_ghost > max_temperature && indic > 0.5)
        Cerr << "Ghost temperature: " << temperature_ghost << " is higher than the maximum temperature:" << max_temperature << finl;

      temperature(ijk_indices_i, ijk_indices_j, ijk_indices_k) = temperature_ghost;
    }
}

void Corrige_flux_FT_temperature_subresolution::set_zero_temperature_increment(IJK_Field_double& d_temperature) const
{
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      const int ijk_indices_i = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 0);
      const int ijk_indices_j = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 1);
      const int ijk_indices_k = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 2);
      d_temperature(ijk_indices_i, ijk_indices_j, ijk_indices_k) = 0.;
    }
}

void Corrige_flux_FT_temperature_subresolution::compute_temperature_face_centre()
{
  int faces_dir[6] = FACES_DIR;
  convective_fluxes_.set_smart_resize(1);
  convective_fluxes_.reset();
  dist_.set_smart_resize(1);
  dist_.reset();
  const DoubleTab dist_interf = intersection_ijk_cell_->dist_pure_faces_interf();
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      double surf_face = 0.;
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      for (int l=0; l<6; l++)
        {
          const int is_neighbour_pure_liquid = intersection_ijk_cell_->get_ijk_pure_face_neighbours(intersection_ijk_cell_index, l);
          if (is_neighbour_pure_liquid)
            {
              for (int c = 0; c < 3; c++)
                if (c!= faces_dir[l])
                  surf_face *= splitting_->get_grid_geometry().get_constant_delta(c);
              const double dist = dist_interf(intersection_ijk_cell_index, l);
              const double temperature_face = thermal_subproblems_->get_temperature_profile_at_point(i, dist);
              const double flux_face = temperature_face * surf_face;
              convective_fluxes_.append_array(flux_face);
              dist_.append_array(dist);
            }
        }
    }
  /*
   * Useless if a treat a sub-problem per mixed cells
   * May be useful if a treat one subproblem per interface portion
   */
  // check_pure_fluxes_duplicates(convective_fluxes_, convective_fluxes_unique_, pure_face_unique_, 0);
}

void Corrige_flux_FT_temperature_subresolution::compute_thermal_fluxes_face_centre()
{
  int faces_dir[6] = FACES_DIR;
  diffusive_fluxes_.set_smart_resize(1);
  diffusive_fluxes_.reset();
  const DoubleTab dist_interf = intersection_ijk_cell_->dist_pure_faces_interf();
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      double surf_face = 0.;
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      for (int l=0; l<6; l++)
        {
          const int is_neighbour_pure_liquid = intersection_ijk_cell_->get_ijk_pure_face_neighbours(intersection_ijk_cell_index, l);
          if (is_neighbour_pure_liquid)
            {
              for (int c = 0; c < 3; c++)
                if (c!= faces_dir[l])
                  surf_face *= splitting_->get_grid_geometry().get_constant_delta(c);
              const double dist = dist_interf(intersection_ijk_cell_index, l);
              const double temperature_gradient_face = thermal_subproblems_->get_temperature_gradient_profile_at_point(i, dist, faces_dir[l]);
              const double flux_face = temperature_gradient_face * surf_face;
              diffusive_fluxes_.append_array(flux_face);
            }
        }
    }
  /*
   * Useless if a treat a sub-problem per mixed cells
   * May be useful if a treat one subproblem per interface portion
   */
  // check_pure_fluxes_duplicates(convective_fluxes_, convective_fluxes_unique_, pure_face_unique_, 1);
}

void Corrige_flux_FT_temperature_subresolution::compute_ijk_pure_faces_indices()
{
  /*
   * Be careful, the ijk_intersection class is not sorting the faces the same way
   */
//	FixedVector<DoubleVect, 3>& ijk_faces_to_correct = intersection_ijk_cell_->get_set_ijk_pure_face_to_correct();
  int faces_dir[6] = FACES_DIR;
  int neighbours_faces_i[6] = NEIGHBOURS_FACES_I;
  int neighbours_faces_j[6] = NEIGHBOURS_FACES_J;
  int neighbours_faces_k[6] = NEIGHBOURS_FACES_K;
  int nb_faces_to_correct = 0;
  for (int i=0; i<ijk_intersections_subproblems_indices_.size_array(); i++)
    {
      const int intersection_ijk_cell_index = ijk_intersections_subproblems_indices_[i];
      const int ijk_indices_i = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 0);
      const int ijk_indices_j = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 1);
      const int ijk_indices_k = (*intersection_ijk_cell_)(intersection_ijk_cell_index, 2);
      for (int l=0; l<6; l++)
        {
          const int is_neighbour_pure_liquid = intersection_ijk_cell_->get_ijk_pure_face_neighbours(intersection_ijk_cell_index, l);
          if (is_neighbour_pure_liquid)
            {
              const int ii_f = neighbours_faces_i[l];
              const int jj_f = neighbours_faces_j[l];
              const int kk_f = neighbours_faces_k[l];
              IntVect& i_pure_face_to_correct = ijk_faces_to_correct_[0];
              IntVect& j_pure_face_to_correct = ijk_faces_to_correct_[1];
              IntVect& k_pure_face_to_correct = ijk_faces_to_correct_[2];
              IntVect& dir_pure_face_to_correct = ijk_faces_to_correct_[3];
              i_pure_face_to_correct[nb_faces_to_correct] = (ijk_indices_i + ii_f);
              j_pure_face_to_correct[nb_faces_to_correct] = (ijk_indices_j + jj_f);
              k_pure_face_to_correct[nb_faces_to_correct] = (ijk_indices_k + kk_f);
              dir_pure_face_to_correct[nb_faces_to_correct] = faces_dir[l];
              nb_faces_to_correct++;
            }
        }
    }
  const int convective_fluxes_size = convective_fluxes_.size();
  const int diffusive_fluxes_size = diffusive_fluxes_.size();
  if (convective_fluxes_size > 0)
    assert(convective_fluxes_size ==  ijk_faces_to_correct_[0].size());
  if (diffusive_fluxes_size > 0)
    assert(diffusive_fluxes_size == ijk_faces_to_correct_[0].size());
}

void Corrige_flux_FT_temperature_subresolution::sort_ijk_intersections_subproblems_indices_by_k_layers()
{
  const int nb_k_layer = ref_ijk_ft_->itfce().I().nk();
  index_face_i_sorted_.resize(nb_k_layer);
  index_face_j_sorted_.resize(nb_k_layer);
  convective_flux_x_sorted_.resize(nb_k_layer);
  convective_flux_y_sorted_.resize(nb_k_layer);
  convective_flux_z_sorted_.resize(nb_k_layer);
  diffusive_flux_x_sorted_.resize(nb_k_layer);
  diffusive_flux_y_sorted_.resize(nb_k_layer);
  diffusive_flux_z_sorted_.resize(nb_k_layer);
  FixedVector<std::vector<DoubleVect>*,3> convective_fluxes;
  convective_fluxes[0] = &convective_flux_x_sorted_;
  convective_fluxes[1] = &convective_flux_y_sorted_;
  convective_fluxes[2] = &convective_flux_z_sorted_;
  FixedVector<std::vector<DoubleVect>*,3> diffusive_fluxes;
  diffusive_fluxes[0] = &diffusive_flux_x_sorted_;
  diffusive_fluxes[1] = &diffusive_flux_y_sorted_;
  diffusive_fluxes[2] = &diffusive_flux_z_sorted_;

  IntVect& i_pure_face_to_correct = ijk_faces_to_correct_[0];
  IntVect& j_pure_face_to_correct = ijk_faces_to_correct_[1];
  IntVect& k_pure_face_to_correct = ijk_faces_to_correct_[2];
  IntVect& dir_pure_face_to_correct = ijk_faces_to_correct_[3];
  const int nb_fluxes = ijk_faces_to_correct_[0].size();
  for (int i_flux=0; i_flux < nb_fluxes; i_flux++)
    {
      const int k = k_pure_face_to_correct[i_flux];
      const int dir = dir_pure_face_to_correct[i_flux];
      std::vector<DoubleVect>* flux = convective_fluxes[dir];
      (*flux)[k].append_array(convective_fluxes_[i_flux]);
      index_face_i_sorted_[k].append_array(i_pure_face_to_correct[i_flux]);
      index_face_j_sorted_[k].append_array(j_pure_face_to_correct[i_flux]);
    }
}


void Corrige_flux_FT_temperature_subresolution::check_pure_fluxes_duplicates(const DoubleVect& fluxes,
                                                                             DoubleVect& fluxes_unique,
                                                                             IntVect& pure_face_unique,
                                                                             const int known_unique)
{
  // FixedVector<DoubleVect, 3>& pure_face_to_correct = intersection_ijk_cell_->get_set_ijk_pure_face_to_correct();
  FixedVector<IntVect, 4>& pure_face_to_correct = ijk_faces_to_correct_;
  const int nb_fluxes = fluxes.size();
  fluxes_unique.set_smart_resize(1);
  fluxes_unique.reset();
  int i;
  if (known_unique)
    {
      const int size_face_unique = pure_face_unique.size();
      for (i=0; i<size_face_unique; i++)
        {
          fluxes_unique.append_array(fluxes(pure_face_unique(i)));
        }
    }
  else
    {
      DoubleVect shared_face;
      shared_face.set_smart_resize(1);
      shared_face.reset();
      for (i=0; i<nb_fluxes; i++)
        {
          const int i_f = pure_face_to_correct[0](i);
          const int j_f = pure_face_to_correct[1](i);
          const int k_f = pure_face_to_correct[2](i);
          for (int j=i; j<nb_fluxes; j++)
            if (i != j)
              {
                const int i_ff = pure_face_to_correct[0](j);
                const int j_ff = pure_face_to_correct[1](j);
                const int k_ff = pure_face_to_correct[2](j);
                if ((i_f==i_ff) && (j_f==j_ff) && (k_f==k_ff))
                  {
                    if (shared_face.size() == 0)
                      shared_face.append_array(i);
                    shared_face.append_array(j);
                  }
              }
          /*
           * Take the closest portion
           */
          const int nb_duplicates = shared_face.size();
          int min_dist_index = 0;
          double min_dist = dist_[i];
          for (int j=1; j<nb_duplicates; j++)
            if (min_dist > dist_[j])
              {
                min_dist = dist_[j];
                min_dist_index = j;
              }
          pure_face_unique.append_array(min_dist_index);
          fluxes_unique.append_array(fluxes(min_dist_index));
        }
    }
}
