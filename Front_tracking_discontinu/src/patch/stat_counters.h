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
// File:        stat_counters.h
// Directory:   $TRUST_ROOT/src/Kernel/Utilitaires
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Stat_counters_included
#define Stat_counters_included
#include <Statistiques.h>

extern void print_statistics_analyse(const char * message, int mode_append);
extern void declare_stat_counters();

extern Stat_Counter_Id mpi_sendrecv_counter_;
extern Stat_Counter_Id mpi_send_counter_;
extern Stat_Counter_Id mpi_recv_counter_;
extern Stat_Counter_Id mpi_bcast_counter_;
extern Stat_Counter_Id mpi_alltoall_counter_;
extern Stat_Counter_Id mpi_partialsum_counter_;
extern Stat_Counter_Id mpi_sumdouble_counter_;
extern Stat_Counter_Id mpi_mindouble_counter_;
extern Stat_Counter_Id mpi_maxdouble_counter_;
extern Stat_Counter_Id mpi_sumint_counter_;
extern Stat_Counter_Id mpi_minint_counter_;
extern Stat_Counter_Id mpi_maxint_counter_;
extern Stat_Counter_Id mpi_barrier_counter_;

extern Stat_Counter_Id mpi_sendrecv_io_counter_;
extern Stat_Counter_Id echange_vect_counter_;
extern Stat_Counter_Id solv_sys_counter_;
extern Stat_Counter_Id solv_sys_petsc_counter_;
extern Stat_Counter_Id diffusion_implicite_counter_;
extern Stat_Counter_Id dt_counter_;
extern Stat_Counter_Id nut_counter_;
extern Stat_Counter_Id convection_counter_;
extern Stat_Counter_Id diffusion_counter_;
extern Stat_Counter_Id decay_counter_;
extern Stat_Counter_Id source_counter_;
extern Stat_Counter_Id source_acc_counter_;
extern Stat_Counter_Id source_interf_counter_;
extern Stat_Counter_Id divergence_counter_;
extern Stat_Counter_Id gradient_counter_;
extern Stat_Counter_Id postraitement_counter_;
extern Stat_Counter_Id divers_counter_;
extern Stat_Counter_Id sauvegarde_counter_;
extern Stat_Counter_Id temporary_counter_;
extern Stat_Counter_Id assemblage_sys_counter_;
extern Stat_Counter_Id update_vars_counter_;
extern Stat_Counter_Id update_fields_counter_;
extern Stat_Counter_Id mettre_a_jour_counter_;
extern Stat_Counter_Id timestep_counter_;
extern Stat_Counter_Id interprete_scatter_counter_;
extern Stat_Counter_Id temps_total_execution_counter_;
extern Stat_Counter_Id initialisation_calcul_counter_;
extern Stat_Counter_Id m1;
extern Stat_Counter_Id m2;
extern Stat_Counter_Id m3;

extern Stat_Counter_Id probleme_fluide_;
extern Stat_Counter_Id probleme_combustible_;

extern Stat_Counter_Id euler_rk3_counter_ ;
extern Stat_Counter_Id calcul_dv_counter_;
extern Stat_Counter_Id projection_counter_ ;
extern Stat_Counter_Id multigrille_counter_ ;
extern Stat_Counter_Id jacobi_residu_counter_;
extern Stat_Counter_Id coarsen_counter_;
extern Stat_Counter_Id interpolate_counter_;

extern Stat_Counter_Id calculer_rho_mu_indicatrice_counter_;
extern Stat_Counter_Id calculer_indicatrice_counter_;
extern Stat_Counter_Id parcours_maillage_counter_;
extern Stat_Counter_Id search_connex_components_counter_;
extern Stat_Counter_Id calculs_drapeaux_counter_;
extern Stat_Counter_Id boucle_indicatrice_counter_;

extern Stat_Counter_Id deplacement_interf_counter_;
extern Stat_Counter_Id redistribute_dplct_interf_counter_;

extern Stat_Counter_Id barycentre_lissage_sys_counter_;
extern Stat_Counter_Id remaillage_loc_interf_counter_;
extern Stat_Counter_Id barycentre_lissage_apres_counter_;
extern Stat_Counter_Id sup_div_aretes_counter_;

extern Stat_Counter_Id updtstat_counter_;
extern Stat_Counter_Id bilanQdM_counter_;
extern Stat_Counter_Id redistribute_tot_counter_;

#endif
