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
// File:        Transport_Interfaces_FT_Disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/38
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Transport_Interfaces_FT_Disc_included
#define Transport_Interfaces_FT_Disc_included

#include <Equation_base.h>
#include <Transport_Interfaces_base.h>
#include <Postraitement_base.h>
#include <Champ_Inc.h>
#include <Remaillage_FT.h>
#include <Parcours_interface.h>
#include <Marching_Cubes.h>
#include <Connectivite_frontieres.h>
#include <Topologie_Maillage_FT.h>
#include <Algorithmes_Transport_FT_Disc.h>
#include <Champ_Fonc.h>
#include <Navier_Stokes_FT_Disc.h>
#include <Proprietes_part_vol.h>
#include <TRUSTTabFT_forward.h>
#include <TRUST_Ref.h>
// EB
#include <Modele_Collision_FT.h>
#include <Particule_Solide.h>
//#include <Ref_Joint.h> // EB
#include <Postraitement_Forces_Interfaces_FT.h>
//#include <Ref_Operateur_Diff.h>
#include <Convection_Diffusion_Temperature_FT_Disc.h>
// fin EB
class Probleme_base;
class Milieu_base;
class Navier_Stokes_FT_Disc;
class Convection_Diffusion_Temperature_FT_Disc; // EB
class Loi_horaire;

class Transport_Interfaces_FT_Disc_interne;
template <typename titi> class TRUSTTab;
using FloatTab = TRUSTTab<float>;

class Transport_Interfaces_FT_Disc : public Transport_Interfaces_base
{

  friend Convection_Diffusion_Temperature_FT_Disc; // EB

  Declare_instanciable_sans_constructeur(Transport_Interfaces_FT_Disc);
public:

  Transport_Interfaces_FT_Disc();
  //
  void set_param(Param& titi) override;
  long lire_motcle_non_standard(const Motcle&, Entree&) override;
  // Methodes virtuelles pures de Equation_base
  //
  long            nombre_d_operateurs(void) const override; // Zero, y'en a pas.
  const Operateur& operateur(long i) const override;    // Erreur
  Operateur&        operateur(long i) override;         // Erreur
  const Champ_Inc& inconnue(void) const override;         // C'est l'indicatrice
  Champ_Inc&        inconnue(void) override;
  Champ_Inc&        inconnue_face(void); // EB
  //
  // Methodes surchargees de Equation_base
  //
  void                associer_milieu_base(const Milieu_base& milieu) override;

  void                 associer_equation_ns(const Navier_Stokes_FT_Disc& ns);
  void                 associer_equation_temp(const Convection_Diffusion_Temperature_FT_Disc& temp); // EB
  Milieu_base&        milieu() override;       // Erreur
  const Milieu_base& milieu() const override;  // Erreur
  void    associer_pb_base(const Probleme_base& probleme) override;
  void    discretiser(void) override;
  Entree& lire_cond_init(Entree& is) override;
  long  verif_Cl() const override;
  double  calculer_pas_de_temps(void) const override;
  DoubleTab& derivee_en_temps_inco(DoubleTab& derivee) override;
  void assembler( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) override ;

  void    mettre_a_jour(double temps) override;
  long  sauvegarder(Sortie& ) const override;
  long  reprendre(Entree&) override;
  long impr(Sortie& os) const override;
  virtual long impr_fpi(Sortie& os) const override; // EB
  void update_critere_statio();

  //
  // Nouvelles methodes
  //
  virtual void                            lire_maillage_ft_cao(Entree& is);
  long                          preparer_calcul() override;
  virtual void                            preparer_pas_de_temps();
  const Maillage_FT_Disc&                 maillage_interface() const;
  const Champ_base&               get_update_indicatrice() override;
  virtual const Champ_base&               get_indicatrice_faces();
  virtual const DoubleVect&               get_indicatrice_aretes(); // EB
  virtual const DoubleVect&               get_indicatrice_aretes() const; // EB
  virtual const Champ_base&               get_compute_indicatrice_faces();
  virtual const DoubleTab&                get_compute_indicatrice_aretes_internes(); // EB
  virtual const Parcours_interface&       parcours_interface() const;
  virtual const Marching_Cubes&           marching_cubes() const;
  virtual const Algorithmes_Transport_FT_Disc& algorithmes_transport() const;
  virtual const Connectivite_frontieres& connectivite_frontieres() const;
  Remaillage_FT&                          remaillage_interface();
  const Remaillage_FT&                    remaillage_interface() const;
  const Topologie_Maillage_FT&            topologie_interface() const;
  virtual double calculer_integrale_indicatrice(const DoubleVect& indicatrice, double& v_ph0) const;

  // debut EB
  Modele_Collision_FT& collision_interface_particule();
  const Modele_Collision_FT& collision_interface_particule() const;


  Postraitement_Forces_Interfaces_FT& postraitement_forces_interf();
  const Postraitement_Forces_Interfaces_FT& postraitement_forces_interf() const;

  const double& get_d_to_interf_interp_v() const;
  // fin EB


  virtual DoubleVect calculer_integrale_indicatrice_face(const DoubleVect& indicatrice_face) const; // EB
  virtual DoubleVect calculer_integrale_indicatrice_arete(const DoubleVect& indicatrice_arete) const; // EB
  const Proprietes_part_vol&           proprietes_particules() const;
  const Maillage_FT_Disc&              maillage_inject() const;
  const Proprietes_part_vol&           proprietes_inject() const;

  void nettoyer_proprietes_particules(const ArrOfInt& som_utilises);

  virtual void calculer_vitesse_transport_interpolee(const Champ_base& champ_vitesse,
                                                     const Maillage_FT_Disc& m,
                                                     DoubleTab& vitesse_noeuds,
                                                     long nv_calc) const
  {
    calculer_vitesse_transport_interpolee(champ_vitesse, m, vitesse_noeuds, nv_calc, 1);
  };

  virtual void calculer_vitesse_transport_interpolee(const Champ_base& champ_vitesse,
                                                     const Maillage_FT_Disc&,
                                                     DoubleTab& vitesse_noeuds,
                                                     long nv_calc,
                                                     long standard) const;
  void calculer_scalaire_interpole(const Champ_base& ch_scal,
                                   const Maillage_FT_Disc&,
                                   DoubleTab& ch_scal_noeuds,
                                   long nv_calc) const;

  virtual void remailler_interface();

  //methodes utilisees pour le post-traitement
  virtual long get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, DoubleTab *dtab = 0) const;
  virtual long get_champ_post_FT(const Motcle& champ, Postraitement_base::Localisation loc, IntTab    *itab = 0) const;
  virtual const Maillage_FT_Disc& maillage_interface_pour_post() const;
  long get_mesh_tag() const override
  {
    return maillage_interface_pour_post().get_mesh_tag();
  };

  //Methode d acces au probleme
  const Probleme_base& get_probleme_base() const;

  //Modifie vpoint (pour N_S) pour imposer au fluide la vitesse de l interface
  void modifier_vpoint_pour_imposer_vit(const DoubleTab& inco_val,DoubleTab& vpoint0,
                                        DoubleTab& vpoint,const DoubleTab& rho_faces,
                                        DoubleTab& terme_source,const double temps, const double dt,
                                        const long is_explicite,const double eta) override;

  //Methode outil utilisee par modifier_vpoint_pour_imposer_vit(...)
  void calcul_source(const DoubleTab& inco_val,
                     const DoubleTab& vpoint,
                     const DoubleTab& rho_faces,
                     DoubleTab& source_val,
                     const DoubleTab& vit_imposee,
                     const DoubleTab& indicatrice_faces,
                     const long is_QC,
                     const double dt,
                     const long is_explicite,
                     const double eta);
  void modifie_source(DoubleTab& so_modif,const DoubleTab& so_val,const DoubleTab& rho_faces,
                      const long n,const long m, const long is_QC,
                      const DoubleVect& vol_entrelaces,const Solveur_Masse& solv_masse);

  void calcul_effort_fluide_interface(const DoubleTab& vpoint,const DoubleTab& rho_faces,
                                      DoubleTab& source_val,const long is_explicite,const double eta);

  void impr_effort_fluide_interface( DoubleTab& source_val, DoubleTab& pressure_part, DoubleTab& friction_part ) ;

  void impr_profil_compo_rms_vitesse(Sortie&,DoubleTab& forces_solide, DoubleTab& moy, DoubleTab& moy_carre, DoubleTab& rms) ; // EB

  //Calcul la vitesse imposee a l interface a partir de expression_vitesse_imposee
  virtual void calcul_vitesse(DoubleTab& vitesse_imp, const DoubleTab& champ_vitesse,
                              const DoubleTab& vpoint, const double temps, const double dt);
  virtual void get_expression_vitesse_imposee(DoubleTab& vitesse_imp);
  //Effectue l integration d un ensemble de points (sans notion de facettes)
  void integrer_ensemble_lagrange(const double temps) override;

  virtual void interpoler_vitesse_face(const DoubleTab& distance_interface,
                                       const long phase, const long stencil_width,
                                       DoubleTab& champ, DoubleTab& gradient,
                                       const double t, const double dt ) ;

  void interpoler_simple_vitesse_face(const DoubleTab& distance_interface,
                                      const long phase, const long stencil_width,
                                      DoubleTab& champ, DoubleTab& gradient,
                                      const double t, const double dt ) ;

  virtual void calcul_nb_traverse(  const DoubleTab& xe, const double dx,
                                    const long dim, const long ori,
                                    Maillage_FT_Disc& maillage, long elem,
                                    long& traverse ) ;
  virtual void calcul_OutElemFa7( Maillage_FT_Disc& maillage,
                                  const DoubleTab& indicatrice,
                                  const long nb_elem,
                                  long& nb_fa7_accepted,
                                  IntList& OutElem, IntList& OutFa7 ) ;
  virtual void PPP_face_interface( Maillage_FT_Disc& maillage, const DoubleTab& indicatrice,
                                   const DoubleTab& indicatrice_face, DoubleTab& Vertex ) ;

  virtual void PPP_face_interface_voisin( const DoubleTab& indicatrice, const DoubleTab& indicatrice_face,
                                          DoubleTab& Vertex, DoubleTab& PPP ) ;
  virtual void PPP_face_voisin( const DoubleTab& indicatrice, const DoubleTab& indicatrice_face, DoubleTab& PPP ) ;

  virtual void calcul_maxfa7( Maillage_FT_Disc& maillage, const DoubleTab& indicatrice,
                              const long nb_elem, long& max_fa7, const long exec_planfa7existan ) ;
  virtual void RenumFa7( DoubleTab& Barycentre, DoubleTab& Tab110,DoubleTab& Tab111,
                         DoubleTab& Tab112, IntTab& Tab12, IntTab& CptFacette,
                         const long nb_facettes, const long nb_facettes_dim ) ;
  virtual void StockageFa7( Maillage_FT_Disc& maillage, IntTab& CptFacette, DoubleTab& Tab100,
                            DoubleTab& Tab101,DoubleTab& Tab102, DoubleTab& Tab103,
                            DoubleTab& Tab110,DoubleTab& Tab111,DoubleTab& Tab112,
                            IntTab& Tab12, DoubleTab& Barycentre, const DoubleTab& indicatrice,
                            IntList& OutElem, ArrOfBit& fa7, const long exec_planfa7existant) ;
  virtual void StockageFa7( Maillage_FT_Disc& maillage, DoubleTab& Tab100, DoubleTab& Tab101,
                            DoubleTab& Tab102, DoubleTab& Tab103, DoubleTab& Tab110,
                            DoubleTab& Tab111, DoubleTab& Tab112, IntTab& Tab12,
                            DoubleTab& Barycentre, IntList& OutElem, IntTab& TabOutFa7, ArrOfBit& fa7 ) ;
  virtual void BaryFa7( Maillage_FT_Disc& maillage, const long i_facette, DoubleTab& Barycentre ) ;
  virtual void plan_facette_existant( Maillage_FT_Disc& maillage,
                                      DoubleList A, DoubleList B, DoubleList C,
                                      DoubleList D, const long i_facette,
                                      long& test_liste ) ;
  virtual void calcul_eq_plan_facette(Maillage_FT_Disc& maillage, const long i_facette,
                                      double& a, double& b, double& c, double& d);
  virtual void calcul_eq_plan_facette(const long i_facette,
                                      double& a, double& b, double& c, double& d);

  virtual void calcul_tolerance_projete_monophasique( const long i_face, const long ori, const long voisin0,
                                                      const long voisin1, const DoubleTab& indicatrice_face,
                                                      const DoubleTab& indicatrice, double& tol ) ;

  virtual void calcul_tolerance_projete_diphasique( const long i_face, const long ori, const long voisin0,
                                                    const long voisin1, const DoubleTab& indicatrice, double& tol ) ;

  void verifprojete(const long monophasique,const double Lref, double d, const DoubleTab& x,
                    const DoubleTab& V, DoubleTab& coord_projete, long& cpt ) ;


  virtual void uzawa(const double d,const DoubleTab& matrice, const DoubleTab& x,
                     const DoubleTab& secmem, DoubleTab& solution) const ;

  virtual void projete_point_face_fluide( long& nb_proj_modif, const long dim_fa7,
                                          const DoubleTab& indicatrice_face, const DoubleTab& indicatrice,
                                          const DoubleTab& dist_face, const double t, const double dt,
                                          DoubleTab& Tab100, DoubleTab& Tab101,DoubleTab& Tab102,
                                          DoubleTab& Tab103, IntTab& Tab12, IntTab& CptFacette,
                                          DoubleTab& v_imp, DoubleTab& Vertex,
                                          Parser& parser_x, Parser& parser_y,Parser& parser_z );
  virtual void projete_point_face_interface(   long& nb_proj_modif,const long dim_fa7,
                                               const DoubleTab& indicatrice_face,
                                               const DoubleTab& indicatrice,
                                               const DoubleTab& dist_face, const double t,
                                               const double dt, DoubleTab& Tab100,
                                               DoubleTab& Tab101,DoubleTab& Tab102,
                                               DoubleTab& Tab103, IntTab& Tab12,
                                               IntTab& CptFacette, DoubleTab& v_imp, DoubleTab& Vertex,
                                               Parser& parser_x, Parser& parser_y,Parser& parser_z) ;

  virtual void transporter_sans_changement_topologie(DoubleTab& vitesse,
                                                     const double coeff,const double temps);

  virtual long calculer_composantes_connexes_pour_suppression(IntVect& num_compo);
  virtual double suppression_interfaces(const IntVect& num_compo, const ArrOfInt& flags_compo_a_supprimer);

  const long& get_vimp_regul() const;

  virtual const Champ_base& get_update_distance_interface() const;
  virtual const Champ_base& get_update_distance_interface_faces() const;
  virtual const Champ_base& get_update_normale_interface() const;
  // renvoie DoubleTab parce qu'il n'existe pas de champ aux sommets en VDF ! Zut...
  virtual const DoubleTab&   get_update_distance_interface_sommets() const;
  // idem qu'aux sommets mais au cg des aretes
  virtual const DoubleTab&   get_update_distance_interface_aretes() const; // EB
  void ramasse_miettes(const Maillage_FT_Disc& maillage,
                       DoubleVect& flux,
                       DoubleVect& valeurs);
  void nettoyer_maillage()
  {
    maillage_interface().nettoyer_maillage();
  };
  virtual void calculer_vmoy_composantes_connexes(const Maillage_FT_Disc& maillage,
                                                  const ArrOfInt& compo_connexes_facettes,
                                                  const long nb_compo_tot,
                                                  const DoubleTab& vitesse_sommets,
                                                  DoubleTab& vitesses,
                                                  DoubleTab& positions) const;

  DoubleTab& get_vitesses_compo(); // EB
  DoubleTab& get_positions_compo(); // EB
  DoubleVect& get_rayons_compo(); // EB

  const DoubleTab& get_vitesses_compo() const; // EB
  const DoubleTab& get_positions_compo() const; // EB
  const DoubleVect& get_rayons_compo() const; // EB

  DoubleTab& get_rms_vitesses_compo(); // EB
  DoubleTab& get_moy_vitesses_compo(); // EB
  DoubleTab& get_moy_vitesses_carre_compo(); // EB

  const  DoubleTab& get_rms_vitesses_compo() const; // EB
  const  DoubleTab& get_moy_vitesses_compo() const; // EB
  const  DoubleTab& get_moy_vitesses_carre_compo() const; // EB

  const long& calcul_precis_indic_faces() const;
  const long& calcul_precis_indic_aretes() const;
  const long& postraiter_indicatrice_aretes() const; // EB

  long is_solid_particle(); // EB
  long is_solid_particle() const; // EB
protected:

  void ajouter_contribution_saut_vitesse(DoubleTab& deplacement) const;
  virtual void deplacer_maillage_ft_v_fluide(const double temps);
  void permuter_positions_particules();
  virtual void calculer_distance_interface(const Maillage_FT_Disc& maillage,
                                           DoubleTab& distance_elements,
                                           DoubleTab& normale_elements,
                                           const long n_iter) const;

  virtual void calculer_distance_interface_sommets(const DoubleTab& dist_elem,
                                                   const DoubleTab& normale_elem,
                                                   DoubleTab&        dist_som) const;
  virtual void calculer_distance_interface_aretes(const DoubleTab& dist_elem,
                                                  const DoubleTab& normale_elem,
                                                  DoubleTab&        dist_arete) const; // EB

  virtual void calculer_vitesse_repere_local(const Maillage_FT_Disc& maillage,
                                             DoubleTab& deplacement,
                                             DoubleTab& Positions,
                                             DoubleTab& Vitesses) const;
  virtual void test_suppression_interfaces_sous_domaine();


  virtual void calculer_distance_interface_faces(const DoubleTab& dist_elem,
                                                 const DoubleTab& normale_elem,
                                                 DoubleTab&        dist_faces) const;
  // Nouvelle methodes
  // Methode outil utilisee par modifier_vpoint_pour_imposer_vit(...)
  // Calcul l'indicatrice sur chaque face
  void calcul_indicatrice_faces(const DoubleTab& indicatrice,
                                const IntTab& face_voisins);

  void calcul_indicatrice_aretes(const DoubleTab& indicatrice); // EB
  void postraiter_forces_interface();  // EB



  REF(Probleme_base) probleme_base_;
  REF(Navier_Stokes_FT_Disc) equation_ns_;
  // L'inconnue du probleme
  Champ_Inc indicatrice_;
  Champ_Inc indicatrice_faces_;
  DoubleTab indicatrice_arete_; // EB DoubleVect car il n'existe pas de champ aux aretes
  // Utiliser ces accesseurs :
  Maillage_FT_Disc& maillage_interface();
  // Utiliser ces accesseurs :
  Marching_Cubes& marching_cubes();
  // Utiliser ces accesseurs :
  Topologie_Maillage_FT& topologie_interface();

  Proprietes_part_vol&         proprietes_particules();
  Maillage_FT_Disc&            maillage_inject();
  Proprietes_part_vol&         proprietes_inject();

  DoubleTab& tableaux_positions();
  IntVect& vecteur_elements();
  DoubleTab&    deplacement_som();

  // On utilise des DERIV() pour ne pas avoir a inclure la definition
  // de ces classes (pour reduire les dependances).
  static void transfert_conservatif_eulerien_vers_lagrangien_sommets(const Maillage_FT_Disc& maillage,
                                                                     const DoubleVect& valeurs_euler,
                                                                     ArrOfDouble& valeurs_lagrange);

  Nom suppression_interfaces_sous_domaine_;

  Champ_Fonc vitesse_imp_interp_;
  long calcul_precis_indicatrice_face_; // EB
  long calcul_precis_indicatrice_arete_; // EB
  long postraiter_indicatrice_arete_; // EB
  long get_radius_; // EB

  long get_nb_compo_tot ();


private:

  void init_positions_vitesses_FT();
  void set_nb_compo_tot(const long nb_compo);
  // Variables internes a la methode de transport
  Transport_Interfaces_FT_Disc_interne *variables_internes_;


  double temps_debut_;

  REF(Milieu_base) ref_milieu_;
  REF(Joint) joint_; // EB
  long interpolation_repere_local_;
  long transport_vitesse_cg_HMS_;
  ArrOfDouble force_;
  ArrOfDouble moment_;
};

class Transport_Interfaces_FT_Disc_interne : public Objet_U
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Transport_Interfaces_FT_Disc_interne);
public:
  Transport_Interfaces_FT_Disc_interne() :
    indicatrice_cache_tag(-1),
    indicatrice_face_cache_tag(-1),
    indicatrice_arete_cache_tag(-1),
    fichier_reprise_collision_FT_(""),
    is_solid_particle_(0),
    iterations_correction_volume(0),
    VOFlike_correction_volume(0),
    nb_lissage_correction_volume(0),
    nb_iterations_correction_volume(3),
    volume_impose_phase_1(-1.),
    n_iterations_distance(3),
    n_iterations_interpolation_ibc(5),
    distance_normale_cache_tag(-1),
    distance_sommets_cache_tag(-1),
    distance_faces_cache_tag(-1),
    distance_aretes_cache_tag(-1),
    methode_transport(INDEFINI),
    methode_interpolation_v(VALEUR_A_ELEM),
    d_to_interf_interp_v_(0.5),
    statut_calcul_forces_(0),
    statut_calcul_flux_thermique_(0),
    injection_interfaces_last_time_(0.),
    interpolation_champ_face(BASE),
    vf_explicite(0),
    is_extra_diphasique_solide(0),
    is_extra_solide(0),
    is_distance_projete_face(0),
    nomb_fa7_accepted(3),
    seuil_uzawa(1.e-8),
    nb_iter_uzawa(30),
    vimp_regul(1),
    type_indic_faces_(STANDARD),
    modified_indic_faces_position(0.),
    modified_indic_faces_thickness(1.),
    type_vitesse_imposee(UNIFORME),
    type_distance_calculee(DIST_INITIALE),
    type_projete_calcule(PROJETE_INITIAL),
    expression_vitesse_imposee(Objet_U::dimension),
    precision_impr_(0), // EB
    nb_compo_tot_(0) // EB
  {};
  ~Transport_Interfaces_FT_Disc_interne() override
  {};
private:
  Transport_Interfaces_FT_Disc_interne(const Transport_Interfaces_FT_Disc_interne& a): Objet_U(a)
  {}; // Interdit
  const Transport_Interfaces_FT_Disc_interne& operator=(const Transport_Interfaces_FT_Disc_interne& a)
  {
    return *this;
  }; // Interdit

public:
  long sauvegarder(Sortie& os) const override;
  long reprendre(Entree& is) override;

  // Les membres suivantes sont sauvegardes et repris:
  Champ_Inc        indicatrice_cache;     // L'indicatrice calculee par get_update_indicatrice
  long           indicatrice_cache_tag; // Le tag du maillage correspondant
  Champ_Inc        indicatrice_face_cache; // EB
  long 			indicatrice_face_cache_tag; // EB
  DoubleTab        indicatrice_arete_cache; // EB
  long 			indicatrice_arete_cache_tag; // EB
  Maillage_FT_Disc maillage_interface;          // Objet qui peut se reduire a un ensemble de sommets
  // quand il represente les positions de particules
  Remaillage_FT    remaillage_interface_;
  // Fin des membres sauvegardes / repris

  Proprietes_part_vol proprietes_particules_; //Proprietes physiques de particules
  Nom fichier_reprise_collision_FT_; // EB
  long is_solid_particle_; // EB
  Modele_Collision_FT collision_interface_particule_; // EB

  Postraitement_Forces_Interfaces_FT postraitement_forces_interf_; // EB

  Maillage_FT_Disc maillage_inject_;              //Ensemble de particules a injecter periodiquement
  Proprietes_part_vol proprietes_inject_;     //Proprietes physiques des particules injectees

  // FORMER Keyword (obsolete):
  // Si iterations_correction_volume == 0, le maillage est transporte par le fluide
  // a l'aide du champ de vitesse L2 interpole (pas de conservation du volume).
  // Si iterations_correction_volume > 0, on calcule une correction de volume aux
  // sommets de l'interface et on l'etale par autant d'iterations d'un lisseur
  // Voir Transport_Interfaces_FT_Disc::mettre_a_jour
  long iterations_correction_volume;

  // NEW Keywords/parameters for volume preserving correction in agreement with phase change :
  long VOFlike_correction_volume;
  // Si VOFlike_correction_volume == 0, le maillage est transporte par le fluide
  // a l'aide du champ de vitesse L2 interpole (pas de conservation du volume).
  // Si VOFlike_correction_volume > 0, on calcule une correction de volume aux
  // sommets de l'interface et on l'etale par nb_lissage_correction_volume iterations d'un lisseur,
  // Voir Transport_Interfaces_FT_Disc::mettre_a_jour
  // pour eviter l'apparition de pic aux interfaces (ie, on lisse legement la correction de volume
  // (uniquement s'il y en a une))
  long nb_lissage_correction_volume;
  // La correction est iterative car on ne corrige pas exactement du volume demande en deplacant les
  // noeuds sequentiellement. On fait nb_iterations_correction_volume ou jusqu'a ce que l'erreur
  // soit inferieure au seuil de correction de volume de Remaillage_FT
  long nb_iterations_correction_volume;

  // ADDITIONAL GLOBAL mass conservation with-out phase change:
  // Rustine introduite pour corriger les pertes de masse:
  // Si cette valeur est positive, on deplace toute l'interface d'une certaine distance
  // pour que le volume de la phase 1 reste toujours egal a la valeur prescrite
  double volume_impose_phase_1;
  // Si non nul, le calcul de l'integrale pour le volume impose porte sur cette sous-domaine
  Nom    nom_domaine_volume_impose_;

  Champ_Inc vitesse_filtree;
  DoubleTab doubletab_pos;
  DoubleTab doubletab_vitesses;
  IntVect   intvect_elements;

  Champ_Fonc distance_interface; // Distance a l'interface (aux elements)
  Champ_Fonc normale_interface;  // Une normale etalee
  Champ_Fonc surface_interface;  // GB : La surface d'interface dans chaque element
  Champ_Fonc tmp_flux;           // Tableau temporaire pour le ramasse-miettes
  Champ_Fonc distance_interface_faces; // CF : Distance a l'interface (aux faces)
  DoubleTab  distance_interface_sommets; // Distance a l'interface (aux sommets)
  DoubleTab  distance_interface_aretes; // Distance a l'interface (aux sommets) // EB
  Champ_Fonc distance_interface_faces_corrigee; // CI : Distance a l'interface corrigee (aux faces)
  Champ_Fonc distance_interface_faces_difference; // CI : Distance a l'interface corrigee - Distance a l'interface calculee (aux faces)
  Champ_Fonc index_element; // CI : indexation des elements
  Champ_Fonc nelem_par_direction; // CI : nombre d'elements par direction
  // Note de B.M. : zut, y'a pas de champ aux sommets en VDF... donc DoubleTab et
  // du coup on ne peut pas postraiter facilement. C'est trop con...
  long     n_iterations_distance;
  long     n_iterations_interpolation_ibc;
  long     distance_normale_cache_tag;
  long     distance_sommets_cache_tag;
  long     distance_faces_cache_tag;
  long     distance_aretes_cache_tag; // EB
  // Le maillage postraite n'est pas forcement le maillage
  Maillage_FT_Disc maillage_pour_post;
  DoubleTabFT   deplacement_sommets;

  enum Methode_transport { INDEFINI, VITESSE_IMPOSEE, LOI_HORAIRE, VITESSE_INTERPOLEE };
  Methode_transport      methode_transport;
  REF(Navier_Stokes_std) refequation_vitesse_transport;
  REF(Convection_Diffusion_Temperature_FT_Disc) refequation_temperature_;
  enum Methode_interpolation_v { VALEUR_A_ELEM, VDF_LINEAIRE, VITESSE_SOLIDE_MOYENNE, VITESSE_SOLIDE_SOMMETS }; // EB : rajout de VITESSE_SOLIDE_MOYENNE et VITESSE_SOLIDE_SOMMETS
  Methode_interpolation_v methode_interpolation_v;

  double d_to_interf_interp_v_; // EB : distance a l'interface pour l'interpolation de la vitesse DANS la particule, en nombre de rayon de la particule  0 on calcul, 1 on ne calcule pas
  long statut_calcul_forces_; // EB : indique si il faut recalculer les forces hydrodynamiques sur les facettes du maillage lagrangien : 0 on calcul, 1 on ne calcule pas
  long statut_calcul_flux_thermique_; // EB idem pour le flux

  // Injecteur d'interfaces au cours du temps
  ArrOfDouble injection_interfaces_temps_;
  ArrOfInt    injection_interfaces_phase_;
  Noms        injection_interfaces_expressions_;
  // temps physique reel auquel a eu lieu la derniere injection (ce n'est pas le temps demande)
  double      injection_interfaces_last_time_;

  enum Interpolation_champ_face { BASE, LINEAIRE };
  Interpolation_champ_face interpolation_champ_face;

  long vf_explicite, is_extra_diphasique_solide, is_extra_solide, is_distance_projete_face;
  long nomb_fa7_accepted ;
  double seuil_uzawa ;
  long nb_iter_uzawa ;
  long vimp_regul ;

  enum Type_indic_faces { STANDARD, MODIFIEE, AI_BASED };
  Type_indic_faces type_indic_faces_;
  double modified_indic_faces_position ;
  double modified_indic_faces_thickness ;

  enum Type_vitesse_imposee { UNIFORME, ANALYTIQUE };
  Type_vitesse_imposee type_vitesse_imposee;

  enum Type_distance_calculee { DIST_INITIALE, DIST_MODIFIEE };
  Type_distance_calculee type_distance_calculee;

  enum Type_projete_calcule { PROJETE_INITIAL, PROJETE_MODIFIE, PROJETE_SIMPLIFIE };
  Type_projete_calcule type_projete_calcule;

  Noms                      expression_vitesse_imposee;
  // Reference a une loi horaire eventuelle
  REF(Loi_horaire)         loi_horaire_;

  // Integration de la vitesse a partir du point de depart (x,y,z)
  //  pendant un temps dt.
  void integrer_vitesse_imposee(double temps, double dt, double& x, double& y, double& z) const;

  // Les objets-algorithmes :
  Parcours_interface      parcours_interface_;
  Marching_Cubes          marching_cubes_;
  Connectivite_frontieres connectivite_frontieres_;
  Topologie_Maillage_FT   topologie_interface_;
  // Cet objet est type en fonction de la discretisation:
  DERIV(Algorithmes_Transport_FT_Disc) algorithmes_transport_;

// HMS : variable pour stocker les vitesses et les positions du centre de gravite pour chaque inclusion
  DoubleTab vitesses_compo;
  DoubleTab positions_compo;
  DoubleVect rayons_compo;


// EB
  long precision_impr_;
  long nb_compo_tot_;
  DoubleTab rms_vitesse_compo_;
  DoubleTab moy_vitesse_compo_;
  DoubleTab moy_vitesse_solide_carre_;


};
#endif
