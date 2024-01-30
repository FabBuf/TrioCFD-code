/////////////////////////////////////////////////////////////////////////////
//
// File      : Modele_Collision_FT.cpp
// Directory : $FPI184_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#ifndef Modele_Collision_FT_included
#define Modele_Collision_FT_included

#include <TRUSTTabFT_forward.h>
#include <Objet_U.h>
#include <TRUST_Deriv.h>
#include <TRUST_Ref.h>
#include <FTd_tools.h>
#include <Pave.h>
#include <Particule_Solide.h>
#include <Fluide_Diphasique.h>
//#include <Transport_Interfaces_FT_Disc.h>
#include <Domaine_VF.h>
#include <Domaine_VDF.h>
#include <MD_Vector.h>
#include <Domaine.h>
#include <TRUSTTabFT.h>
class Param;
class Maillage_FT_Disc;
class Transport_Interfaces_FT_Disc;
// ====================================================================
// .DESCRIPTION        : class Modele_Collision_FT
//  Cette classe implemente le calcul de la force de contact
//
// .SECTION voir aussi
//  Transport_Interfaces_FT_Disc
#include <MD_Vector_tools.h>
//#include <iostream> // EB
//#include <vector> // EB

class Modele_Collision_FT : public Objet_U
{
  Declare_instanciable_sans_constructeur(Modele_Collision_FT);

public:
  Modele_Collision_FT();
  void set_param(Param& p);
  long lire_motcle_non_standard(const Motcle&, Entree&) override;
  long reprendre(Entree& is) override;
  long sauvegarder(Sortie& os) const override;
  void reset(); // il faut que les tableaux aient les bonnes dimensions pour etre lu lors de la reorise
  void calculer_force_contact(DoubleTab& force_contact, long& isFirstStepOfCollision, double& dist_int, double& next_dist_int, DoubleTab& norm, DoubleTab& dUn, double& masse_eff, long& compo, long& voisin, double& Stb, double& ed, double& vitesseRelNorm, double& dt, double& prod_scal);
  const double& sigma() const;
  const double& tau_coll() const;
  const double& delta_n() const;
  const double& get_raideur_cst() const;
  const double& get_amortissement_cst() const;

  void set_s_Verlet(double s_Verlet);
  const long& is_detection_Verlet() const;
  const long& is_LC_activated() const;
  double& get_s_Verlet();
  const long& get_Px() const;
  const long& get_Py() const;
  const long& get_Pz() const;
  long& get_nb_dt_Verlet();
  long& get_dt_compute_Verlet();
  long& get_nb_pas_dt_max_Verlet();
  const double& get_d_act_lub() const;
  const double& get_d_sat_lub() const;
  ArrOfIntFT& get_liste_zone_sup();
  const ArrOfIntFT& get_liste_zone_sup() const;
  ArrOfIntFT& get_liste_zone_inf();
  const ArrOfIntFT& get_liste_zone_inf() const;


  DoubleTab& get_raideur();
  DoubleTab& get_e_eff();
  DoubleTab& get_F_old();
  DoubleTab& get_F_now();
  DoubleTab& get_forces_solide();

  DoubleVect& get_collisions_detected();

  long& compteur_collisions();
  double& rayon_particule();
  const long& modele_lubrification() const;
  const long& force_elem_diphasique() const;
  const DoubleTab& position_bords() const;
  static void set_resize_parametres_geometriques();
  static void set_longueur(DoubleVect& Longueurs);
  static void set_longueur(double Lx, double Ly, double Lz);
  void set_param_geom(Domaine_VDF& domaine_vdf);
  static void set_origin(DoubleVect& Origin);
  static void set_origin(double Ox, double Oy, double Oz);

  static DoubleVect& get_origin();
  static DoubleVect& get_longueurs();
  static IntVect& get_nb_noeuds();
  void calculer_positions_bords(const DoubleVect& rayon_particule);
  long checkForDuplicates(ArrOfInt& vector);

  void associer_equation_transport(const Equation_base& equation);
  void set_nb_compo_tot(long nb_compo_tot);
  void set_nom_fichier_reprise_FT(Nom fichier_reprise_FT);
  void set_d_act_lub(double d_act_lub);
  void set_d_sat_lub(double d_sat_lub);

protected:
  double tau_coll_;
  long decalage_bords_;
  long f_elem_diph_;
  long modele_lubrification_;
  double sigma_;
  double delta_n_;
  double raideur_cst_;
  double amortissement_cst_;
  double d_act_lub_;
  double d_sat_lub_;
  long nb_compo_tot_;
  //long sauvegarde_reprise_fichier_unique_;
  Nom fichier_reprise_FT_;
  DoubleTab raideur_;
  DoubleTab e_eff_;
  DoubleTab F_old_;
  DoubleTab F_now_;
  DoubleTab forces_solide_;
  DoubleTab positions_bords_;
  DoubleVect collision_detected_;

  double s_Verlet_;
  long nb_dt_Verlet_;
  long dt_compute_Verlet_;
  long nb_pas_dt_max_Verlet_;
  long is_detection_Verlet_;
  long activate_linked_cell_;
  long Px_; // nb procs suivant x
  long Py_; // nb procs suivant y
  long Pz_; // nb procs suivant z
  ArrOfIntFT liste_zone_sup_;
  ArrOfIntFT liste_zone_inf_;


  DoubleVect valeurs_decalage;
  long compteur_collisions_;

  REF(Transport_Interfaces_FT_Disc) refequation_transport_;

private:
  REF(Domaine) ref_domaine;
  static DoubleVect Longueurs_modele_collision;
  static IntVect Nb_Noeuds_modele_collision;
  static DoubleVect Origine_modele_collision;

  enum Modele_collision { RESSORT_AMORTI_VVA, RESSORT_AMORTI_ESI, RESSORT_AMORTI_EE, MOHAGHEGH, HYBRID_EE, HYBRID_ESI, HYBRID_VVA, BREUGEM };
  Modele_collision modele_collision_;

};




#endif
