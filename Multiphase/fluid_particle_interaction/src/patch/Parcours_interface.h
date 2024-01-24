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
// File:        Parcours_interface.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Parcours_interface_included
#define Parcours_interface_included

#include <ArrOfBit.h>
#include <FTd_tools.h>
#include <Maillage_FT_Disc.h>
#include <TRUSTTabs_forward.h>
#include <TRUST_Ref.h>
#include <Domaine_VDF.h> // EB

class Zone_VF;
class Maillage_FT_Disc;
class Domaine;
class Domaine_dis;
class Connectivite_frontieres;

class Parcours_interface : public Objet_U
{
  Declare_instanciable_sans_constructeur(Parcours_interface);
public:
  Parcours_interface();
  void associer_domaine_dis(const Domaine_dis& domaine_dis);
  void associer_connectivite_frontieres(const Connectivite_frontieres& connect);
  long calculer_face_sortie_element(const Domaine_VF& domaine_vf,
                                    const long num_element,
                                    double x0, double y0, double z0,
                                    double x1, double y1, double z1,
                                    double& pos_intersection) const;
  long calculer_sortie_face_bord(const long face_0,
                                 const long num_element,
                                 double x0, double y0, double z0,
                                 double x1, double y1, double z1,
                                 double& x, double& y, double& z) const;
  double distance_sommet_faces(const Domaine_VF& domaine_vf,
                               const long num_element,
                               double x, double y, double z) const;

  double distance_sommet_faces_face(const Domaine_VF& domaine_vf, const long num_face, double x, double y, double z) const; // EB

  double get_erreur_geometrique() const;
  long get_correction_parcours_thomas() const
  {
    return correction_parcours_thomas_ ;
  };
  void projeter_vecteur_sur_face(const long num_face, double& x_, double& y_, double& z_) const;

  void calculer_normale_face_bord(long num_face, double x, double y, double z,
                                  double& nx_, double& ny_, double& nz_) const;
  void remplir_equation_plan_faces_aretes_internes(const Domaine_dis& domaine_dis); // EB

protected:
  // Ceci est le point d'entree pour le parcours de l'interface.
  // Appel autorise uniquement a Maillage_FT_Disc::parcourir_maillage()
  void parcourir(Maillage_FT_Disc& maillage) const;
  friend void Maillage_FT_Disc::parcourir_maillage();

  Parcours_interface(const Parcours_interface&);  // Interdit !
  const Parcours_interface& operator=(const Parcours_interface&);   // Interdit !

  void parcours_facette_face_x(const Domaine_VF& domaine_vf,
                               Maillage_FT_Disc& maillage,
                               ArrOfInt& echange_facettes_numfacette,
                               ArrOfInt& echange_facettes_numface,
                               long num_facette,
                               long face_depart, Intersections_Face_Facettes& intersections_face,
                               ArrOfBit& drapeaux_faces_parcourues) const;   // EB

  void parcours_facette_face_y(const Domaine_VF& domaine_vf,
                               Maillage_FT_Disc& maillage,
                               ArrOfInt& echange_facettes_numfacette,
                               ArrOfInt& echange_facettes_numface,
                               long num_facette,
                               long face_depart, Intersections_Face_Facettes& intersections_face,
                               ArrOfBit& drapeaux_faces_parcourues) const;   // EB
  void parcours_facette_face_z(const Domaine_VF& domaine_vf,
                               Maillage_FT_Disc& maillage,
                               ArrOfInt& echange_facettes_numfacette,
                               ArrOfInt& echange_facettes_numface,
                               long num_facette,
                               long face_depart, Intersections_Face_Facettes& intersections_face,
                               ArrOfBit& drapeaux_faces_parcourues) const;   // EB

  void parcours_facette_arete_x(const Domaine_VF& domaine_vf,
                                Maillage_FT_Disc& maillage,
                                ArrOfInt& echange_facettes_numfacette,
                                ArrOfInt& echange_facettes_numarete,
                                long num_facette,
                                long arete_depart, Intersections_Arete_Facettes& intersections_arete,
                                ArrOfBit& drapeaux_aretes_parcourues) const;   // EB
  void parcours_facette_arete_y(const Domaine_VF& domaine_vf,
                                Maillage_FT_Disc& maillage,
                                ArrOfInt& echange_facettes_numfacette,
                                ArrOfInt& echange_facettes_numarete,
                                long num_facette,
                                long arete_depart, Intersections_Arete_Facettes& intersections_arete,
                                ArrOfBit& drapeaux_aretes_parcourues) const;   // EB
  void parcours_facette_arete_z(const Domaine_VF& domaine_vf,
                                Maillage_FT_Disc& maillage,
                                ArrOfInt& echange_facettes_numfacette,
                                ArrOfInt& echange_facettes_numarete,
                                long num_facette,
                                long arete_depart, Intersections_Arete_Facettes& intersections_arete,
                                ArrOfBit& drapeaux_aretes_parcourues) const;   // EB
  void parcours_facette(const Domaine_VF& domaine_vf,
                        Maillage_FT_Disc& maillage,
                        ArrOfInt& echange_facettes_numfacette,
                        ArrOfInt& echange_facettes_numelement,
                        long num_facette,
                        long element_depart) const;

  long calcul_intersection_facelem_2D(const Domaine_VF& domaine_vf,
                                      Maillage_FT_Disc& maillage,
                                      long num_facette,
                                      long num_element) const;

  long calcul_intersection_facelem_3D(const Domaine_VF& domaine_vf,
                                      Maillage_FT_Disc& maillage,
                                      long num_facette,
                                      long num_element) const;
  long calcul_intersection_facette_face_3D(const Domaine_VF& domaine_vf,
                                           Maillage_FT_Disc& maillage,
                                           long num_facette,
                                           long num_face, Intersections_Face_Facettes& intersections_face) const; // EB
  long calcul_intersection_facette_arete_3D(const Domaine_VF& domaine_vf,
                                            Maillage_FT_Disc& maillage,
                                            long num_facette,
                                            long num_arete, Intersections_Arete_Facettes& intersections_face) const; // EB

  double calcul_eq_plan(const Domaine_VF& domaine_vf,
                        const long num_element, const long num_face_element,
                        double& a, double& b, double& c, double& d) const;
  double calcul_eq_plan_face(const Domaine_VF& domaine_vf, const long num_face, const long num_face_face, double& a, double& b, double& c, double& d) const; // EB
  double calcul_eq_aretes(const Domaine_VDF& domaine_vdf, const long num_arete, const long num_face_arete, double& a, double& b, double& c, double& d) const; // EB

  double volume_rectangle(const Domaine_VF& domaine_vf, long num_element,
                          double x0, double y0, double x1, double y1,
                          double epsilon) const;

  // New function to get phase-barycenter in 2D or 2D-axi calculations.
  // 3D equivalent not yet available (2020/10/26)
  double volume_barycentre_rectangle(const Domaine_VF& domaine_vf, long num_element,
                                     double x0, double y0, double x1, double y1,
                                     double epsilon, double liquid_barycentre[3]) const;

  double volume_triangle(const Domaine_VF& domaine_vf, long num_element,
                         double x0, double y0, double x1, double y1,
                         long plan_coupe0, long plan_coupe1) const;

  double volume_tetraedre(const Domaine_VF& domaine_vf,
                          long num_element,
                          long num_facette,
                          const Maillage_FT_Disc& maillage,
                          const DoubleTab& poly_reelles,
                          const FTd_vecteur3& centre_de_gravite,
                          double epsilon) const;

  double volume_tetraedre_reference(const DoubleTab& poly_reelles_ref,
                                    const FTd_vecteur3& norme_ref,
                                    const FTd_vecteur3& centre_de_gravite_ref,
                                    double epsilon) const;

  double volume_hexaedre(const Domaine_VF& domaine_vf, long num_element,
                         const DoubleTab& poly_reelles,
                         const FTd_vecteur3& norme,
                         const FTd_vecteur3& centre_de_gravite,
                         const ArrOfInt& polygone_plan_coupe,
                         double epsilon) const;
// debut EB
  double volume_hexaedre_face(const Domaine_VDF& domaine_vdf, long num_face,
                              const DoubleTab& poly_reelles,
                              const FTd_vecteur3& norme,
                              const FTd_vecteur3& centre_de_gravite,
                              const ArrOfInt& polygone_plan_coupe,
                              double epsilon) const;

  double volume_hexaedre_arete_interne(const Domaine_VDF& domaine_vdf, long num_arete,
                                       const DoubleTab& poly_reelles,
                                       const FTd_vecteur3& norme,
                                       const FTd_vecteur3& centre_de_gravite,
                                       const ArrOfInt& polygone_plan_coupe,
                                       double epsilon) const;
// fin EB

  void matrice_triangle(long num_element,
                        FTd_vecteur2& origine,
                        FTd_matrice22& matrice,
                        double& surface) const;

  void transformation_2d(const FTd_vecteur2& origine,
                         const FTd_matrice22& matrice,
                         double x, double y,
                         double& u, double& v) const;

  static void calcul_inverse_matrice33(const FTd_matrice33& matrice, FTd_matrice33& matrice_inv);
  static void calcul_produit_matrice33_vecteur(const FTd_matrice33& matrice, const FTd_vecteur3& vect, FTd_vecteur3& res);

  // Variables persistantes de la classe :
  REF(Domaine_VF) refdomaine_vf_;
  REF(Connectivite_frontieres) refconnect_front_;
  long nb_faces_elem_;
  long nb_faces_reelles_;
  long nb_aretes_reelles_; // EB
  long nb_elements_reels_;
  long nb_sommets_par_face_;
  enum { TRIANGLE, RECTANGLE, TETRA, HEXA } type_element_;

  // Une tableau de taille nb_faces * 4 contenant les coefficients du
  // plan contenant la face sous la forme
  //  a = equation_plans_faces_(num_face, 0)
  //  b =                                 1
  //  c =                                 2  = 0. en 2D.
  //  d =                                 3
  //  f(x,y,z) = a*x + b*y + c*z + d
  // et le vecteur (a,b,c) unitaire colineaire a la normale, dirige
  // de l'element face_voisins(num_face, 0) vers l'element face_voisins(num_face, 1)
  DoubleTabFT equations_plans_faces_;
  DoubleTabFT equations_plans_faces_face_;
  DoubleTabFT equations_plans_faces_arete_;
  //
  // Variables temporaires utilisees dans l'algorithme de parcours
  //
  // A des fins de statistiques (voir "calcul_intersection_facelem_*"
  // et "parcourir" )
  mutable long compteur_erreur_grossiere;

  // Raccourcis vers les elements et les sommets du maillage eulerien
  mutable const IntTab * domaine_elem_ptr;
  mutable const DoubleTab * domaine_sommets_ptr;

  // Marqueurs des elements deja visites :
  mutable ArrOfBit drapeaux_elements_parcourus_;
  // Marqueurs des faces deja visites : // EB
  mutable ArrOfBit drapeaux_faces_parcourues_x_; // EB
  mutable ArrOfBit drapeaux_faces_parcourues_y_; // EB
  mutable ArrOfBit drapeaux_faces_parcourues_z_; // EB

  mutable ArrOfBit drapeaux_aretes_parcourues_x_; // EB
  mutable ArrOfBit drapeaux_aretes_parcourues_y_; // EB
  mutable ArrOfBit drapeaux_aretes_parcourues_z_; // EB

  // Si on constate un probleme de precision numerique, la correction
  // suppose que les calculs geometriques ont la precision relative
  // suivante:
  static const double Erreur_relative_maxi_;
  // On suppose que toutes les coordonnees du domaine sont inferieures a cette valeur:
  double Valeur_max_coordonnees_;
  // Cette valeur est egale a Erreur_relative_maxi_ * valeur_max_coordonnees_ :
  double Erreur_max_coordonnees_;

  // Drapeau d'activation de la correction du parcours (par Thomas Fortin,
  //  pour corriger les problemes lies aux sommets qui tombent sur les faces du maillages eulerien
  long correction_parcours_thomas_;
  long eloigner_sommets_des_faces(Maillage_FT_Disc& maillage) const;
  double uzawa2(const Domaine_VF& domaine_vf, const long elem,
                double& x, double& y, double& z,  const long is_face=0) const;
};

#endif