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

#ifndef Maillage_FT_Disc_included
#define Maillage_FT_Disc_included

#include <Ensemble_Lagrange_base.h>
#include <TRUSTTabFT.h>
#include <Descripteur_FT.h>
#include <Intersections_Elem_Facettes_Data.h>
#include <Intersections_Face_Facettes_Data.h> // EB
#include <Intersections_Arete_Facettes_Data.h> // EB
#include <TRUSTTabs_forward.h>
#include <TRUST_Deriv.h>
#include <TRUST_Ref.h>

class Domaine_dis;
class Transport_Interfaces_FT_Disc;
class Remaillage_FT;
class Topologie_Maillage_FT;
class Parcours_interface;
class Maillage_Echange;
class Domaine_VF;
class Maillage_FT_Disc_Data_Cache;

/*! @brief : class Maillage_FT_Disc Cette classe decrit un maillage:
 *
 *    un tableau de coordonnees des sommets,
 *    un tableau de facettes,
 *    drapeaux,
 *    intersections facettes / elements
 *    ...
 *
 */
class Maillage_FT_Disc : public Ensemble_Lagrange_base
{
  Declare_instanciable_sans_constructeur(Maillage_FT_Disc);
public:
  enum Statut_Maillage { RESET = 0, MINIMAL = 1, PARCOURU = 2, COMPLET = 3};
  //
  // Methodes reimplementees de Objet_U
  //
  long reprendre(Entree&) override;
  long sauvegarder(Sortie&) const override;
  //
  // Nouvelles methodes
  //
  Maillage_FT_Disc();
  virtual Entree& lire_param_maillage(Entree& is);
  virtual void ajouter_maillage(const Maillage_FT_Disc& maillage_tmp,long skip_facettes=0);
  virtual void recopie(const Maillage_FT_Disc& maillage_source, Statut_Maillage niveau_copie);
  const Equation_base& equation_associee() const override;

  // Proprietes pour le schema en temps...
  double temps() const;
  double changer_temps(double t);

  // Acces au tag :
  long get_mesh_tag() const;

  // Acces aux elements du maillage
  const DoubleTab& sommets() const;
  long              nb_sommets() const;  // Egal a sommets().dimension(0)
  const DoubleTab& cg_fa7() const; // EB
  const IntTab&     facettes() const;
  long              nb_facettes() const; // Egal a facettes().dimension(0)
  const ArrOfInt& drapeaux_sommets() const;  // pour postraitement uniquement
  const ArrOfInt& sommet_PE_owner() const;   // pour postraitement uniquement
  const ArrOfInt& sommet_num_owner() const;   // pour postraitement uniquement
  void               facette_PE_owner(ArrOfInt& facette_pe) const; // pour postraitement uniquement
  const ArrOfInt& sommet_elem() const;       // pour postraitement uniquement
  const IntTab& sommet_face() const;       // pour postraitement uniquement // EB
  const ArrOfInt& sommet_face_bord() const;  // pour postraitement uniquement

  const Desc_Structure_FT& desc_sommets() const;
  const Desc_Structure_FT& desc_facettes() const;

  const Intersections_Elem_Facettes& intersections_elem_facettes() const;

  const ArrOfInt& som_init_util() const;
  // Ces fonctions renvoient 1 si le test est vrai, 0 sinon
  inline long sommet_virtuel(long i) const;
  inline long sommet_virtuel_face(long i) const; // EB
  inline long sommet_ligne_contact(long i) const;
  inline long sommet_face_bord(long i) const;
  inline long facette_virtuelle(long i) const;
  inline long facette_virtuelle_face(long i) const; // EB
  long        facettes_voisines(long fa70, long fa71, long& iarete0, long& iarete1) const;

  long calculer_voisinage_facettes(IntTab& fa7Voisines,
                                   const Intersections_Elem_Facettes* ief=NULL) const;

  void associer_equation_transport(const Equation_base& equation) override;
  void associer_domaine_dis_parcours(const Domaine_dis& domaine_dis, const Parcours_interface& parcours);
  Transport_Interfaces_FT_Disc& equation_transport();
  const Transport_Interfaces_FT_Disc& equation_transport() const;
  long type_statut() const;

  // Vide le maillage
  void reset();
  // Amene le maillage dans l'etat PARCOURU
  void parcourir_maillage();
  void remplir_equation_plan_faces_aretes_internes(Domaine_dis& zone_dis); // EB
  // Amene le maillage dans l'etat COMPLET
  void completer_maillage();
  // Calcul de l'indicatrice de phase
  void calcul_indicatrice(DoubleVect& indicatrice, const DoubleVect& indicatrice_precedente);
  void calcul_indicatrice_face(const DoubleVect& indicatrice, DoubleVect& indicatrice_face, const DoubleVect& indicatrice_face_precedente);
  void calcul_indicatrice_arete(const DoubleVect& indicatrice, DoubleVect& indicatrice_arete, const DoubleVect& indicatrice_arete_precedente);
  // Deplacement des sommets du maillage selon un vecteur donne
  virtual void transporter(const DoubleTab& deplacement);

  // Deplacement simple des sommets du maillage selon un vecteur donne
  // (sans consideration sur les facettes)
  void transporter_simple(const DoubleTab& deplacement);
  //Nettoyage des sommets virtuels et des sommets sur frontiere ouverte
  void nettoyer_noeuds_virtuels_et_frontieres();
  void nettoyer_phase(const Nom& nom_eq, const long phase);

  void nettoyer_elements_virtuels();
  virtual void nettoyer_maillage();
  void supprimer_facettes(const ArrOfInt& liste_facettes);

  //ecrit le maillage de facon a etre relu sous gnuplot
  void ecrire_plot(const Nom& nom,double temps, long niveau_requete) const;
  //fonctions d'affichage standard des facettes et sommets
  Sortie& printFa7(long fa7,long affsom, Sortie& os) const;
  Sortie& printSom(long som,Sortie& os) const;

  enum AjoutPhase { AJOUTE_TOUTES_PHASES = -1, AJOUTE_PHASE0 = 0 ,AJOUTE_PHASE1  = 1 };

  // Acces aux grandeurs calculees a partir du maillage:
  virtual const ArrOfDouble& get_update_surface_facettes() const;
  virtual const DoubleTab& get_update_normale_facettes() const;
  virtual double compute_normale_element(const long elem, const bool normalize, ArrOfDouble&  normale) const;
  virtual double compute_surface_and_normale_element(const long elem, const bool normalize, double surface, double normal[3]) const;
  virtual const ArrOfDouble& get_update_courbure_sommets() const;

  virtual const ArrOfDouble& get_surface_facettes() const;
  virtual const DoubleTab& get_normale_facettes() const;

  inline long set_niveau_plot(long niv);


  //fonctions statiques (pour ne pas utiliser des membres propres au maillage
  static long deplacer_un_point(double& x, double& y, double& z,
                                double x1, double y1, double z1,
                                long& element_suivant,
                                long& face_suivante,
                                const Parcours_interface& parcours,
                                const Domaine_VF& domaine_vf,
                                const IntTab& face_voisins,
                                long skip_facettes=0);

  static long deplacer_un_sommet(double& x, double& y, double& z,
                                 double x1, double y1, double z1,
                                 long& element,
                                 long& face_bord,
                                 const long num_sommet,
                                 const Parcours_interface& parcours,
                                 const Domaine_VF& domaine_vf,
                                 const IntTab& face_voisins,
                                 ArrOfInt& sommets_envoyes,
                                 ArrOfInt& element_virtuel_arrivee,
                                 ArrOfInt& face_virtuelle_arrivee,
                                 DoubleTab& deplacement_restant,
                                 long skip_facettes=0);

  void update_sommet_face(); // EB
  void update_sommet_arete(); // EB
  // Utiliser ces deux fonctions pour recuperer des champs valides apres avoir
  // transporte le maillage :
  void preparer_tableau_avant_transport(ArrOfDouble& tableau,
                                        const Desc_Structure_FT& descripteur) const;
  void preparer_tableau_avant_transport(ArrOfInt& tableau,
                                        const Desc_Structure_FT& descripteur) const;
  void preparer_tableau_avant_transport(DoubleTab& tableau,
                                        const Desc_Structure_FT& descripteur) const;
  void update_tableau_apres_transport(ArrOfDouble& tableau, long nb_elements,
                                      const Desc_Structure_FT& descripteur) const;
  void update_tableau_apres_transport(ArrOfInt& tableau, long nb_elements,
                                      const Desc_Structure_FT& descripteur) const;
  void update_tableau_apres_transport(DoubleTab& tableau, long nb_elements,
                                      const Desc_Structure_FT& descripteur) const;

  static double angle_bidim_axi();

  Maillage_FT_Disc(const Maillage_FT_Disc&);
  void generer_structure();
  void remplir_structure(const DoubleTab& soms);

  void creer_tableau_sommets(Array_base&, Array_base::Resize_Options opt = Array_base::COPY_INIT) const;
  void creer_tableau_elements(Array_base&, Array_base::Resize_Options opt = Array_base::COPY_INIT) const;
  double calcul_normale_3D(long num_facette, double norme[3]) const;
  void calcul_cg_fa7(); //EB
  virtual void   calculer_costheta_minmax(DoubleTab& costheta) const;
  Schema_Comm_FT get_schema_comm_FT() const;
  void set_is_solid_particle(long is_solid_particle); // EB

protected:
  void pre_lissage_courbure(ArrOfDouble& store_courbure_sommets, const long niter) const;
  void correction_costheta(const double c, const long s0, const long facette,
                           /* const ArrOfDouble& s0s1, const ArrOfDouble& s0s2, */
                           double ps) const;
  double calculer_costheta_objectif(const long som, const long facette, const long call, const double c,
                                    const DoubleTabFT& tab_cos_theta, ArrOfBit& drapeau_angle_in_range) const;

  const Maillage_FT_Disc& operator=(const Maillage_FT_Disc&);   // Interdit !

  // Cette methode permet de declarer que le maillage a ete
  // modifie (et invalider le cache de valeurs calculees)
  void maillage_modifie(Statut_Maillage nouveau_statut);

  void calculer_voisins();

  void calcul_surface_normale(ArrOfDouble& surface, DoubleTab& normale) const;
  void calcul_courbure_sommets(ArrOfDouble& courbure, const long call_number) const;

  //fonction qui cree un nouveau sommet par copie d'une existant
  //utilise dans Remailleur_Collision_FT_Collision_Seq
  long copier_sommet(long som);
  //idem, mais rend le sommet interne (sommet_face_bord_=-1)
  long copier_sommet_interne(long som);

  //renvoie les coefficients du plan passant par l'arete iarete de la facette fa7 et normal a la facette
  long calcule_equation_plan_areteFa7(
    long fa7, long iarete,
    double& a, double& b, double& c, double& d) const;

  void creer_sommets_virtuels(const ArrOfInt& liste_sommets,
                              const ArrOfInt& liste_pe,
                              const Schema_Comm_FT& comm);

  void creer_sommets_virtuels_numowner(const ArrOfInt& request_sommets_pe,
                                       const ArrOfInt& request_sommets_num);

  void echanger_sommets_PE(const ArrOfInt& liste_sommets,
                           const ArrOfInt& liste_elem_virtuel_arrivee,
                           const ArrOfInt& liste_face_virtuelle_arrivee,
                           const DoubleTab& deplacement,
                           ArrOfInt& liste_nouveaux_sommets,
                           DoubleTab& deplacement_restant,long skip_facettes=0);

  void corriger_proprietaires_facettes();

  virtual void creer_facettes_virtuelles(const ArrOfInt& liste_facettes,
                                         const ArrOfInt& liste_pe,
                                         const ArrOfInt& facettes_send_pe_list,
                                         const ArrOfInt& facettes_recv_pe_list);

  void buffer_envoyer_facette_PE(long num_facette,
                                 long PE,
                                 long num_element_sur_PE);

  void echanger_facettes(const ArrOfInt& liste_facettes,
                         const ArrOfInt& liste_elem_arrivee,
                         ArrOfInt& facettes_recues_numfacettes,
                         ArrOfInt& facettes_recues_numelement);

  void echanger_facettes_face_x(const ArrOfInt& liste_facettes,
                                const ArrOfInt& liste_face_arrivee,
                                ArrOfInt& facettes_recues_numfacettes,
                                ArrOfInt& facettes_recues_numface);

  void echanger_facettes_face_y(const ArrOfInt& liste_facettes,
                                const ArrOfInt& liste_face_arrivee,
                                ArrOfInt& facettes_recues_numfacettes,
                                ArrOfInt& facettes_recues_numface);

  void echanger_facettes_face_z(const ArrOfInt& liste_facettes,
                                const ArrOfInt& liste_arete_arrivee,
                                ArrOfInt& facettes_recues_numfacettes,
                                ArrOfInt& facettes_recues_numarete);

  void echanger_facettes_arete_x(const ArrOfInt& liste_facettes,
                                 const ArrOfInt& liste_arete_arrivee,
                                 ArrOfInt& facettes_recues_numfacettes,
                                 ArrOfInt& facettes_recues_numarete);

  void echanger_facettes_arete_y(const ArrOfInt& liste_facettes,
                                 const ArrOfInt& liste_arete_arrivee,
                                 ArrOfInt& facettes_recues_numfacettes,
                                 ArrOfInt& facettes_recues_numarete);

  void echanger_facettes_arete_z(const ArrOfInt& liste_facettes,
                                 const ArrOfInt& liste_arete_arrivee,
                                 ArrOfInt& facettes_recues_numfacettes,
                                 ArrOfInt& facettes_recues_numface);

  void convertir_numero_distant_local(const Desc_Structure_FT& descripteur,
                                      const ArrOfInt& element_num_owner,
                                      const ArrOfInt& numeros_distants,
                                      const ArrOfInt& pe_distant,
                                      ArrOfInt& numeros_locaux) const;
  void convertir_numero_distant_local(const Desc_Structure_FT& descripteur,
                                      const ArrOfInt& element_num_owner,
                                      const long numero_distant,
                                      const long pe_distant,
                                      long& numero_local) const;

  virtual void deplacer_sommets(const ArrOfInt& liste_sommets_initiale,
                                const DoubleTab& deplacement_initial,
                                ArrOfInt& liste_sommets_sortis,
                                ArrOfInt& numero_face_sortie, long skip_facettes=0);

  virtual long check_sommets(long error_is_fatal = 1) const;
  virtual long check_mesh(long error_is_fatal = 1, long skip_facette_owner = 0, long skip_facettes = 0) const;

  void construire_noeuds(IntTab& def_noeud,const DoubleTab& soms);
  void calculer_coord_noeuds(const IntTab& def_noeud,const DoubleTab& soms,IntTab& renum);

  REF(Transport_Interfaces_FT_Disc) refequation_transport_;
  // Pour pouvoir utiliser le maillage_FT_IJK sans equation de transport, j'ajoute une ref
  // au domaine_dis, et on l'utilise directement chaque fois que c'est possible au lieu de l'eq. transp.:
  // C'est initialise dans associer_domaine_dis.
  REF(Domaine_dis) refdomaine_dis_;
  // Pour la meme raison, ajout d'une ref au parcours de l'interface:
  REF(Parcours_interface) refparcours_interface_;

  ArrOfIntFT som_init_util_; //Indique les positions de sommets qui sont effectivement
  //la propriete du domaine (ou sous domaine en parallele)

  // Schema de communication permettant des echanges bidirectionnels avec
  // n'importe quel processeur voisin par un joint du maillage eulerien
  Schema_Comm_FT schema_comm_domaine_;

  // Cette variable indique l'etat de l'objet :
  // * RESET :   Maillage vide, non initialise.
  // * MINIMAL : les structures suivantes sont remplies de facon coherente,
  //             et respecte les conventions de maillage sur les differents
  //             processeurs :
  //             - statut_ (le meme pour tous les processeurs)
  //             - sommets_
  //		 - cg_fa7_ // EB
  //             - facettes_
  //             - sommet_elem_
  //		 - sommet_face_ // EB
  //		 - sommet_arete_ // EB
  //             - sommet_face_bord_
  //             - sommet_PE_owner_
  //             - sommet_num_owner_
  //             - desc_sommets_
  //             - desc_facettes_
  //             - drapeaux_sommets_
  //             La liste des facettes peut etre incomplete : les processeurs
  //             pauvres ne connaissent pas encore les facettes dont tous
  //             les sommets sont virtuels
  // * PARCOURU : En plus du maillage minimal, la structure
  //              intersections_elem_facettes_ est remplie et la liste des
  //              facettes est complete (y compris les processeurs pauvres)
  // * COMPLET : Toutes les structures de donnees sont remplies:
  //             - surface_facette_
  //             - normale_facette_
  //             - voisins_
  //             - ...
  Statut_Maillage statut_;

  // Ce compteur est incremente chaque fois que le maillage change
  // (soit les noeuds sont deplaces, soit le nombre de noeuds change, etc...)
  // Cela permet de verifier qu'un tableau de valeurs calculees a partir du maillage
  // est bien a jour (exemple: courbure, ...)
  long mesh_state_tag_;

  // Le temps physique associe a cette interface
  double temps_physique_;

  //////////////////////////////////////////////////////////////////////////////////
  //                MEMBRES DEFINISSANT L'ETAT MINIMAL                            //
  //                                                                              //
  // mes_sommets(i,j)  = j-ieme composante de la position du i-eme sommet du maillage
  // La liste des sommets connus contient au minimum tous les sommets reels et
  // virtuels utilises dans la liste de facettes.
  // Un sommet est reel s'il est situe a l'interieur de ma domaine (c'est a dire
  // si sommet_elem_ >= 0). Un sommet est reel pour exactement un processeur,
  // il est virtuel pour tous les autres. Pour les sommets situes a proximite d'un
  // joint (a epsilon pres), le choix du PE proprietaire est arbitraire.
  DoubleTabFT sommets_;
  DoubleTab sommets2_; // EB est-ce encore utile ???
  // mes_facettes(i,j) = indice du j-ieme sommet de la i-ieme facette du maillage
  //                     En 2D : j=0..1, en 3D j=0..2
  // Conventions :
  //  * le premier sommet de la facette determine le PE proprietaire de la facette:
  //    si ce sommet m'appartient, la facette est reelle, sinon elle est virtuelle.
  //    La facette est dans l'espace virtuel des facettes si et seulement si je ne suis pas
  //    proprietaire du premier sommet. Si elle est reelle, l'espace distant indique
  //    pour quels processeurs elle est virtuelle.
  //  * la normale a la facette pointe vers la phase 1 (l'autre est la phase 0)
  //  * les facettes sont numerotees de facon identique sur tous les processeurs
  //    (ainsi, le premier sommet est toujours le meme).
  //  Une facette est dite "virtuelle pure" si le processeur ne possede aucun des
  //  noeuds de la facette.
  IntTabFT    facettes_;

  // Pour chaque sommet du maillage, numero de l'element eulerien qui
  // le contient et -1 si le sommet est virtuel.
  // Par definition, le sommet m'appartient ssi sommet_elem_ >= 0
  // (si sommet_elem_ < 0, le sommet est virtuel)
  ArrOfIntFT sommet_elem_;
  // EB : on a besoin d'un sommet_face_ pour le calcul de l'indicatrice aux faces
  // Pour chaque sommet du maillage, numero de la face euleriene qui
  // le contient et -1 si le sommet est virtuel.
  // Par definition, le sommet m'appartient ssi sommet_face_ >= 0
  // (si sommet_face_ < 0, le sommet est virtuel)
  IntTabFT sommet_face_;
  // EB : idem pour les aretes
  IntTabFT sommet_arete_;
  // Pour chaque sommet du maillage, numero de la face de bord eulerienne
  // ou se trouve le sommet si le sommet est de type "ligne de contact",
  // -1 si le sommet n'est pas sur le bord. Si le sommet est virtuel et sur
  // un bord, le numero est positif ou nul.
  ArrOfIntFT sommet_face_bord_;
  // Prises ensembles, ces deux donnees forment une cle qui designe un noeud
  // de facon unique.
  // Le numero du PE qui possede le sommet
  ArrOfIntFT sommet_PE_owner_;
  // Le numero du sommet sur le PE proprietaire
  ArrOfIntFT sommet_num_owner_;
  // Le numero de la facette sur le proprietaire de celle-ci
  ArrOfIntFT facette_num_owner_;
  DoubleTabFT cg_fa7_; // EB
  // Descripteur d'elements distants et virtuels pour les sommets et les facettes.
  // Definition : on dit que les espaces distant et virtuels sont coherents pour
  // un element (sommet ou facette) si :
  // * l'element se trouve soit
  //    - dans l'espace virtuel, associe a un unique PE,
  //    - dans l'espace distant, associe a un ou plusieurs PEs,
  //    - ni dans l'espace distant, ni dans l'espace virtuel.
  // * sur le proc B l'element est virtuel, associe au processeur A si et seulement si
  //   cet element est dans l'espace distant du proc A, associe au proc B.
  // La methode Desc_Structure_FT::check() permet de verifier en partie que
  // ces regles sont verifiees.
  Desc_Structure_FT desc_sommets_;
  Desc_Structure_FT desc_facettes_;

  // Attributs des sommets (peut-etre que plusieurs ArrOfBits seraient mieux?)
  // Avec IntTab, echange espace virtuel est plus facile.
  ArrOfIntFT drapeaux_sommets_;
  // Drapeaux attributs des sommets
  static const long EFFACE = 1;
  //                                                                              //
  // FIN DES MEMBRES DEFINISSANT L'ETAT MINIMAL                                   //
  //////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////////
  // LES MEMBRES SUIVANTS PEUVENT ETRE RECONSTRUITS A PARTIR DE L'ETAT MINIMAL    //
  //                                                                              //
  // Intersection entre les facettes du maillage lagrangien et les elements
  // euleriens (ce champ est rempli lors du parcours de l'interface)
  //
  long is_solid_particle_;
  Intersections_Elem_Facettes intersections_elem_facettes_;
  Intersections_Face_Facettes intersections_face_facettes_x_; // EB
  Intersections_Face_Facettes intersections_face_facettes_y_; // EB
  Intersections_Face_Facettes intersections_face_facettes_z_; // EB
  Intersections_Arete_Facettes intersections_arete_facettes_x_; // EB
  Intersections_Arete_Facettes intersections_arete_facettes_y_; // EB
  Intersections_Arete_Facettes intersections_arete_facettes_z_; // EB
  // voisins_(i,j)  = indice de la facette voisine de la i-eme facette du
  // maillage par l'arete j.
  IntTabFT    voisins_;

  /////////////////////////////////////////////////////////////////////////////////
  // Cet objet contient les valeurs calculees a partir du maillage
  // du base (normale, surface, courbure, etc...).
  // Les membres de ce tableau sont calcules lors de l'appel a get_update_...
  // si ce n'est pas deja fait (il se comporte comme un cache).
  // L'appel a un membre public non const du maillage invalide tous les
  // tableaux par increment du "mesh_state_tag".
  // L'appel a get_update_xxx modifie cet objet, bien que la methode
  // get_update_xxx soit "const", mais le principe du "const" est bien
  // respecte: pas d'effet de bord visible de l'exterieur.

  // attention: la methode d'acces est const, et renvoie un objet non const :
  Maillage_FT_Disc_Data_Cache& mesh_data_cache() const;
  DERIV(Maillage_FT_Disc_Data_Cache) mesh_data_cache_;

  // Classes amies (ces classes construisent ou modifient le maillage,
  // les algorithmes dependent fortement des details de mise en oeuvre
  // internes du maillage... il vaut mieux qu'ils soient bons amis...:)
  friend class Marching_Cubes;
  friend class Parcours_interface;
  friend class Remaillage_FT;
  friend class Remaillage_FT_IJK;
  friend class Topologie_Maillage_FT;
  friend class Remailleur_Collision_FT_Collision_Seq;
  friend class Sauvegarde_Reprise_Maillage_FT;

  //parametre de gestion des sorties au format gnuplot
  //-1 : aucunes sorties au format gnuplot
  // 1 : 1 sortie par pas de temps (dans mettre_a_jour)
  long niveau_plot_;
  double correction_contact_courbure_coeff_;
  long calcul_courbure_iterations_;
  long niter_pre_lissage_;
  enum enum_methode_calcul_courbure_contact_line_ { STANDARD=0, MIRROR=1, IMPROVED=2, none=3, WEIGHTED=4, HYSTERESIS=5 };
  long methode_calcul_courbure_contact_line_;
  double weight_CL_; // Le poids de la CL


};

// Description :
//  Renvoie 0 si le sommet m'appartient, 1 sinon.
//  Si le sommet i m'appartient, alors il se trouve a l'interieur de l'element
//  sommet_elem_[i]. On garantit qu'un sommet est reel sur exactement un
//  processeur et virtuel sur tous les autres.
inline long Maillage_FT_Disc::sommet_virtuel(long i) const
{
  return (sommet_elem_[i] < 0) ? 1 : 0;
}
// debut EB
// Description :
//  Renvoie 0 si le sommet m'appartient, 1 sinon.
//  Si le sommet i m'appartient, alors il se trouve a l'interieur de la face
//  sommet_face_[i]. On garantit qu'un sommet est reel sur exactement un
//  processeur et virtuel sur tous les autres.
inline long Maillage_FT_Disc::sommet_virtuel_face(long i) const { return (sommet_face_[i] < 0) ? 1 : 0 ; }
// fin EB

// Description :
//  Renvoie 1 si le sommet se trouve sur un bord du domaine, 0 sinon.
//  Si le sommet est virtuel mais qu'il se trouve sur un bord du domaine
//  du proceseur voisin, on renvoie 1 aussi.
inline long Maillage_FT_Disc::sommet_ligne_contact(long i) const
{
  return (sommet_face_bord_[i] >= 0) ? 1 : 0;
}

// Description :
//  Renvoie le numero de la face de bord.
inline long Maillage_FT_Disc::sommet_face_bord(long i) const
{
  return sommet_face_bord_[i];
}



/*! @brief Renvoie 0 si la facette m'appartient, 1 sinon.
 *
 * (le test est "le premier sommet de la facette m'appartient-t-il ?")
 *
 */
inline long Maillage_FT_Disc::facette_virtuelle(long i) const
{
  const long sommet = facettes_(i, 0);
  return (sommet_elem_[sommet] < 0) ? 1 : 0;
}
// debut EB
// Description:
//  Renvoie 0 si la facette m'appartient, 1 sinon.
//  (le test est "le premier sommet de la facette m'appartient-t-il ?")
inline long Maillage_FT_Disc::facette_virtuelle_face(long i) const
{
  const long sommet = facettes_(i, 0);
  return (sommet_face_[sommet] < 0) ? 1 : 0;
}
// fin EB
inline long Maillage_FT_Disc::set_niveau_plot(long niv)
{
  niveau_plot_ = niv;
  return niveau_plot_;
}


class Maillage_FT_Disc_Data_Cache : public Objet_U
{
  Declare_instanciable(Maillage_FT_Disc_Data_Cache);
public:
  void clear();

  long        tag_surface_;
  ArrOfDoubleFT surface_facettes_;
  long        tag_normale_;
  DoubleTabFT   normale_facettes_;
  long        tag_courbure_;
  ArrOfDoubleFT courbure_sommets_;
};
#endif