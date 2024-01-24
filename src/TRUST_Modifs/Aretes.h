/****************************************************************************
* Copyright (c) 2022, CEA
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

#ifndef Aretes_included
#define Aretes_included

#include <TRUSTTab.h>

class Domaine_VDF;

//////////////////////////////////////////////////////////////////////////////
//
// Classe Aretes : represente un ensemble d' aretes
//
//////////////////////////////////////////////////////////////////////////////

class Aretes: public Objet_U
{
  Declare_instanciable(Aretes);
public:

  inline Aretes(int n);
  void affecter(int&, int, int, int, int, int, int, int, const ArrOfInt&);
  int affecter_aretes(int&, int, int, int,int, int, int, int, const ArrOfInt&); // EB
  int affecter_aretes_virtuelle(int&, int, int, int, int, int, int, int); // EB
  inline IntTab& faces();
  inline IntVect& type1(); // EB
  inline IntVect& type2(); // EB
  void dimensionner(int);
  void trier(int&, int&, int&, int&);
  void trier(int&, int&, int&, int&, const int nb_aretes_reelles,const int nb_elem_reels,IntTab& Aretes_Som, IntTab& Elem_Aretes); // EB
  void trier_pour_debog(int&, int&, int&, int&, const DoubleTab&);
  void calculer_centre_de_gravite(Domaine_VDF& domaine);

private:
  IntTab faces_;
  IntVect type1_, type2_;
  void swap(int, int);
  void swap_Aretes_Som(int a1, int a2, IntTab& Aretes_Som); // EB
  void swap_Elem_Arete(int a1, int a2, IntTab& Elem_Aretes); // EB
};

/*! @brief Constructeur : dimensionne les tableaux
 *
 */
inline Aretes::Aretes(int n) : faces_(n,4),type1_(n),type2_(n)
{
}

/*! @brief retourne le tableau (IntTab) des faces de dimension nb_aretes,4
 *
 *  faces(num_arete,i) 0 <= i <= 3:
 *  numeros des 4 faces qui se partagent l'arete num_arete
 *  Pour une arete interne les 4 faces existent
 *  Pour une arete bord un ou deux des numeros valent -1.
 *
 */
inline IntTab& Aretes::faces() { return faces_; }

// EB orientation des aretes
/*! @brief retourne l'orientation des aretes selon la definition dans Aretes::affecter_aretes
 *
 */
inline IntVect& Aretes::type1() { return type1_; }

// EB type des aretes
/*! @brief retourne le type des aretes selon la definition dans Aretes::affecter_aretes
 *
 */
inline IntVect& Aretes::type2() { return type2_; }

#endif  // Aretes_H
