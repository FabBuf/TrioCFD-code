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
// File:        Source_Gravite_PF_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Phase_field/src/VDF
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Gravite_PF_VDF.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Les_Cl.h>
#include <Zone_Cl_dis.h>
#include <Navier_Stokes_phase_field.h>

Implemente_instanciable(Source_Gravite_PF_VDF,"Source_Gravite_PF_VDF",Source_base);


// Description:
//    Imprime la source sur un flot de sortie.
// Precondition:
// Parametre: Sortie& os
//    Signification: le flot de sortie pour l'impression
//    Valeurs par defaut:
//    Contraintes:
//    Acces: sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord: le flot de sortie est modifie
// Postcondition: la methode ne modifie pas l'objet
Sortie& Source_Gravite_PF_VDF::printOn(Sortie& os) const
{
  os <<que_suis_je()<< finl;
  return os;
}

// Description:
//    Lecture de la source sur un flot d'entree.
// Precondition:
// Parametre: Entree& is
//    Signification: le flot d'entree pour la lecture des parametres
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& Source_Gravite_PF_VDF::readOn(Entree& is)
{
  return is;
}

// Description:
//    Remplit le tableau volumes
// Precondition:
// Parametre: Entree& is
//    Signification: le flot d'entree pour la lecture des parametres
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Source_Gravite_PF_VDF::associer_zones(const Zone_dis& zone,const Zone_Cl_dis& zone_cl)
{
  la_zone = ref_cast(Zone_VDF,zone.valeur());
  la_zone_Cl = ref_cast(Zone_Cl_VDF,zone_cl.valeur());
}



// Description:
//    Ajoute les termes sources
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Entree&
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
DoubleTab& Source_Gravite_PF_VDF::ajouter(DoubleTab& resu) const
{
  int face, nb_faces = la_zone->nb_faces();
  int premiere_face_interne = la_zone->premiere_face_int();

  const IntVect& orientation = la_zone->orientation();
  const DoubleVect& volumes_entrelaces = la_zone->volumes_entrelaces();
  DoubleVect porosite_surf ;          // porosites surfaciques
  porosite_surf.ref(la_zone->porosite_face());
  const Navier_Stokes_phase_field& eq_NS_PF = ref_cast(Navier_Stokes_phase_field,equation());
  const DoubleVect& g = eq_NS_PF.get_g_();

  int num_cl;
  for (num_cl=0 ; num_cl<la_zone->nb_front_Cl() ; num_cl++)
    {
      const Cond_lim& la_cl = la_zone_Cl->les_conditions_limites(num_cl);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();
      if (sub_type(Dirichlet,la_cl.valeur()) || sub_type(Dirichlet_homogene,la_cl.valeur()))
        ;
      else
        for (face=ndeb ; face<nfin ; face++)
          {
            resu(face) += g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);

          }
    }

  for (face=premiere_face_interne ; face<nb_faces; face++)
    {
      resu(face) += g(orientation(face)) * porosite_surf(face) * volumes_entrelaces(face);
    }


  return resu;
}


// Description:
//    Calcule la contribution de cette source
// Precondition:
// Parametre: DoubleTab& resu
//    Signification: flux
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: DoubleTab&
//    Signification: le flux
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
DoubleTab& Source_Gravite_PF_VDF::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}
