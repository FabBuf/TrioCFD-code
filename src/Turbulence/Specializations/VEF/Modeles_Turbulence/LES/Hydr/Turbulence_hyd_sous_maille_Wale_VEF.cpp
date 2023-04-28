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
// File:        Turbulence_hyd_sous_maille_Wale_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Turbulence_hyd_sous_maille_Wale_VEF.h>
#include <Champ_P1NC.h>
#include <TRUSTTrav.h>
#include <Debog.h>
#include <Schema_Temps_base.h>
#include <Param.h>
#include <Equation_base.h>
#include <Domaine_VEF.h>

Implemente_instanciable_sans_constructeur(Turbulence_hyd_sous_maille_Wale_VEF,"Modele_turbulence_hyd_sous_maille_Wale_VEF",Mod_turb_hyd_ss_maille_VEF);

Turbulence_hyd_sous_maille_Wale_VEF::Turbulence_hyd_sous_maille_Wale_VEF()
{
  declare_support_masse_volumique(1);
  cw=0.5;
}
//// printOn
//

Sortie& Turbulence_hyd_sous_maille_Wale_VEF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Turbulence_hyd_sous_maille_Wale_VEF::readOn(Entree& is )
{
  Mod_turb_hyd_ss_maille_VEF::readOn(is);
  return is;
}

void Turbulence_hyd_sous_maille_Wale_VEF::set_param(Param& param)
{
  Mod_turb_hyd_ss_maille_VEF::set_param(param);
  param.ajouter("cw",&cw);
  param.ajouter_condition("value_of_cw_ge_0","sous_maille_Wale model constant must be positive.");
}

///////////////////////////////////////////////////////////////////////////////
//
//  Implementation de fonctions de la classe Turbulence_hyd_sous_maille_Wale_VEF
//
//////////////////////////////////////////////////////////////////////////////


Champ_Fonc& Turbulence_hyd_sous_maille_Wale_VEF::calculer_viscosite_turbulente()
{
  // cw est la constante du modele WALE qui correspond a une correction
  //  de la constante Cs du modele de Smagorinsky.
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  double temps = mon_equation->inconnue().temps();
  DoubleTab& visco_turb = la_viscosite_turbulente.valeurs();
  const int nb_elem = domaine_VEF.nb_elem();
  const DoubleTab& la_vitesse = mon_equation->inconnue().valeurs();
  const Domaine_Cl_VEF& domaine_Cl_VEF = le_dom_Cl_VEF.valeur();

  if (visco_turb.size() != nb_elem)
    {
      Cerr << "Size error for the array containing the values of the turbulent viscosity." << finl;
      exit();
    }
  Debog::verifier("Turbulence_hyd_sous_maille_Wale_VEF::calculer_viscosite_turbulente visco_turb 0",visco_turb);

  // Patrick : on travaille sur le champ filtre.
  //const Champ_P1NC& ch=(const Champ_P1NC&) mon_equation->inconnue().valeur();
  DoubleTab ubar(la_vitesse);
  //  ch.filtrer_L2(ubar);
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  DoubleTab duidxj(nb_elem_tot,dimension,dimension);
  Champ_P1NC::calcul_gradient(ubar,duidxj,domaine_Cl_VEF);

  // OpenMP ToDo:
  // l_ passe une seule fois dans une region
  // Utiliser copyToDevice pour duidxj_addr et l_addr
  double* visco_turb_addr=visco_turb.addr();
  const double* l_addr=l_.addr();
  const double* duidxj_addr=duidxj.addr();
  #pragma omp target teams map(to:duidxj_addr[0:duidxj.size_array()],l_addr[0:l_.size_array()]) map(from:visco_turb_addr[0:visco_turb.size_array()])
  {
    double gij2[3][3] {};
    double sd[3][3] {};
    #pragma omp distribute parallel for private(gij2, sd)
    for (int elem=0; elem<nb_elem; elem ++)
      {

        // Calcul du terme gij2.
        for(int i=0; i<dimension; i++)
          for(int j=0; j<dimension; j++)
            {
              gij2[i][j]=0;
              for(int k=0; k<dimension; k++)
                gij2[i][j]+=duidxj_addr[elem*dimension*dimension+i*dimension+k] * duidxj_addr[elem*dimension*dimension+k*dimension+j];
            }

        // Calcul du terme gkk2.
        double gkk2=0;
        for(int k=0; k<dimension; k++)
          gkk2 += gij2[k][k];

        // Calcul de sd.
        for(int i=0; i<dimension; i++)
          for(int j=0; j<dimension; j++)
            {
              sd[i][j]=0.5*(gij2[i][j]+gij2[j][i]);
              if(i==j)
                {
                  sd[i][j]-=gkk2/3.; // Terme derriere le tenseur de Kronecker.
                }
            }

        // Calcul de sd2 et Sij2.
        double sd2=0.;
        double Sij2=0.;
        for (int i=0 ; i<dimension ; i++)
          for (int j=0 ; j<dimension ; j++)
            {
              sd2+=sd[i][j]*sd[i][j];
              //Deplacement du calcul de Sij
              double Sij=0.5*(duidxj_addr[elem*dimension*dimension+i*dimension+j] + duidxj_addr[elem*dimension*dimension+j*dimension+i]);
              Sij2+=Sij*Sij;
            }

        // Calcul de OP1 et OP2.
        // Replace pow by sqrt and multiply, faster
        //OP1=pow(sd2,1.5);
        double OP1=sd2*sqrt(sd2);
        //OP2=pow(Sij2,2.5)+pow(sd2,1.25);
        double OP2=Sij2*Sij2*sqrt(Sij2)+sd2*sqrt(sqrt(sd2));

        if (OP1!=0.) // donc sd2 et OP2 par voie de consequence sont differents de zero
          visco_turb_addr[elem]=cw*cw*l_addr[elem]*l_addr[elem]*OP1/OP2;
        else
          visco_turb_addr[elem]=0;
      }// fin de la boucle sur les elements
  }

  Debog::verifier("Turbulence_hyd_sous_maille_Wale_VEF::calculer_viscosite_turbulente visco_turb 1",visco_turb);

  la_viscosite_turbulente.changer_temps(temps);
  return la_viscosite_turbulente;
}
