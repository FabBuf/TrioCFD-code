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
// File:        Ch_front_input.h
// Directory:   $TRUST_ROOT/src/Kernel/Champs
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Ch_front_input_included
#define Ch_front_input_included


#include <Ch_front_var_instationnaire_dep.h>
#include <Champ_Input_Proto.h>
#include <IntTab.h>

//////////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//     class Ch_front_input
//
//     Cette classe represente un champ accessible par setInputField
//     defini sur une frontiere avec une valeur par face.
//
// .SECTION voir aussi
//   Champ_Input_Proto
/////////////////////////////////////////////////////////////////////////////////

class Ch_front_input : public Ch_front_var_instationnaire_dep, public Champ_Input_Proto
{
  Declare_instanciable(Ch_front_input);

public:

  virtual Champ_front_base& affecter_(const Champ_front_base&)
  {
    return *this;
  }
  virtual void getTemplate(TrioField& afield) const;
  virtual void setValue(const TrioField& afield);

  virtual int initialiser(double temps, const Champ_Inc_base& inco);
protected:

  // Factorisation function between several input field classes
  virtual void set_nb_comp(int i); // calls fixer_nb_comp
  virtual void set_name(const Nom& ); // calls nommer
  virtual const Nom& get_name() const; // calls le_nom

  void buildSommetsFaces() const; // const because used in Ch_Front_input_ALE::getTemplate() which is const - actually updates the 2 members below due to ALE mesh movement:
  mutable DoubleTab sommets_;
  mutable IntTab faces_;


};

#endif

