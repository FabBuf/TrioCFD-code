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

#include <Op_Dift_VDF_Face_leaves.h>

Implemente_instanciable_sans_constructeur(Op_Dift_VDF_Face,"Op_Dift_VDF_Face",Op_Dift_VDF_Face_base);
Sortie& Op_Dift_VDF_Face::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Op_Dift_VDF_Face::readOn(Entree& s ) { return s ; }
Op_Dift_VDF_Face::Op_Dift_VDF_Face() : Op_Dift_VDF_Face_base(Iterateur_VDF_Face<Eval_Dift_VDF_Face>())
{
  declare_support_masse_volumique(1);
}
Implemente_instanciable_sans_constructeur(Op_Dift_VDF_Face_FT,"Op_Dift_VDF_Face_FT",Op_Dift_VDF_Face_base);
Sortie& Op_Dift_VDF_Face_FT::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Op_Dift_VDF_Face_FT::readOn(Entree& s ) { return s ; }
Op_Dift_VDF_Face_FT::Op_Dift_VDF_Face_FT() : Op_Dift_VDF_Face_base(Iterateur_VDF_Face<Eval_Dift_VDF_Face_FT>())
{
  declare_support_masse_volumique(1);
}

Implemente_instanciable(Op_Dift_VDF_Face_Axi,"Op_Dift_VDF_Face_Axi",Op_Dift_VDF_Face_Axi_base);
Sortie& Op_Dift_VDF_Face_Axi::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Op_Dift_VDF_Face_Axi::readOn(Entree& s ) { return s ; }
