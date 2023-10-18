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
/////////////////////////////////////////////////////////////////////////////
//
// File      : IJK_Lata_writer.h
// Directory : $IJK_ROOT/src/IJK
//
/////////////////////////////////////////////////////////////////////////////
#include <IJK_Field.h>
//
// WARNING: DO NOT EDIT THIS FILE! Only edit the template file IJK_Lata_writer.h.P
//
// #include <FixedVector.h>


void dumplata_header(const char *filename);
void dumplata_header(const char *filename, const IJK_Field_double& f);
void dumplata_add_geometry(const char *filename, const IJK_Field_double& f);

void dumplata_newtime(const char *filename, double time);

void dumplata_vector(const char *filename, const char *fieldname,
                     const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                     int step);

void dumplata_vector_parallele_plan(const char *filename, const char *fieldname,
                                    const IJK_Field_double& vx, const IJK_Field_double& vy, const IJK_Field_double& vz,
                                    int step);

void dumplata_cellvector(const char *filename, const char *fieldname,
                         const FixedVector<IJK_Field_double, 3>& v,
                         int step);

void dumplata_scalar(const char *filename, const char *fieldname,
                     const IJK_Field_double& f,
                     int step);

void dumplata_scalar_parallele_plan(const char *filename, const char *fieldname,
                                    const IJK_Field_double& f,
                                    int step);

void lire_dans_lata(const char *filename, int tstep, const char *geometryname, const char *fieldname,
                    IJK_Field_double& f);

void lire_dans_lata(const char *filename, int tstep, const char *geometryname, const char *fieldname,
                    IJK_Field_double& vx, IJK_Field_double& vy, IJK_Field_double& vz);


void dumplata_header(const char *filename);
void dumplata_header(const char *filename, const IJK_Field_float& f);
void dumplata_add_geometry(const char *filename, const IJK_Field_float& f);

void dumplata_newtime(const char *filename, double time);

void dumplata_vector(const char *filename, const char *fieldname,
                     const IJK_Field_float& vx, const IJK_Field_float& vy, const IJK_Field_float& vz,
                     int step);

void dumplata_vector_parallele_plan(const char *filename, const char *fieldname,
                                    const IJK_Field_float& vx, const IJK_Field_float& vy, const IJK_Field_float& vz,
                                    int step);

void dumplata_cellvector(const char *filename, const char *fieldname,
                         const FixedVector<IJK_Field_float, 3>& v,
                         int step);

void dumplata_scalar(const char *filename, const char *fieldname,
                     const IJK_Field_float& f,
                     int step);

void dumplata_scalar_parallele_plan(const char *filename, const char *fieldname,
                                    const IJK_Field_float& f,
                                    int step);

void lire_dans_lata(const char *filename, int tstep, const char *geometryname, const char *fieldname,
                    IJK_Field_float& f);

void lire_dans_lata(const char *filename, int tstep, const char *geometryname, const char *fieldname,
                    IJK_Field_float& vx, IJK_Field_float& vy, IJK_Field_float& vz);


#if 0
// Forward declaration de la classe Maillage_FT_IJK
// (permet d'utiliser des Maillage_FT_IJK &)
class Maillage_FT_IJK;
void dumplata_ft_mesh(const char *filename, const char *meshname,
                      const Maillage_FT_IJK& mesh, int step);
#endif
void dumplata_ft_field(const char *filename, const char *meshname,
                       const char *field_name, const char *localisation,
                       const ArrOfInt& field, int step);
void dumplata_ft_field(const char *filename, const char *meshname,
                       const char *field_name, const char *localisation,
                       const ArrOfDouble& field, int step);