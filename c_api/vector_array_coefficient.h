#ifndef CMFEM_VECTOR_ARRAY_COEFFICIENT_H
#define CMFEM_VECTOR_ARRAY_COEFFICIENT_H

#include "common.h"

CMFEM_STORAGE(VectorArrayCoefficient, 88);

CMFEM_BEGIN_EXTERN_C

CMFEM_VectorArrayCoefficient *CMFEM_VectorArrayCoefficient_New(int dim);
void CMFEM_VectorArrayCoefficient_Delete(CMFEM_VectorArrayCoefficient
                                         *coefficient);
void CMFEM_VectorArrayCoefficient_SetConstantCoefficient(
   CMFEM_VectorArrayCoefficient *coefficient, int component,
   CMFEM_ConstantCoefficient *entry);
void CMFEM_VectorArrayCoefficient_SetPWConstCoefficient(
   CMFEM_VectorArrayCoefficient *coefficient, int component,
   CMFEM_PWConstCoefficient *entry);

CMFEM_END_EXTERN_C

#endif
