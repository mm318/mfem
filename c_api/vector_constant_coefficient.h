#ifndef CMFEM_VECTOR_CONSTANT_COEFFICIENT_H
#define CMFEM_VECTOR_CONSTANT_COEFFICIENT_H

#include "common.h"

CMFEM_STORAGE(VectorConstantCoefficient, 64);

CMFEM_BEGIN_EXTERN_C

CMFEM_VectorConstantCoefficient *CMFEM_VectorConstantCoefficient_New(
   const CMFEM_Vector *value);
void CMFEM_VectorConstantCoefficient_Delete(
   CMFEM_VectorConstantCoefficient *coefficient);

CMFEM_END_EXTERN_C

#endif
