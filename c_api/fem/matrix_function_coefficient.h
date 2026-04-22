#ifndef CMFEM_MATRIX_FUNCTION_COEFFICIENT_H
#define CMFEM_MATRIX_FUNCTION_COEFFICIENT_H

#include "../common.h"

CMFEM_STORAGE(MatrixFunctionCoefficient, 240);

typedef void (*CMFEM_MatrixFunctionCallback)(const CMFEM_Vector *x,
                                             double *values,
                                             int dim,
                                             void *context);

CMFEM_BEGIN_EXTERN_C

CMFEM_MatrixFunctionCoefficient *CMFEM_MatrixFunctionCoefficient_New(
   int dim,
   CMFEM_MatrixFunctionCallback callback,
   void *context);
void CMFEM_MatrixFunctionCoefficient_Delete(
   CMFEM_MatrixFunctionCoefficient *coefficient);

CMFEM_END_EXTERN_C

#endif
