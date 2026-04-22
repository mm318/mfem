#ifndef CMFEM_VECTOR_FUNCTION_COEFFICIENT_H
#define CMFEM_VECTOR_FUNCTION_COEFFICIENT_H

#include "../common.h"

CMFEM_STORAGE(VectorFunctionCoefficient, 144);

typedef void (*CMFEM_VectorFunctionCallback)(const CMFEM_Vector *x,
                                             CMFEM_Vector *value,
                                             void *context);

CMFEM_BEGIN_EXTERN_C

CMFEM_VectorFunctionCoefficient *CMFEM_VectorFunctionCoefficient_New(int dim,
                                                                     CMFEM_VectorFunctionCallback callback, void *context);
void CMFEM_VectorFunctionCoefficient_Delete(CMFEM_VectorFunctionCoefficient
                                            *coefficient);

CMFEM_END_EXTERN_C

#endif
