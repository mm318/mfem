#ifndef CMFEM_FUNCTION_COEFFICIENT_H
#define CMFEM_FUNCTION_COEFFICIENT_H

#include "common.h"

CMFEM_STORAGE(FunctionCoefficient, 112);

typedef double (*CMFEM_FunctionCallback)(const CMFEM_Vector *x, void *context);

CMFEM_BEGIN_EXTERN_C

CMFEM_FunctionCoefficient *CMFEM_FunctionCoefficient_New(
   CMFEM_FunctionCallback callback,
   void *context);
void CMFEM_FunctionCoefficient_Delete(CMFEM_FunctionCoefficient *coefficient);

CMFEM_END_EXTERN_C

#endif
