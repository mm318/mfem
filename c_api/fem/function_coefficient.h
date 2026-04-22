#ifndef CMFEM_FUNCTION_COEFFICIENT_H
#define CMFEM_FUNCTION_COEFFICIENT_H

#include "../common.h"

CMFEM_STORAGE(FunctionCoefficient, 112);

typedef double (*CMFEM_FunctionCallback)(const CMFEM_Vector *x, void *context);
typedef double (*CMFEM_TimeDependentFunctionCallback)(const CMFEM_Vector *x,
                                                      double time,
                                                      void *context);

CMFEM_BEGIN_EXTERN_C

CMFEM_FunctionCoefficient *CMFEM_FunctionCoefficient_New(
   CMFEM_FunctionCallback callback,
   void *context);
CMFEM_FunctionCoefficient *CMFEM_FunctionCoefficient_NewTimeDependent(
   CMFEM_TimeDependentFunctionCallback callback,
   void *context);
void CMFEM_FunctionCoefficient_SetTime(CMFEM_FunctionCoefficient *coefficient,
                                       double time);
void CMFEM_FunctionCoefficient_Delete(CMFEM_FunctionCoefficient *coefficient);

CMFEM_END_EXTERN_C

#endif
