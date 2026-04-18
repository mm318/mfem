#ifndef CMFEM_GRID_FUNCTION_COEFFICIENT_H
#define CMFEM_GRID_FUNCTION_COEFFICIENT_H

#include "common.h"

CMFEM_STORAGE(GridFunctionCoefficient, 32);

CMFEM_BEGIN_EXTERN_C

CMFEM_GridFunctionCoefficient *CMFEM_GridFunctionCoefficient_NewGf(
   const CMFEM_GridFunction *grid_function);
void CMFEM_GridFunctionCoefficient_Delete(
   CMFEM_GridFunctionCoefficient *coefficient);

CMFEM_END_EXTERN_C

#endif
