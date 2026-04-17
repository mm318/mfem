#ifndef CMFEM_GRID_FUNCTION_H
#define CMFEM_GRID_FUNCTION_H

#include "common.h"

CMFEM_STORAGE(GridFunction, 104);

CMFEM_BEGIN_EXTERN_C

CMFEM_GridFunction *CMFEM_GridFunction_New(CMFEM_FiniteElementSpace *fespace);
void CMFEM_GridFunction_Delete(CMFEM_GridFunction *grid_function);
void CMFEM_GridFunction_Assign(CMFEM_GridFunction *grid_function, double value);
void CMFEM_GridFunction_ProjectVectorFunctionCoefficient(CMFEM_GridFunction *grid_function, const CMFEM_VectorFunctionCoefficient *coefficient);
double CMFEM_GridFunction_ComputeL2ErrorVectorFunctionCoefficient(const CMFEM_GridFunction *grid_function, const CMFEM_VectorFunctionCoefficient *coefficient);
void CMFEM_GridFunction_Add(CMFEM_GridFunction *grid_function, const CMFEM_GridFunction *other);
void CMFEM_GridFunction_Scale(CMFEM_GridFunction *grid_function, double scale);
void CMFEM_GridFunction_Save(const CMFEM_GridFunction *grid_function, const char *path, int precision);

CMFEM_END_EXTERN_C

#endif
