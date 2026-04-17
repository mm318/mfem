#ifndef CMFEM_GRID_FUNCTION_H
#define CMFEM_GRID_FUNCTION_H

#include "common.h"

CMFEM_STORAGE(GridFunction, 104);

CMFEM_BEGIN_EXTERN_C

CMFEM_GridFunction *CMFEM_GridFunction_New(CMFEM_FiniteElementSpace *fespace);
void CMFEM_GridFunction_Delete(CMFEM_GridFunction *grid_function);
void CMFEM_GridFunction_Assign(CMFEM_GridFunction *grid_function, double value);
void CMFEM_GridFunction_ProjectBdrCoefficientCc(
   CMFEM_GridFunction *grid_function,
   const CMFEM_ConstantCoefficient *coefficient,
   const CMFEM_ArrayInt *attributes);
void CMFEM_GridFunction_ProjectBdrCoefficientFc(
   CMFEM_GridFunction *grid_function,
   const CMFEM_FunctionCoefficient *coefficient,
   const CMFEM_ArrayInt *attributes);
void CMFEM_GridFunction_ProjectBdrCoefficientVcc(
   CMFEM_GridFunction *grid_function,
   const CMFEM_VectorConstantCoefficient *coefficient,
   const CMFEM_ArrayInt *attributes);
void CMFEM_GridFunction_ProjectCoefficientVfc(
   CMFEM_GridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *coefficient);
double CMFEM_GridFunction_ComputeL2ErrorFc(
   const CMFEM_GridFunction *grid_function,
   const CMFEM_FunctionCoefficient *coefficient);
double CMFEM_GridFunction_ComputeL2ErrorVfc(
   const CMFEM_GridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *coefficient);
void CMFEM_GridFunction_Add(CMFEM_GridFunction *grid_function,
                            const CMFEM_GridFunction *other);
void CMFEM_GridFunction_Scale(CMFEM_GridFunction *grid_function, double scale);
void CMFEM_GridFunction_Update(CMFEM_GridFunction *grid_function);
void CMFEM_GridFunction_Save(const CMFEM_GridFunction *grid_function,
                             const char *path, int precision);

CMFEM_END_EXTERN_C

#endif
