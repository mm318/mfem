#ifndef CMFEM_GRID_FUNCTION_H
#define CMFEM_GRID_FUNCTION_H

#include "common.h"

CMFEM_STORAGE(GridFunction, 104);

CMFEM_BEGIN_EXTERN_C

CMFEM_GridFunction *CMFEM_GridFunction_New(CMFEM_FiniteElementSpace *fespace);
void CMFEM_GridFunction_Delete(CMFEM_GridFunction *grid_function);
void CMFEM_GridFunction_Assign(CMFEM_GridFunction *grid_function, double value);
double CMFEM_GridFunction_Get(const CMFEM_GridFunction *grid_function,
                              int index);
void CMFEM_GridFunction_Set(CMFEM_GridFunction *grid_function,
                            int index,
                            double value);
int CMFEM_GridFunction_Size(const CMFEM_GridFunction *grid_function);
void CMFEM_GridFunction_ProjectCoefficientFc(
   CMFEM_GridFunction *grid_function,
   const CMFEM_FunctionCoefficient *coefficient);
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
void CMFEM_GridFunction_ComputeFluxDi(
   const CMFEM_GridFunction *grid_function,
   CMFEM_DiffusionIntegrator *integrator,
   CMFEM_GridFunction *flux);
void CMFEM_GridFunction_SetTrueVector(CMFEM_GridFunction *grid_function);
void CMFEM_GridFunction_SetFromTrueVector(CMFEM_GridFunction *grid_function);
double CMFEM_GridFunction_ComputeL2ErrorFc(
   const CMFEM_GridFunction *grid_function,
   const CMFEM_FunctionCoefficient *coefficient);
double CMFEM_GridFunction_ComputeL2ErrorFcOrder(
   const CMFEM_GridFunction *grid_function,
   const CMFEM_FunctionCoefficient *coefficient,
   int quadrature_order);
double CMFEM_GridFunction_ComputeL2ErrorVfc(
   const CMFEM_GridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *coefficient);
double CMFEM_GridFunction_ComputeHCurlErrorVfcVfc(
   const CMFEM_GridFunction *grid_function,
   CMFEM_VectorFunctionCoefficient *exact_solution,
   CMFEM_VectorFunctionCoefficient *exact_curl);
void CMFEM_GridFunction_CopyToVec(const CMFEM_GridFunction *grid_function,
                                  CMFEM_Vector *vector);
void CMFEM_GridFunction_SetFromVec(CMFEM_GridFunction *grid_function,
                                   const CMFEM_Vector *vector);
void CMFEM_GridFunction_Add(CMFEM_GridFunction *grid_function,
                            const CMFEM_GridFunction *other);
void CMFEM_GridFunction_Scale(CMFEM_GridFunction *grid_function, double scale);
void CMFEM_GridFunction_GetTrueDofs(const CMFEM_GridFunction *grid_function,
                                    CMFEM_Vector *true_dofs);
void CMFEM_GridFunction_SetFromTrueDofs(CMFEM_GridFunction *grid_function,
                                        const CMFEM_Vector *true_dofs);
void CMFEM_GridFunction_Update(CMFEM_GridFunction *grid_function);
void CMFEM_GridFunction_Save(const CMFEM_GridFunction *grid_function,
                             const char *path, int precision);

CMFEM_END_EXTERN_C

#endif
