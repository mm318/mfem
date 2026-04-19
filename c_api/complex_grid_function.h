#ifndef CMFEM_COMPLEX_GRID_FUNCTION_H
#define CMFEM_COMPLEX_GRID_FUNCTION_H

#include "common.h"

CMFEM_STORAGE(ComplexGridFunction, 80);

CMFEM_BEGIN_EXTERN_C

CMFEM_ComplexGridFunction *CMFEM_ComplexGridFunction_New(
   CMFEM_FiniteElementSpace *fespace);
void CMFEM_ComplexGridFunction_Delete(CMFEM_ComplexGridFunction *grid_function);
void CMFEM_ComplexGridFunction_Assign(CMFEM_ComplexGridFunction *grid_function,
                                      double real_value,
                                      double imag_value);
void CMFEM_ComplexGridFunction_ProjectCoefficientFcFc(
   CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_FunctionCoefficient *real_coefficient,
   const CMFEM_FunctionCoefficient *imag_coefficient);
void CMFEM_ComplexGridFunction_ProjectCoefficientVfcVfc(
   CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *real_coefficient,
   const CMFEM_VectorFunctionCoefficient *imag_coefficient);
void CMFEM_ComplexGridFunction_ProjectBdrCoefficientFcFcAi(
   CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_FunctionCoefficient *real_coefficient,
   const CMFEM_FunctionCoefficient *imag_coefficient,
   const CMFEM_ArrayInt *attributes);
void CMFEM_ComplexGridFunction_ProjectBdrCoefficientTangentVfcVfcAi(
   CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *real_coefficient,
   const CMFEM_VectorFunctionCoefficient *imag_coefficient,
   const CMFEM_ArrayInt *attributes);
void CMFEM_ComplexGridFunction_ProjectBdrCoefficientNormalVfcVfcAi(
   CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *real_coefficient,
   const CMFEM_VectorFunctionCoefficient *imag_coefficient,
   const CMFEM_ArrayInt *attributes);
double CMFEM_ComplexGridFunction_ComputeL2ErrorRealFc(
   const CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_FunctionCoefficient *coefficient);
double CMFEM_ComplexGridFunction_ComputeL2ErrorImagFc(
   const CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_FunctionCoefficient *coefficient);
double CMFEM_ComplexGridFunction_ComputeL2ErrorRealVfc(
   const CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *coefficient);
double CMFEM_ComplexGridFunction_ComputeL2ErrorImagVfc(
   const CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *coefficient);
double CMFEM_ComplexGridFunction_ComputeL2ErrorRealVfcAiQo(
   const CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *coefficient,
   const CMFEM_ArrayInt *element_marker,
   int quadrature_order);
double CMFEM_ComplexGridFunction_ComputeL2ErrorImagVfcAiQo(
   const CMFEM_ComplexGridFunction *grid_function,
   const CMFEM_VectorFunctionCoefficient *coefficient,
   const CMFEM_ArrayInt *element_marker,
   int quadrature_order);
void CMFEM_ComplexGridFunction_CopyRealToGf(
   const CMFEM_ComplexGridFunction *grid_function,
   CMFEM_GridFunction *real_grid_function);
void CMFEM_ComplexGridFunction_CopyImagToGf(
   const CMFEM_ComplexGridFunction *grid_function,
   CMFEM_GridFunction *imag_grid_function);
void CMFEM_ComplexGridFunction_Save(const CMFEM_ComplexGridFunction
                                    *grid_function,
                                    const char *path,
                                    int precision);
void CMFEM_ComplexGridFunction_SaveReal(
   const CMFEM_ComplexGridFunction *grid_function,
   const char *path,
   int precision);
void CMFEM_ComplexGridFunction_SaveImag(
   const CMFEM_ComplexGridFunction *grid_function,
   const char *path,
   int precision);

CMFEM_END_EXTERN_C

#endif
