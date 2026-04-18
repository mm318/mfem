#ifndef CMFEM_DISCRETE_LINEAR_OPERATOR_H
#define CMFEM_DISCRETE_LINEAR_OPERATOR_H

#include "common.h"

CMFEM_STORAGE(DiscreteLinearOperator, 504);

CMFEM_BEGIN_EXTERN_C

CMFEM_DiscreteLinearOperator *CMFEM_DiscreteLinearOperator_New(
   CMFEM_FiniteElementSpace *domain_fespace,
   CMFEM_FiniteElementSpace *range_fespace);
void CMFEM_DiscreteLinearOperator_Delete(
   CMFEM_DiscreteLinearOperator *discrete_linear_operator);
void CMFEM_DiscreteLinearOperator_AddDomainInterpolatorGi(
   CMFEM_DiscreteLinearOperator *discrete_linear_operator);
void CMFEM_DiscreteLinearOperator_AddDomainInterpolatorCi(
   CMFEM_DiscreteLinearOperator *discrete_linear_operator);
void CMFEM_DiscreteLinearOperator_AddDomainInterpolatorDi(
   CMFEM_DiscreteLinearOperator *discrete_linear_operator);
void CMFEM_DiscreteLinearOperator_Assemble(
   CMFEM_DiscreteLinearOperator *discrete_linear_operator);
void CMFEM_DiscreteLinearOperator_Mult(
   const CMFEM_DiscreteLinearOperator *discrete_linear_operator,
   const CMFEM_Vector *x,
   CMFEM_Vector *y);

CMFEM_END_EXTERN_C

#endif
