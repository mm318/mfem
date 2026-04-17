#ifndef CMFEM_SECOND_ORDER_TIME_DEPENDENT_OPERATOR_H
#define CMFEM_SECOND_ORDER_TIME_DEPENDENT_OPERATOR_H

#include "common.h"

CMFEM_STORAGE(SecondOrderTimeDependentOperator, 72);

typedef void (*CMFEM_SecondOrderTimeDependentOperatorMultCallback)(
   const CMFEM_Vector *u,
   const CMFEM_Vector *du_dt,
   CMFEM_Vector *d2udt2,
   double time,
   void *context);
typedef void (*CMFEM_SecondOrderTimeDependentOperatorImplicitSolveCallback)(
   double fac0,
   double fac1,
   const CMFEM_Vector *u,
   const CMFEM_Vector *du_dt,
   CMFEM_Vector *d2udt2,
   double time,
   void *context);

CMFEM_BEGIN_EXTERN_C

CMFEM_SecondOrderTimeDependentOperator *
CMFEM_SecondOrderTimeDependentOperator_New(
   int size,
   CMFEM_SecondOrderTimeDependentOperatorMultCallback mult,
   void *mult_context,
   CMFEM_SecondOrderTimeDependentOperatorImplicitSolveCallback implicit_solve,
   void *implicit_context);
void CMFEM_SecondOrderTimeDependentOperator_Delete(
   CMFEM_SecondOrderTimeDependentOperator *oper);

CMFEM_END_EXTERN_C

#endif
