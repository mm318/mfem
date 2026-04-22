#ifndef CMFEM_TIME_DEPENDENT_OPERATOR_H
#define CMFEM_TIME_DEPENDENT_OPERATOR_H

#include "../common.h"

CMFEM_STORAGE(TimeDependentOperator, 72);

enum CMFEM_TimeDependentOperatorType
{
   CMFEM_TIME_DEPENDENT_OPERATOR_EXPLICIT = 0,
   CMFEM_TIME_DEPENDENT_OPERATOR_IMPLICIT = 1,
   CMFEM_TIME_DEPENDENT_OPERATOR_HOMOGENEOUS = 2,
};

enum CMFEM_ImplicitVariableType
{
   CMFEM_IMPLICIT_VARIABLE_SLOPE = 0,
   CMFEM_IMPLICIT_VARIABLE_STATE = 1,
};

typedef void (*CMFEM_TimeDependentOperatorMultCallback)(
   const CMFEM_Vector *u,
   CMFEM_Vector *du_dt,
   double time,
   void *context);
typedef void (*CMFEM_TimeDependentOperatorImplicitSolveCallback)(
   double gamma,
   const CMFEM_Vector *u,
   CMFEM_Vector *k,
   double time,
   void *context);

CMFEM_BEGIN_EXTERN_C

CMFEM_TimeDependentOperator *CMFEM_TimeDependentOperator_New(
   int size,
   enum CMFEM_TimeDependentOperatorType type,
   CMFEM_TimeDependentOperatorMultCallback mult,
   void *mult_context,
   CMFEM_TimeDependentOperatorImplicitSolveCallback implicit_solve,
   void *implicit_context);
CMFEM_TimeDependentOperator *CMFEM_TimeDependentOperator_NewImex(
   int size,
   CMFEM_TimeDependentOperatorMultCallback mult_term_1,
   void *mult_term_1_context,
   CMFEM_TimeDependentOperatorImplicitSolveCallback implicit_solve_term_2,
   void *implicit_solve_term_2_context);
void CMFEM_TimeDependentOperator_Delete(CMFEM_TimeDependentOperator *oper);
void CMFEM_TimeDependentOperator_SetImplicitVariableType(
   CMFEM_TimeDependentOperator *oper,
   enum CMFEM_ImplicitVariableType variable_type);

CMFEM_END_EXTERN_C

#endif
