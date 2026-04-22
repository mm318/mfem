#ifndef CMFEM_SECOND_ORDER_ODE_SOLVER_H
#define CMFEM_SECOND_ORDER_ODE_SOLVER_H

#include "../common.h"

CMFEM_STORAGE(SecondOrderODESolver, 8);

CMFEM_BEGIN_EXTERN_C

CMFEM_SecondOrderODESolver *CMFEM_SecondOrderODESolver_NewType(int type);
void CMFEM_SecondOrderODESolver_Delete(CMFEM_SecondOrderODESolver *solver);
void CMFEM_SecondOrderODESolver_Init(
   CMFEM_SecondOrderODESolver *solver,
   CMFEM_SecondOrderTimeDependentOperator *oper);
void CMFEM_SecondOrderODESolver_Step(CMFEM_SecondOrderODESolver *solver,
                                     CMFEM_Vector *u,
                                     CMFEM_Vector *du_dt,
                                     double *time,
                                     double *dt);

CMFEM_END_EXTERN_C

#endif
