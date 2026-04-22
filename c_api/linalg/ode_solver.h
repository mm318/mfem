#ifndef CMFEM_ODE_SOLVER_H
#define CMFEM_ODE_SOLVER_H

#include "../common.h"

CMFEM_STORAGE(ODESolver, 8);

CMFEM_BEGIN_EXTERN_C

CMFEM_ODESolver *CMFEM_ODESolver_NewType(int type);
CMFEM_ODESolver *CMFEM_ODESolver_NewImexType(int type);
void CMFEM_ODESolver_Delete(CMFEM_ODESolver *solver);
void CMFEM_ODESolver_Init(CMFEM_ODESolver *solver,
                          CMFEM_TimeDependentOperator *oper);
void CMFEM_ODESolver_Step(CMFEM_ODESolver *solver,
                          CMFEM_Vector *u,
                          double *time,
                          double *dt);

CMFEM_END_EXTERN_C

#endif
