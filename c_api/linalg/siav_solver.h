#ifndef CMFEM_SIAV_SOLVER_H
#define CMFEM_SIAV_SOLVER_H

#include "../common.h"

CMFEM_STORAGE(SIAVSolver, 176);

typedef void (*CMFEM_OperatorMultCallback)(const CMFEM_Vector *x,
                                           CMFEM_Vector *y,
                                           void *context);
typedef void (*CMFEM_TimeDependentOperatorMultCallback)(const CMFEM_Vector *x,
                                                        CMFEM_Vector *y,
                                                        double time,
                                                        void *context);

CMFEM_BEGIN_EXTERN_C

CMFEM_SIAVSolver *CMFEM_SIAVSolver_NewOrder(int order);
void CMFEM_SIAVSolver_Delete(CMFEM_SIAVSolver *solver);
void CMFEM_SIAVSolver_Init(CMFEM_SIAVSolver *solver,
                           int size,
                           CMFEM_OperatorMultCallback p_mult,
                           void *p_context,
                           CMFEM_TimeDependentOperatorMultCallback f_mult,
                           void *f_context);
void CMFEM_SIAVSolver_Step(CMFEM_SIAVSolver *solver,
                           CMFEM_Vector *q,
                           CMFEM_Vector *p,
                           double *time,
                           double *dt);

CMFEM_END_EXTERN_C

#endif
