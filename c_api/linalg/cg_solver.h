#ifndef CMFEM_CG_SOLVER_H
#define CMFEM_CG_SOLVER_H

#include "../common.h"

CMFEM_STORAGE(CGSolver, 232);

CMFEM_BEGIN_EXTERN_C

CMFEM_CGSolver *CMFEM_CGSolver_New(void);
void CMFEM_CGSolver_Delete(CMFEM_CGSolver *solver);
void CMFEM_CGSolver_SetOperatorSm(CMFEM_CGSolver *solver,
                                  const CMFEM_SparseMatrix *matrix);
void CMFEM_CGSolver_SetPrintLevel(CMFEM_CGSolver *solver, int print_level);
void CMFEM_CGSolver_SetRelTol(CMFEM_CGSolver *solver, double rel_tol);
void CMFEM_CGSolver_SetMaxIter(CMFEM_CGSolver *solver, int max_iter);
void CMFEM_CGSolver_SetIterativeMode(CMFEM_CGSolver *solver,
                                     int iterative_mode);

CMFEM_END_EXTERN_C

#endif
