#ifndef CMFEM_SCHUR_CONSTRAINED_SOLVER_H
#define CMFEM_SCHUR_CONSTRAINED_SOLVER_H

#include "common.h"

CMFEM_STORAGE(SchurConstrainedSolver, 360);

CMFEM_BEGIN_EXTERN_C

CMFEM_SchurConstrainedSolver *CMFEM_SchurConstrainedSolver_NewSmSmGs(
   CMFEM_SparseMatrix *primal_operator,
   CMFEM_SparseMatrix *constraint_operator,
   CMFEM_GSSmoother *primal_pc);
void CMFEM_SchurConstrainedSolver_Delete(
   CMFEM_SchurConstrainedSolver *solver);
void CMFEM_SchurConstrainedSolver_SetRelTol(
   CMFEM_SchurConstrainedSolver *solver,
   double rel_tol);
void CMFEM_SchurConstrainedSolver_SetMaxIter(
   CMFEM_SchurConstrainedSolver *solver,
   int max_iter);
void CMFEM_SchurConstrainedSolver_SetPrintLevel(
   CMFEM_SchurConstrainedSolver *solver,
   int print_level);
void CMFEM_SchurConstrainedSolver_Mult(
   const CMFEM_SchurConstrainedSolver *solver,
   const CMFEM_Vector *rhs,
   CMFEM_Vector *solution);

CMFEM_END_EXTERN_C

#endif
