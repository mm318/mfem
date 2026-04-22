#ifndef CMFEM_INCOMPRESSIBLE_RUBBER_H
#define CMFEM_INCOMPRESSIBLE_RUBBER_H

#include "../../../c_api/common.h"

CMFEM_BEGIN_EXTERN_C

void *CMFEM_Ex19RubberOperator_New(
   CMFEM_FiniteElementSpace *displacement_space,
   CMFEM_FiniteElementSpace *pressure_space,
   const CMFEM_ArrayInt *ess_bdr_u,
   const CMFEM_ArrayInt *ess_bdr_p,
   double rel_tol,
   double abs_tol,
   int max_iter,
   double mu);
void CMFEM_Ex19RubberOperator_Delete(void *oper);
void CMFEM_Ex19RubberOperator_Solve(void *oper,
                                    CMFEM_GridFunction *position,
                                    CMFEM_GridFunction *pressure);

CMFEM_END_EXTERN_C

#endif
