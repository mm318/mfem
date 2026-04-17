#ifndef CMFEM_OPERATOR_JACOBI_SMOOTHER_H
#define CMFEM_OPERATOR_JACOBI_SMOOTHER_H

#include "common.h"

CMFEM_STORAGE(OperatorJacobiSmoother, 144);

CMFEM_BEGIN_EXTERN_C

CMFEM_OperatorJacobiSmoother *CMFEM_OperatorJacobiSmoother_NewBilinearForm(
   const CMFEM_BilinearForm *bilinear_form, const CMFEM_ArrayInt *ess_tdof_list);
void CMFEM_OperatorJacobiSmoother_Delete(CMFEM_OperatorJacobiSmoother
                                         *smoother);

CMFEM_END_EXTERN_C

#endif
