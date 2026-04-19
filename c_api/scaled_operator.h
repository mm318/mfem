#ifndef CMFEM_SCALED_OPERATOR_H
#define CMFEM_SCALED_OPERATOR_H

#include "common.h"

CMFEM_STORAGE(ScaledOperator, 32);

CMFEM_BEGIN_EXTERN_C

CMFEM_ScaledOperator *CMFEM_ScaledOperator_NewGs(
   const CMFEM_GSSmoother *operator_to_scale,
   double scale);
CMFEM_ScaledOperator *CMFEM_ScaledOperator_NewDs(
   const CMFEM_DSmoother *operator_to_scale,
   double scale);
CMFEM_ScaledOperator *CMFEM_ScaledOperator_NewOjs(
   const CMFEM_OperatorJacobiSmoother *operator_to_scale,
   double scale);
void CMFEM_ScaledOperator_Delete(CMFEM_ScaledOperator *scaled_operator);

CMFEM_END_EXTERN_C

#endif
