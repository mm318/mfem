#ifndef CMFEM_GS_SMOOTHER_H
#define CMFEM_GS_SMOOTHER_H

#include "common.h"

CMFEM_STORAGE(GSSmoother, 56);

CMFEM_BEGIN_EXTERN_C

CMFEM_GSSmoother *CMFEM_GSSmoother_NewSm(CMFEM_SparseMatrix *matrix);
CMFEM_GSSmoother *CMFEM_GSSmoother_NewOp(CMFEM_OperatorPtr *op);
void CMFEM_GSSmoother_Delete(CMFEM_GSSmoother *smoother);

CMFEM_END_EXTERN_C

#endif
