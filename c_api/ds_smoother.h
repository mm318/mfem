#ifndef CMFEM_DS_SMOOTHER_H
#define CMFEM_DS_SMOOTHER_H

#include "common.h"

CMFEM_STORAGE(DSmoother, 112);

CMFEM_BEGIN_EXTERN_C

CMFEM_DSmoother *CMFEM_DSmoother_NewBf(const CMFEM_BilinearForm *bilinear_form);
CMFEM_DSmoother *CMFEM_DSmoother_NewSm(CMFEM_SparseMatrix *matrix);
void CMFEM_DSmoother_Delete(CMFEM_DSmoother *smoother);

CMFEM_END_EXTERN_C

#endif
