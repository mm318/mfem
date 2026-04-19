#ifndef CMFEM_L2_FE_COLLECTION_H
#define CMFEM_L2_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(L2FeCollection, 488);

CMFEM_BEGIN_EXTERN_C

CMFEM_L2FeCollection *CMFEM_L2FeCollection_NewOrderDim(int order, int dim);
CMFEM_L2FeCollection *CMFEM_L2FeCollection_NewOrderDimGl(int order, int dim);
void CMFEM_L2FeCollection_Delete(CMFEM_L2FeCollection *fec);

CMFEM_END_EXTERN_C

#endif
