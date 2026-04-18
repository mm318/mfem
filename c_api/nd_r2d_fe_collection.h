#ifndef CMFEM_ND_R2D_FE_COLLECTION_H
#define CMFEM_ND_R2D_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(NdR2dFeCollection, 200);

CMFEM_BEGIN_EXTERN_C

CMFEM_NdR2dFeCollection *CMFEM_NdR2dFeCollection_NewOrderDim(int order,
                                                             int dim);
void CMFEM_NdR2dFeCollection_Delete(CMFEM_NdR2dFeCollection *fec);

CMFEM_END_EXTERN_C

#endif
