#ifndef CMFEM_ND_FE_COLLECTION_H
#define CMFEM_ND_FE_COLLECTION_H

#include "../common.h"

CMFEM_STORAGE(NdFeCollection, 320);

CMFEM_BEGIN_EXTERN_C

CMFEM_NdFeCollection *CMFEM_NdFeCollection_NewOrderDim(int order, int dim);
void CMFEM_NdFeCollection_Delete(CMFEM_NdFeCollection *fec);

CMFEM_END_EXTERN_C

#endif
