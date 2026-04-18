#ifndef CMFEM_ND_R1D_FE_COLLECTION_H
#define CMFEM_ND_R1D_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(NdR1dFeCollection, 184);

CMFEM_BEGIN_EXTERN_C

CMFEM_NdR1dFeCollection *CMFEM_NdR1dFeCollection_NewOrderDim(int order,
                                                             int dim);
void CMFEM_NdR1dFeCollection_Delete(CMFEM_NdR1dFeCollection *fec);

CMFEM_END_EXTERN_C

#endif
