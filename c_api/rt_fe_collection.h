#ifndef CMFEM_RT_FE_COLLECTION_H
#define CMFEM_RT_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(RtFeCollection, 320);

CMFEM_BEGIN_EXTERN_C

CMFEM_RtFeCollection *CMFEM_RtFeCollection_NewOrderDim(int order, int dim);
void CMFEM_RtFeCollection_Delete(CMFEM_RtFeCollection *fec);

CMFEM_END_EXTERN_C

#endif
