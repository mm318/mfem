#ifndef CMFEM_RT_FE_COLLECTION_H
#define CMFEM_RT_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(RT_FECollection, 320);

CMFEM_BEGIN_EXTERN_C

CMFEM_RT_FECollection *CMFEM_RT_FECollection_NewOrderDim(int order, int dim);
void CMFEM_RT_FECollection_Delete(CMFEM_RT_FECollection *fec);

CMFEM_END_EXTERN_C

#endif
