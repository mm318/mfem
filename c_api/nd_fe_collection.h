#ifndef CMFEM_ND_FE_COLLECTION_H
#define CMFEM_ND_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(ND_FECollection, 320);

CMFEM_BEGIN_EXTERN_C

CMFEM_ND_FECollection *CMFEM_ND_FECollection_NewOrderDim(int order, int dim);
void CMFEM_ND_FECollection_Delete(CMFEM_ND_FECollection *fec);

CMFEM_END_EXTERN_C

#endif
