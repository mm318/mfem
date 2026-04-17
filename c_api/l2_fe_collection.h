#ifndef CMFEM_L2_FE_COLLECTION_H
#define CMFEM_L2_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(L2_FECollection, 488);

CMFEM_BEGIN_EXTERN_C

CMFEM_L2_FECollection *CMFEM_L2_FECollection_NewOrderDim(int order, int dim);
void CMFEM_L2_FECollection_Delete(CMFEM_L2_FECollection *fec);

CMFEM_END_EXTERN_C

#endif
