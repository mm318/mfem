#ifndef CMFEM_H1_FE_COLLECTION_H
#define CMFEM_H1_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(H1_FECollection, 512);

CMFEM_BEGIN_EXTERN_C

CMFEM_H1_FECollection *CMFEM_H1_FECollection_NewOrderDim(int order, int dim);
void CMFEM_H1_FECollection_Delete(CMFEM_H1_FECollection *fec);

CMFEM_END_EXTERN_C

#endif
