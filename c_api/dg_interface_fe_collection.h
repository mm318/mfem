#ifndef CMFEM_DG_INTERFACE_FE_COLLECTION_H
#define CMFEM_DG_INTERFACE_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(DgInterfaceFeCollection, 320);

CMFEM_BEGIN_EXTERN_C

CMFEM_DgInterfaceFeCollection *CMFEM_DgInterfaceFeCollection_NewOrderDim(
   int order, int dim);
void CMFEM_DgInterfaceFeCollection_Delete(CMFEM_DgInterfaceFeCollection
                                          *fec);

CMFEM_END_EXTERN_C

#endif
