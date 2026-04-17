#ifndef CMFEM_DG_INTERFACE_FE_COLLECTION_H
#define CMFEM_DG_INTERFACE_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(DG_Interface_FECollection, 320);

CMFEM_BEGIN_EXTERN_C

CMFEM_DG_Interface_FECollection *CMFEM_DG_Interface_FECollection_NewOrderDim(
   int order, int dim);
void CMFEM_DG_Interface_FECollection_Delete(CMFEM_DG_Interface_FECollection
                                            *fec);

CMFEM_END_EXTERN_C

#endif
