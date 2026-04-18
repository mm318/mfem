#ifndef CMFEM_RT_TRACE_FE_COLLECTION_H
#define CMFEM_RT_TRACE_FE_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(RtTraceFeCollection, 320);

CMFEM_BEGIN_EXTERN_C

CMFEM_RtTraceFeCollection *CMFEM_RtTraceFeCollection_NewOrderDim(int order,
                                                                 int dim);
void CMFEM_RtTraceFeCollection_Delete(CMFEM_RtTraceFeCollection *fec);

CMFEM_END_EXTERN_C

#endif
