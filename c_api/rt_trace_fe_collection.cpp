#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_RtTraceFeCollection, mfem::RT_Trace_FECollection);

} // namespace

extern "C" {

   CMFEM_RtTraceFeCollection *CMFEM_RtTraceFeCollection_NewOrderDim(int order,
                                                                    int dim)
   {
      return reinterpret_cast<CMFEM_RtTraceFeCollection *>(
                new mfem::RT_Trace_FECollection(order, dim));
   }

   void CMFEM_RtTraceFeCollection_Delete(CMFEM_RtTraceFeCollection *fec)
   {
      delete cmfem::As<mfem::RT_Trace_FECollection>(fec);
   }

} // extern "C"
