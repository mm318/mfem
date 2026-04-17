#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_RtFeCollection, mfem::RT_FECollection);

} // namespace

extern "C" {

   CMFEM_RtFeCollection *CMFEM_RtFeCollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_RtFeCollection *>(
                new mfem::RT_FECollection(order, dim));
   }

   void CMFEM_RtFeCollection_Delete(CMFEM_RtFeCollection *fec)
   {
      delete cmfem::As<mfem::RT_FECollection>(fec);
   }

} // extern "C"
