#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_RT_FECollection, mfem::RT_FECollection);

} // namespace

extern "C" {

   CMFEM_RT_FECollection *CMFEM_RT_FECollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_RT_FECollection *>(
                new mfem::RT_FECollection(order, dim));
   }

   void CMFEM_RT_FECollection_Delete(CMFEM_RT_FECollection *fec)
   {
      delete cmfem::As<mfem::RT_FECollection>(fec);
   }

} // extern "C"
