#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_NdFeCollection, mfem::ND_FECollection);

} // namespace

extern "C" {

   CMFEM_NdFeCollection *CMFEM_NdFeCollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_NdFeCollection *>(
                new mfem::ND_FECollection(order, dim));
   }

   void CMFEM_NdFeCollection_Delete(CMFEM_NdFeCollection *fec)
   {
      delete cmfem::As<mfem::ND_FECollection>(fec);
   }

} // extern "C"
