#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_NdR2dFeCollection, mfem::ND_R2D_FECollection);

} // namespace

extern "C" {

   CMFEM_NdR2dFeCollection *CMFEM_NdR2dFeCollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_NdR2dFeCollection *>(
                new mfem::ND_R2D_FECollection(order, dim));
   }

   void CMFEM_NdR2dFeCollection_Delete(CMFEM_NdR2dFeCollection *fec)
   {
      delete cmfem::As<mfem::ND_R2D_FECollection>(fec);
   }

} // extern "C"
