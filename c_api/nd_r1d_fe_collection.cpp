#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_NdR1dFeCollection, mfem::ND_R1D_FECollection);

} // namespace

extern "C" {

   CMFEM_NdR1dFeCollection *CMFEM_NdR1dFeCollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_NdR1dFeCollection *>(
                new mfem::ND_R1D_FECollection(order, dim));
   }

   void CMFEM_NdR1dFeCollection_Delete(CMFEM_NdR1dFeCollection *fec)
   {
      delete cmfem::As<mfem::ND_R1D_FECollection>(fec);
   }

} // extern "C"
