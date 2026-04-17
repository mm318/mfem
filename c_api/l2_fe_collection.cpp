#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_L2_FECollection, mfem::L2_FECollection);

} // namespace

extern "C" {

   CMFEM_L2_FECollection *CMFEM_L2_FECollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_L2_FECollection *>(
                new mfem::L2_FECollection(order, dim));
   }

   void CMFEM_L2_FECollection_Delete(CMFEM_L2_FECollection *fec)
   {
      delete cmfem::As<mfem::L2_FECollection>(fec);
   }

} // extern "C"
