#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ND_FECollection, mfem::ND_FECollection);

} // namespace

extern "C" {

   CMFEM_ND_FECollection *CMFEM_ND_FECollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_ND_FECollection *>(
                new mfem::ND_FECollection(order, dim));
   }

   void CMFEM_ND_FECollection_Delete(CMFEM_ND_FECollection *fec)
   {
      delete cmfem::As<mfem::ND_FECollection>(fec);
   }

} // extern "C"
