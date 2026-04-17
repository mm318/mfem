#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_DG_Interface_FECollection,
                  mfem::DG_Interface_FECollection);

} // namespace

extern "C" {

   CMFEM_DG_Interface_FECollection *CMFEM_DG_Interface_FECollection_NewOrderDim(
      int order,
      int dim)
   {
      return reinterpret_cast<CMFEM_DG_Interface_FECollection *>(
                new mfem::DG_Interface_FECollection(order, dim));
   }

   void CMFEM_DG_Interface_FECollection_Delete(CMFEM_DG_Interface_FECollection
                                               *fec)
   {
      delete cmfem::As<mfem::DG_Interface_FECollection>(fec);
   }

} // extern "C"
