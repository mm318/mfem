#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_DgInterfaceFeCollection,
                  mfem::DG_Interface_FECollection);

} // namespace

extern "C" {

   CMFEM_DgInterfaceFeCollection *CMFEM_DgInterfaceFeCollection_NewOrderDim(
      int order,
      int dim)
   {
      return reinterpret_cast<CMFEM_DgInterfaceFeCollection *>(
                new mfem::DG_Interface_FECollection(order, dim));
   }

   void CMFEM_DgInterfaceFeCollection_Delete(CMFEM_DgInterfaceFeCollection
                                             *fec)
   {
      delete cmfem::As<mfem::DG_Interface_FECollection>(fec);
   }

} // extern "C"
