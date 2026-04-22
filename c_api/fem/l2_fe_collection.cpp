#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_L2FeCollection, mfem::L2_FECollection);

} // namespace

extern "C" {

   CMFEM_L2FeCollection *CMFEM_L2FeCollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_L2FeCollection *>(
                new mfem::L2_FECollection(order, dim));
   }

   CMFEM_L2FeCollection *CMFEM_L2FeCollection_NewOrderDimGl(int order, int dim)
   {
      return reinterpret_cast<CMFEM_L2FeCollection *>(
                new mfem::L2_FECollection(order,
                                          dim,
                                          mfem::BasisType::GaussLobatto));
   }

   void CMFEM_L2FeCollection_Delete(CMFEM_L2FeCollection *fec)
   {
      delete cmfem::As<mfem::L2_FECollection>(fec);
   }

} // extern "C"
