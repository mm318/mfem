#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_H1FeCollection, mfem::H1_FECollection);

} // namespace

extern "C" {

   CMFEM_H1FeCollection *CMFEM_H1FeCollection_NewOrderDim(int order, int dim)
   {
      return reinterpret_cast<CMFEM_H1FeCollection *>(
                new mfem::H1_FECollection(order, dim));
   }

   void CMFEM_H1FeCollection_Delete(CMFEM_H1FeCollection *fec)
   {
      delete cmfem::As<mfem::H1_FECollection>(fec);
   }

} // extern "C"
