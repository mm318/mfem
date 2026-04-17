#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_H1_FECollection, mfem::H1_FECollection);

} // namespace

extern "C" {

CMFEM_H1_FECollection *CMFEM_H1_FECollection_NewOrderDim(int order, int dim)
{
   return reinterpret_cast<CMFEM_H1_FECollection *>(
      new mfem::H1_FECollection(order, dim));
}

void CMFEM_H1_FECollection_Delete(CMFEM_H1_FECollection *fec)
{
   delete cmfem::As<mfem::H1_FECollection>(fec);
}

} // extern "C"
