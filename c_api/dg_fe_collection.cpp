#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_DG_FECollection, mfem::DG_FECollection);

int OrGaussLegendre(int basis_type)
{
   return basis_type == CMFEM_BASIS_GAUSS_LOBATTO ?
          mfem::BasisType::GaussLobatto :
          mfem::BasisType::GaussLegendre;
}

} // namespace

extern "C" {

   CMFEM_DG_FECollection *CMFEM_DG_FECollection_NewOrderDimBasis(int order,
                                                                 int dim,
                                                                 int basis_type)
   {
      return reinterpret_cast<CMFEM_DG_FECollection *>(
                new mfem::DG_FECollection(order, dim, OrGaussLegendre(basis_type)));
   }

   void CMFEM_DG_FECollection_Delete(CMFEM_DG_FECollection *fec)
   {
      delete cmfem::As<mfem::DG_FECollection>(fec);
   }

} // extern "C"
