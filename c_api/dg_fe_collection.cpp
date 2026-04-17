#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_DgFeCollection, mfem::DG_FECollection);

int OrGaussLegendre(int basis_type)
{
   return basis_type == CMFEM_BASIS_GAUSS_LOBATTO ?
          mfem::BasisType::GaussLobatto :
          mfem::BasisType::GaussLegendre;
}

} // namespace

extern "C" {

   CMFEM_DgFeCollection *CMFEM_DgFeCollection_NewOrderDimBasis(int order,
                                                               int dim,
                                                               int basis_type)
   {
      return reinterpret_cast<CMFEM_DgFeCollection *>(
                new mfem::DG_FECollection(order, dim, OrGaussLegendre(basis_type)));
   }

   void CMFEM_DgFeCollection_Delete(CMFEM_DgFeCollection *fec)
   {
      delete cmfem::As<mfem::DG_FECollection>(fec);
   }

} // extern "C"
