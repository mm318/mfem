#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_GSSmoother, mfem::GSSmoother);

} // namespace

extern "C" {

   CMFEM_GSSmoother *CMFEM_GSSmoother_NewSm(CMFEM_SparseMatrix *matrix)
   {
      return reinterpret_cast<CMFEM_GSSmoother *>(
                new mfem::GSSmoother(*cmfem::As<mfem::SparseMatrix>(matrix)));
   }

   CMFEM_GSSmoother *CMFEM_GSSmoother_NewOp(CMFEM_OperatorPtr *op)
   {
      return reinterpret_cast<CMFEM_GSSmoother *>(
                new mfem::GSSmoother(cmfem::SparseMatrixFromOperator(op)));
   }

   void CMFEM_GSSmoother_Delete(CMFEM_GSSmoother *smoother)
   {
      delete cmfem::As<mfem::GSSmoother>(smoother);
   }

} // extern "C"
