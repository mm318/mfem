#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_DSmoother, mfem::DSmoother);

} // namespace

extern "C" {

   CMFEM_DSmoother *CMFEM_DSmoother_NewSm(CMFEM_SparseMatrix *matrix)
   {
      return reinterpret_cast<CMFEM_DSmoother *>(
                new mfem::DSmoother(*cmfem::As<mfem::SparseMatrix>(matrix)));
   }

   void CMFEM_DSmoother_Delete(CMFEM_DSmoother *smoother)
   {
      delete cmfem::As<mfem::DSmoother>(smoother);
   }

} // extern "C"
