#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_DSmoother, mfem::DSmoother);

} // namespace

extern "C" {

   CMFEM_DSmoother *CMFEM_DSmoother_NewBf(const CMFEM_BilinearForm *bilinear_form)
   {
      auto &matrix = const_cast<mfem::SparseMatrix &>(
                        cmfem::As<const mfem::BilinearForm>(bilinear_form)->SpMat());
      return reinterpret_cast<CMFEM_DSmoother *>(new mfem::DSmoother(matrix));
   }

   CMFEM_DSmoother *CMFEM_DSmoother_NewSm(CMFEM_SparseMatrix *matrix)
   {
      return reinterpret_cast<CMFEM_DSmoother *>(
                new mfem::DSmoother(*cmfem::As<mfem::SparseMatrix>(matrix)));
   }

   void CMFEM_DSmoother_SetIterativeMode(CMFEM_DSmoother *smoother,
                                         int iterative_mode)
   {
      cmfem::As<mfem::DSmoother>(smoother)->iterative_mode =
         iterative_mode != 0;
   }

   void CMFEM_DSmoother_Delete(CMFEM_DSmoother *smoother)
   {
      delete cmfem::As<mfem::DSmoother>(smoother);
   }

} // extern "C"
