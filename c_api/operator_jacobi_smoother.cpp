#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_OperatorJacobiSmoother, mfem::OperatorJacobiSmoother);

} // namespace

extern "C" {

   CMFEM_OperatorJacobiSmoother *
   CMFEM_OperatorJacobiSmoother_NewBilinearForm(const CMFEM_BilinearForm
                                                *bilinear_form,
                                                const CMFEM_ArrayInt *ess_tdof_list)
   {
      return reinterpret_cast<CMFEM_OperatorJacobiSmoother *>(
                new mfem::OperatorJacobiSmoother(
                   *cmfem::As<const mfem::BilinearForm>(bilinear_form),
                   cmfem::ArrayIntRef(ess_tdof_list)));
   }

   void CMFEM_OperatorJacobiSmoother_Delete(CMFEM_OperatorJacobiSmoother *smoother)
   {
      delete cmfem::As<mfem::OperatorJacobiSmoother>(smoother);
   }

} // extern "C"
