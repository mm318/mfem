#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ScaledOperator, mfem::ScaledOperator);

} // namespace

extern "C" {

   CMFEM_ScaledOperator *CMFEM_ScaledOperator_NewGs(
      const CMFEM_GSSmoother *operator_to_scale,
      double scale)
   {
      return reinterpret_cast<CMFEM_ScaledOperator *>(
                new mfem::ScaledOperator(
                   cmfem::As<const mfem::GSSmoother>(operator_to_scale),
                   static_cast<mfem::real_t>(scale)));
   }

   CMFEM_ScaledOperator *CMFEM_ScaledOperator_NewDs(
      const CMFEM_DSmoother *operator_to_scale,
      double scale)
   {
      return reinterpret_cast<CMFEM_ScaledOperator *>(
                new mfem::ScaledOperator(
                   cmfem::As<const mfem::DSmoother>(operator_to_scale),
                   static_cast<mfem::real_t>(scale)));
   }

   CMFEM_ScaledOperator *CMFEM_ScaledOperator_NewOjs(
      const CMFEM_OperatorJacobiSmoother *operator_to_scale,
      double scale)
   {
      return reinterpret_cast<CMFEM_ScaledOperator *>(
                new mfem::ScaledOperator(
                   cmfem::As<const mfem::OperatorJacobiSmoother>(operator_to_scale),
                   static_cast<mfem::real_t>(scale)));
   }

   void CMFEM_ScaledOperator_Delete(CMFEM_ScaledOperator *scaled_operator)
   {
      delete cmfem::As<mfem::ScaledOperator>(scaled_operator);
   }

} // extern "C"
