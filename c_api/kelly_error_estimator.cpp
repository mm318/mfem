#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_KellyErrorEstimator, mfem::KellyErrorEstimator);

} // namespace

extern "C" {

   CMFEM_KellyErrorEstimator *
   CMFEM_KellyErrorEstimator_NewDiGfFes(
      CMFEM_DiffusionIntegrator *integrator,
      CMFEM_GridFunction *solution,
      CMFEM_FiniteElementSpace *flux_fespace)
   {
      return reinterpret_cast<CMFEM_KellyErrorEstimator *>(
                new mfem::KellyErrorEstimator(
                   *cmfem::As<mfem::DiffusionIntegrator>(integrator),
                   *cmfem::As<mfem::GridFunction>(solution),
                   cmfem::As<mfem::FiniteElementSpace>(flux_fespace)));
   }

   void CMFEM_KellyErrorEstimator_Delete(CMFEM_KellyErrorEstimator *estimator)
   {
      delete cmfem::As<mfem::KellyErrorEstimator>(estimator);
   }

   void CMFEM_KellyErrorEstimator_Reset(CMFEM_KellyErrorEstimator *estimator)
   {
      cmfem::As<mfem::KellyErrorEstimator>(estimator)->Reset();
   }

} // extern "C"
