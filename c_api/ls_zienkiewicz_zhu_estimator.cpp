#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_LSZienkiewiczZhuEstimator,
                  mfem::LSZienkiewiczZhuEstimator);

} // namespace

extern "C" {

   CMFEM_LSZienkiewiczZhuEstimator *
   CMFEM_LSZienkiewiczZhuEstimator_NewDiGf(
      CMFEM_DiffusionIntegrator *integrator,
      CMFEM_GridFunction *solution)
   {
      return reinterpret_cast<CMFEM_LSZienkiewiczZhuEstimator *>(
                new mfem::LSZienkiewiczZhuEstimator(
                   *cmfem::As<mfem::DiffusionIntegrator>(integrator),
                   *cmfem::As<mfem::GridFunction>(solution)));
   }

   void CMFEM_LSZienkiewiczZhuEstimator_Delete(
      CMFEM_LSZienkiewiczZhuEstimator *estimator)
   {
      delete cmfem::As<mfem::LSZienkiewiczZhuEstimator>(estimator);
   }

   void CMFEM_LSZienkiewiczZhuEstimator_SetTichonovRegularization(
      CMFEM_LSZienkiewiczZhuEstimator *estimator,
      double coefficient)
   {
      cmfem::As<mfem::LSZienkiewiczZhuEstimator>(estimator)
      ->SetTichonovRegularization(static_cast<mfem::real_t>(coefficient));
   }

} // extern "C"
