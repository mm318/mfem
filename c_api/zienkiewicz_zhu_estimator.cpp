#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ZienkiewiczZhuEstimator, mfem::ZienkiewiczZhuEstimator);

} // namespace

extern "C" {

   CMFEM_ZienkiewiczZhuEstimator *
   CMFEM_ZienkiewiczZhuEstimator_NewDiGfFes(
      CMFEM_DiffusionIntegrator *integrator,
      CMFEM_GridFunction *solution,
      CMFEM_FiniteElementSpace *flux_fespace)
   {
      return reinterpret_cast<CMFEM_ZienkiewiczZhuEstimator *>(
                new mfem::ZienkiewiczZhuEstimator(
                   *cmfem::As<mfem::DiffusionIntegrator>(integrator),
                   *cmfem::As<mfem::GridFunction>(solution),
                   cmfem::As<mfem::FiniteElementSpace>(flux_fespace)));
   }

   CMFEM_ZienkiewiczZhuEstimator *
   CMFEM_ZienkiewiczZhuEstimator_NewEiGfFes(
      CMFEM_ElasticityIntegrator *integrator,
      CMFEM_GridFunction *solution,
      CMFEM_FiniteElementSpace *flux_fespace)
   {
      return reinterpret_cast<CMFEM_ZienkiewiczZhuEstimator *>(
                new mfem::ZienkiewiczZhuEstimator(
                   *cmfem::As<mfem::ElasticityIntegrator>(integrator),
                   *cmfem::As<mfem::GridFunction>(solution),
                   cmfem::As<mfem::FiniteElementSpace>(flux_fespace)));
   }

   void CMFEM_ZienkiewiczZhuEstimator_Delete(
      CMFEM_ZienkiewiczZhuEstimator *estimator)
   {
      delete cmfem::As<mfem::ZienkiewiczZhuEstimator>(estimator);
   }

   void CMFEM_ZienkiewiczZhuEstimator_SetFluxAveraging(
      CMFEM_ZienkiewiczZhuEstimator *estimator,
      int flux_averaging)
   {
      cmfem::As<mfem::ZienkiewiczZhuEstimator>(estimator)->SetFluxAveraging(
         flux_averaging);
   }

   void CMFEM_ZienkiewiczZhuEstimator_SetAnisotropic(
      CMFEM_ZienkiewiczZhuEstimator *estimator,
      int anisotropic)
   {
      cmfem::As<mfem::ZienkiewiczZhuEstimator>(estimator)->SetAnisotropic(
         anisotropic != 0);
   }

} // extern "C"
