#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ElasticityIntegrator, mfem::ElasticityIntegrator);

} // namespace

extern "C" {

   CMFEM_ElasticityIntegrator *CMFEM_ElasticityIntegrator_NewPwcPwc(
      const CMFEM_PWConstCoefficient *lambda,
      const CMFEM_PWConstCoefficient *mu)
   {
      auto &lambda_ref = const_cast<mfem::PWConstCoefficient &>(
                            *cmfem::As<const mfem::PWConstCoefficient>(lambda));
      auto &mu_ref = const_cast<mfem::PWConstCoefficient &>(
                        *cmfem::As<const mfem::PWConstCoefficient>(mu));
      return reinterpret_cast<CMFEM_ElasticityIntegrator *>(
                new mfem::ElasticityIntegrator(lambda_ref, mu_ref));
   }

   void CMFEM_ElasticityIntegrator_Delete(CMFEM_ElasticityIntegrator *integrator)
   {
      delete cmfem::As<mfem::ElasticityIntegrator>(integrator);
   }

} // extern "C"
