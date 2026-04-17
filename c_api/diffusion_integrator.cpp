#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_DiffusionIntegrator, mfem::DiffusionIntegrator);

} // namespace

extern "C" {

   CMFEM_DiffusionIntegrator *CMFEM_DiffusionIntegrator_NewCc(
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      return reinterpret_cast<CMFEM_DiffusionIntegrator *>(
                new mfem::DiffusionIntegrator(coef));
   }

   void CMFEM_DiffusionIntegrator_Delete(CMFEM_DiffusionIntegrator *integrator)
   {
      delete cmfem::As<mfem::DiffusionIntegrator>(integrator);
   }

} // extern "C"
