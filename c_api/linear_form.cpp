#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_LinearForm, mfem::LinearForm);

} // namespace

extern "C" {

   CMFEM_LinearForm *CMFEM_LinearForm_New(CMFEM_FiniteElementSpace *fespace)
   {
      return reinterpret_cast<CMFEM_LinearForm *>(
                new mfem::LinearForm(cmfem::As<mfem::FiniteElementSpace>(fespace)));
   }

   void CMFEM_LinearForm_Delete(CMFEM_LinearForm *linear_form)
   {
      delete cmfem::As<mfem::LinearForm>(linear_form);
   }

   void CMFEM_LinearForm_AddDomainIntegratorDliCc(
      CMFEM_LinearForm *linear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::LinearForm>(linear_form)->AddDomainIntegrator(
         new mfem::DomainLFIntegrator(coef));
   }

   void CMFEM_LinearForm_AddDomainIntegratorDliFc(
      CMFEM_LinearForm *linear_form,
      const CMFEM_FunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      cmfem::As<mfem::LinearForm>(linear_form)->AddDomainIntegrator(
         new mfem::DomainLFIntegrator(coef));
   }

   void
   CMFEM_LinearForm_AddDomainIntegratorDliCcMarker(
      CMFEM_LinearForm *linear_form,
      const CMFEM_ConstantCoefficient *coefficient,
      const CMFEM_ArrayInt *marker)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::LinearForm>(linear_form)->AddDomainIntegrator(
         new mfem::DomainLFIntegrator(coef),
         const_cast<mfem::Array<int> &>(cmfem::ArrayIntRef(marker)));
   }

   void CMFEM_LinearForm_AddBdrFaceIntegratorBfiFcVfc(
      CMFEM_LinearForm *linear_form,
      const CMFEM_FunctionCoefficient *coefficient,
      const CMFEM_VectorFunctionCoefficient *velocity,
      double alpha)
   {
      auto &coef = const_cast<mfem::FunctionCoefficient &>(
                      *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      auto &vel = const_cast<mfem::VectorFunctionCoefficient &>(
                     *cmfem::As<const mfem::VectorFunctionCoefficient>(
                        velocity));
      cmfem::As<mfem::LinearForm>(linear_form)->AddBdrFaceIntegrator(
         new mfem::BoundaryFlowIntegrator(coef,
                                          vel,
                                          static_cast<mfem::real_t>(alpha)));
   }

   void CMFEM_LinearForm_AddDomainIntegratorVfd(
      CMFEM_LinearForm *linear_form,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      cmfem::As<mfem::LinearForm>(linear_form)->AddDomainIntegrator(
         new mfem::VectorFEDomainLFIntegrator(coef));
   }

   void CMFEM_LinearForm_AddBdrFaceIntegratorDgl(
      CMFEM_LinearForm *linear_form,
      const CMFEM_ConstantCoefficient *u_coeff,
      const CMFEM_ConstantCoefficient *q_coeff,
      double sigma,
      double kappa)
   {
      auto &u_ref = const_cast<mfem::ConstantCoefficient &>(
                       *cmfem::As<const mfem::ConstantCoefficient>(u_coeff));
      auto &q_ref = const_cast<mfem::ConstantCoefficient &>(
                       *cmfem::As<const mfem::ConstantCoefficient>(q_coeff));
      cmfem::As<mfem::LinearForm>(linear_form)->AddBdrFaceIntegrator(
         new mfem::DGDirichletLFIntegrator(
            u_ref, q_ref, static_cast<mfem::real_t>(sigma),
            static_cast<mfem::real_t>(kappa)));
   }

   void CMFEM_LinearForm_AddBdrFaceIntegratorDgeliVfcPwcPwcAi(
      CMFEM_LinearForm *linear_form,
      const CMFEM_VectorFunctionCoefficient *coefficient,
      const CMFEM_PWConstCoefficient *lambda,
      const CMFEM_PWConstCoefficient *mu,
      double alpha,
      double kappa,
      const CMFEM_ArrayInt *marker)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(
                         coefficient));
      auto &lambda_ref = const_cast<mfem::PWConstCoefficient &>(
                            *cmfem::As<const mfem::PWConstCoefficient>(lambda));
      auto &mu_ref = const_cast<mfem::PWConstCoefficient &>(
                        *cmfem::As<const mfem::PWConstCoefficient>(mu));
      cmfem::As<mfem::LinearForm>(linear_form)->AddBdrFaceIntegrator(
         new mfem::DGElasticityDirichletLFIntegrator(
            coef,
            lambda_ref,
            mu_ref,
            static_cast<mfem::real_t>(alpha),
            static_cast<mfem::real_t>(kappa)),
         const_cast<mfem::Array<int> &>(cmfem::ArrayIntRef(marker)));
   }

   void CMFEM_LinearForm_AddBoundaryIntegratorVbl(
      CMFEM_LinearForm *linear_form,
      const CMFEM_VectorArrayCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorArrayCoefficient &>(
                      *cmfem::As<const mfem::VectorArrayCoefficient>(coefficient));
      cmfem::As<mfem::LinearForm>(linear_form)->AddBoundaryIntegrator(
         new mfem::VectorBoundaryLFIntegrator(coef));
   }

   void CMFEM_LinearForm_Assemble(CMFEM_LinearForm *linear_form)
   {
      cmfem::As<mfem::LinearForm>(linear_form)->Assemble();
   }

   void CMFEM_LinearForm_CopyToVector(const CMFEM_LinearForm *linear_form,
                                      CMFEM_Vector *vector)
   {
      cmfem::VectorRef(vector) = *cmfem::As<const mfem::LinearForm>(
                                    linear_form);
   }

   void CMFEM_LinearForm_Update(CMFEM_LinearForm *linear_form)
   {
      cmfem::As<mfem::LinearForm>(linear_form)->Update();
   }

} // extern "C"
