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

   void CMFEM_LinearForm_AddDomainIntegrator_DomainLFIntegrator_ConstantCoefficient(
      CMFEM_LinearForm *linear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::LinearForm>(linear_form)->AddDomainIntegrator(
         new mfem::DomainLFIntegrator(coef));
   }

   void CMFEM_LinearForm_AddDomainIntegrator_VectorFEDomainLFIntegrator(
      CMFEM_LinearForm *linear_form,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      cmfem::As<mfem::LinearForm>(linear_form)->AddDomainIntegrator(
         new mfem::VectorFEDomainLFIntegrator(coef));
   }

   void CMFEM_LinearForm_AddBoundaryIntegrator_VectorBoundaryLFIntegrator(
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

} // extern "C"
