#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_SesquilinearForm, mfem::SesquilinearForm);

mfem::ComplexOperator::Convention ToConvention(
   enum CMFEM_ComplexConvention convention)
{
   return convention == CMFEM_ComplexConventionBlockSymmetric ?
          mfem::ComplexOperator::BLOCK_SYMMETRIC :
          mfem::ComplexOperator::HERMITIAN;
}

} // namespace

extern "C" {

   CMFEM_SesquilinearForm *CMFEM_SesquilinearForm_NewFesConv(
      CMFEM_FiniteElementSpace *fespace,
      enum CMFEM_ComplexConvention convention)
   {
      return reinterpret_cast<CMFEM_SesquilinearForm *>(
                new mfem::SesquilinearForm(
                   cmfem::As<mfem::FiniteElementSpace>(fespace),
                   ToConvention(convention)));
   }

   void CMFEM_SesquilinearForm_Delete(
      CMFEM_SesquilinearForm *sesquilinear_form)
   {
      delete cmfem::As<mfem::SesquilinearForm>(sesquilinear_form);
   }

   void CMFEM_SesquilinearForm_SetAssemblyLevelPartial(
      CMFEM_SesquilinearForm *sesquilinear_form)
   {
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)
      ->SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
   }

   void CMFEM_SesquilinearForm_AddDomainIntegratorDiCcNull(
      CMFEM_SesquilinearForm *sesquilinear_form,
      const CMFEM_ConstantCoefficient *real_coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(real_coefficient));
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)->AddDomainIntegrator(
         new mfem::DiffusionIntegrator(coef),
         nullptr);
   }

   void CMFEM_SesquilinearForm_AddDomainIntegratorMiCcMiCc(
      CMFEM_SesquilinearForm *sesquilinear_form,
      const CMFEM_ConstantCoefficient *real_coefficient,
      const CMFEM_ConstantCoefficient *imag_coefficient)
   {
      auto &real_coef = const_cast<mfem::ConstantCoefficient &>(
                           *cmfem::As<const mfem::ConstantCoefficient>(real_coefficient));
      auto &imag_coef = const_cast<mfem::ConstantCoefficient &>(
                           *cmfem::As<const mfem::ConstantCoefficient>(imag_coefficient));
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)->AddDomainIntegrator(
         new mfem::MassIntegrator(real_coef),
         new mfem::MassIntegrator(imag_coef));
   }

   void CMFEM_SesquilinearForm_AddDomainIntegratorCciCcNull(
      CMFEM_SesquilinearForm *sesquilinear_form,
      const CMFEM_ConstantCoefficient *real_coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(real_coefficient));
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)->AddDomainIntegrator(
         new mfem::CurlCurlIntegrator(coef),
         nullptr);
   }

   void CMFEM_SesquilinearForm_AddDomainIntegratorVmiCcVmiCc(
      CMFEM_SesquilinearForm *sesquilinear_form,
      const CMFEM_ConstantCoefficient *real_coefficient,
      const CMFEM_ConstantCoefficient *imag_coefficient)
   {
      auto &real_coef = const_cast<mfem::ConstantCoefficient &>(
                           *cmfem::As<const mfem::ConstantCoefficient>(real_coefficient));
      auto &imag_coef = const_cast<mfem::ConstantCoefficient &>(
                           *cmfem::As<const mfem::ConstantCoefficient>(imag_coefficient));
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)->AddDomainIntegrator(
         new mfem::VectorFEMassIntegrator(real_coef),
         new mfem::VectorFEMassIntegrator(imag_coef));
   }

   void CMFEM_SesquilinearForm_AddDomainIntegratorDdiCcNull(
      CMFEM_SesquilinearForm *sesquilinear_form,
      const CMFEM_ConstantCoefficient *real_coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(real_coefficient));
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)->AddDomainIntegrator(
         new mfem::DivDivIntegrator(coef),
         nullptr);
   }

   void CMFEM_SesquilinearForm_Assemble(
      CMFEM_SesquilinearForm *sesquilinear_form)
   {
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)->Assemble();
   }

   void CMFEM_SesquilinearForm_AssembleSkipZeros(
      CMFEM_SesquilinearForm *sesquilinear_form,
      int skip_zeros)
   {
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)->Assemble(skip_zeros);
   }

   void CMFEM_SesquilinearForm_FormLinearSystemOp(
      CMFEM_SesquilinearForm *sesquilinear_form,
      const CMFEM_ArrayInt *ess_tdof_list,
      CMFEM_ComplexGridFunction *x,
      CMFEM_ComplexLinearForm *b,
      CMFEM_OperatorPtr *A,
      CMFEM_Vector *X,
      CMFEM_Vector *B)
   {
      cmfem::As<mfem::SesquilinearForm>(sesquilinear_form)->FormLinearSystem(
         cmfem::ArrayIntRef(ess_tdof_list),
         *cmfem::As<mfem::ComplexGridFunction>(x),
         *cmfem::As<mfem::ComplexLinearForm>(b),
         cmfem::OperatorPtrRef(A),
         cmfem::VectorRef(X),
         cmfem::VectorRef(B));
   }

   void CMFEM_SesquilinearForm_RecoverFEMSolution(
      const CMFEM_SesquilinearForm *sesquilinear_form,
      const CMFEM_Vector *X,
      const CMFEM_ComplexLinearForm *b,
      CMFEM_ComplexGridFunction *x)
   {
      const_cast<mfem::SesquilinearForm *>(
         cmfem::As<const mfem::SesquilinearForm>(sesquilinear_form))
      ->RecoverFEMSolution(cmfem::VectorRef(X),
                           *cmfem::As<const mfem::ComplexLinearForm>(b),
                           *cmfem::As<mfem::ComplexGridFunction>(x));
   }

} // extern "C"
