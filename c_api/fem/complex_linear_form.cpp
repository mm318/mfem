#include "common.hpp"

#include <complex>

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ComplexLinearForm, mfem::ComplexLinearForm);

mfem::ComplexOperator::Convention ToConvention(
   enum CMFEM_ComplexConvention convention)
{
   return convention == CMFEM_ComplexConventionBlockSymmetric ?
          mfem::ComplexOperator::BLOCK_SYMMETRIC :
          mfem::ComplexOperator::HERMITIAN;
}

} // namespace

extern "C" {

   CMFEM_ComplexLinearForm *CMFEM_ComplexLinearForm_NewFesConv(
      CMFEM_FiniteElementSpace *fespace,
      enum CMFEM_ComplexConvention convention)
   {
      return reinterpret_cast<CMFEM_ComplexLinearForm *>(
                new mfem::ComplexLinearForm(
                   cmfem::As<mfem::FiniteElementSpace>(fespace),
                   ToConvention(convention)));
   }

   void CMFEM_ComplexLinearForm_Delete(CMFEM_ComplexLinearForm *linear_form)
   {
      delete cmfem::As<mfem::ComplexLinearForm>(linear_form);
   }

   void CMFEM_ComplexLinearForm_Assign(CMFEM_ComplexLinearForm *linear_form,
                                       double real_value,
                                       double imag_value)
   {
      *cmfem::As<mfem::ComplexLinearForm>(linear_form) =
         std::complex<mfem::real_t>(static_cast<mfem::real_t>(real_value),
                                    static_cast<mfem::real_t>(imag_value));
   }

   void CMFEM_ComplexLinearForm_AddDomainIntegratorNullVfdVfc(
      CMFEM_ComplexLinearForm *linear_form,
      const CMFEM_VectorFunctionCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::VectorFunctionCoefficient &>(
                      *cmfem::As<const mfem::VectorFunctionCoefficient>(coefficient));
      cmfem::As<mfem::ComplexLinearForm>(linear_form)->AddDomainIntegrator(
         nullptr,
         new mfem::VectorFEDomainLFIntegrator(coef));
   }

   void CMFEM_ComplexLinearForm_Assemble(CMFEM_ComplexLinearForm *linear_form)
   {
      cmfem::As<mfem::ComplexLinearForm>(linear_form)->Assemble();
   }

} // extern "C"
