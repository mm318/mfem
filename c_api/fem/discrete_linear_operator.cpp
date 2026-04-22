#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_DiscreteLinearOperator, mfem::DiscreteLinearOperator);

} // namespace

extern "C" {

   CMFEM_DiscreteLinearOperator *CMFEM_DiscreteLinearOperator_New(
      CMFEM_FiniteElementSpace *domain_fespace,
      CMFEM_FiniteElementSpace *range_fespace)
   {
      return reinterpret_cast<CMFEM_DiscreteLinearOperator *>(
                new mfem::DiscreteLinearOperator(
                   cmfem::As<mfem::FiniteElementSpace>(domain_fespace),
                   cmfem::As<mfem::FiniteElementSpace>(range_fespace)));
   }

   void CMFEM_DiscreteLinearOperator_Delete(
      CMFEM_DiscreteLinearOperator *discrete_linear_operator)
   {
      delete cmfem::As<mfem::DiscreteLinearOperator>(discrete_linear_operator);
   }

   void CMFEM_DiscreteLinearOperator_AddDomainInterpolatorGi(
      CMFEM_DiscreteLinearOperator *discrete_linear_operator)
   {
      cmfem::As<mfem::DiscreteLinearOperator>(discrete_linear_operator)
      ->AddDomainInterpolator(new mfem::GradientInterpolator());
   }

   void CMFEM_DiscreteLinearOperator_AddDomainInterpolatorCi(
      CMFEM_DiscreteLinearOperator *discrete_linear_operator)
   {
      cmfem::As<mfem::DiscreteLinearOperator>(discrete_linear_operator)
      ->AddDomainInterpolator(new mfem::CurlInterpolator());
   }

   void CMFEM_DiscreteLinearOperator_AddDomainInterpolatorDi(
      CMFEM_DiscreteLinearOperator *discrete_linear_operator)
   {
      cmfem::As<mfem::DiscreteLinearOperator>(discrete_linear_operator)
      ->AddDomainInterpolator(new mfem::DivergenceInterpolator());
   }

   void CMFEM_DiscreteLinearOperator_Assemble(
      CMFEM_DiscreteLinearOperator *discrete_linear_operator)
   {
      cmfem::As<mfem::DiscreteLinearOperator>(discrete_linear_operator)->Assemble();
   }

   void CMFEM_DiscreteLinearOperator_Mult(
      const CMFEM_DiscreteLinearOperator *discrete_linear_operator,
      const CMFEM_Vector *x,
      CMFEM_Vector *y)
   {
      cmfem::As<const mfem::DiscreteLinearOperator>(discrete_linear_operator)->Mult(
         cmfem::VectorRef(x),
         cmfem::VectorRef(y));
   }

} // extern "C"
