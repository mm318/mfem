#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_MixedBilinearForm, mfem::MixedBilinearForm);

} // namespace

extern "C" {

   CMFEM_MixedBilinearForm *CMFEM_MixedBilinearForm_New(
      CMFEM_FiniteElementSpace *trial_fespace,
      CMFEM_FiniteElementSpace *test_fespace)
   {
      return reinterpret_cast<CMFEM_MixedBilinearForm *>(
                new mfem::MixedBilinearForm(
                   cmfem::As<mfem::FiniteElementSpace>(trial_fespace),
                   cmfem::As<mfem::FiniteElementSpace>(test_fespace)));
   }

   void CMFEM_MixedBilinearForm_Delete(
      CMFEM_MixedBilinearForm *mixed_bilinear_form)
   {
      delete cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form);
   }

   void CMFEM_MixedBilinearForm_SetAssemblyLevelPartial(
      CMFEM_MixedBilinearForm *mixed_bilinear_form)
   {
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->SetAssemblyLevel(
         mfem::AssemblyLevel::PARTIAL);
   }

   void CMFEM_MixedBilinearForm_AddDomainIntegratorMvgiCc(
      CMFEM_MixedBilinearForm *mixed_bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->AddDomainIntegrator(
         new mfem::MixedVectorGradientIntegrator(coef));
   }

   void CMFEM_MixedBilinearForm_AddDomainIntegratorDiCc(
      CMFEM_MixedBilinearForm *mixed_bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->AddDomainIntegrator(
         new mfem::DiffusionIntegrator(coef));
   }

   void CMFEM_MixedBilinearForm_AddDomainIntegratorMvciCc(
      CMFEM_MixedBilinearForm *mixed_bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->AddDomainIntegrator(
         new mfem::MixedVectorCurlIntegrator(coef));
   }

   void CMFEM_MixedBilinearForm_AddDomainIntegratorVfdiCc(
      CMFEM_MixedBilinearForm *mixed_bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->AddDomainIntegrator(
         new mfem::VectorFEDivergenceIntegrator(coef));
   }

   void CMFEM_MixedBilinearForm_AddDomainIntegratorMsi(
      CMFEM_MixedBilinearForm *mixed_bilinear_form)
   {
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->AddDomainIntegrator(
         new mfem::MixedScalarMassIntegrator());
   }

   void CMFEM_MixedBilinearForm_AddTraceFaceIntegratorTji(
      CMFEM_MixedBilinearForm *mixed_bilinear_form)
   {
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)
      ->AddTraceFaceIntegrator(new mfem::TraceJumpIntegrator());
   }

   void CMFEM_MixedBilinearForm_Assemble(
      CMFEM_MixedBilinearForm *mixed_bilinear_form)
   {
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->Assemble();
   }

   void CMFEM_MixedBilinearForm_Finalize(
      CMFEM_MixedBilinearForm *mixed_bilinear_form)
   {
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->Finalize();
   }

   void CMFEM_MixedBilinearForm_EliminateTrialEssentialBCAiVecLf(
      CMFEM_MixedBilinearForm *mixed_bilinear_form,
      const CMFEM_ArrayInt *essential_bdr,
      const CMFEM_Vector *solution,
      CMFEM_LinearForm *rhs)
   {
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)
      ->EliminateTrialEssentialBC(cmfem::ArrayIntRef(essential_bdr),
                                  cmfem::VectorRef(solution),
                                  *cmfem::As<mfem::LinearForm>(rhs));
   }

   void CMFEM_MixedBilinearForm_EliminateTrialEssentialBCAiVecVec(
      CMFEM_MixedBilinearForm *mixed_bilinear_form,
      const CMFEM_ArrayInt *essential_bdr,
      const CMFEM_Vector *solution,
      CMFEM_Vector *rhs)
   {
      cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)
      ->EliminateTrialEssentialBC(cmfem::ArrayIntRef(essential_bdr),
                                  cmfem::VectorRef(solution),
                                  cmfem::VectorRef(rhs));
   }

   CMFEM_SparseMatrix *CMFEM_MixedBilinearForm_SpMat(
      CMFEM_MixedBilinearForm *mixed_bilinear_form)
   {
      return reinterpret_cast<CMFEM_SparseMatrix *>(
                &cmfem::As<mfem::MixedBilinearForm>(mixed_bilinear_form)->SpMat());
   }

   void CMFEM_MixedBilinearForm_Mult(
      const CMFEM_MixedBilinearForm *mixed_bilinear_form,
      const CMFEM_Vector *x,
      CMFEM_Vector *y)
   {
      cmfem::As<const mfem::MixedBilinearForm>(mixed_bilinear_form)->Mult(
         cmfem::VectorRef(x),
         cmfem::VectorRef(y));
   }

} // extern "C"
