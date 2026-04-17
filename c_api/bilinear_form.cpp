#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_BilinearForm, mfem::BilinearForm);

} // namespace

extern "C" {

   CMFEM_BilinearForm *CMFEM_BilinearForm_New(CMFEM_FiniteElementSpace *fespace)
   {
      return reinterpret_cast<CMFEM_BilinearForm *>(
                new mfem::BilinearForm(cmfem::As<mfem::FiniteElementSpace>(fespace)));
   }

   void CMFEM_BilinearForm_Delete(CMFEM_BilinearForm *bilinear_form)
   {
      delete cmfem::As<mfem::BilinearForm>(bilinear_form);
   }

   void CMFEM_BilinearForm_SetAssemblyLevelPartial(CMFEM_BilinearForm
                                                   *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->SetAssemblyLevel(
         mfem::AssemblyLevel::PARTIAL);
   }

   void CMFEM_BilinearForm_SetAssemblyLevelFull(CMFEM_BilinearForm *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->SetAssemblyLevel(
         mfem::AssemblyLevel::FULL);
   }

   void CMFEM_BilinearForm_SetAssemblyLevelElement(CMFEM_BilinearForm
                                                   *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->SetAssemblyLevel(
         mfem::AssemblyLevel::ELEMENT);
   }

   void CMFEM_BilinearForm_SetDiagonalPolicyOne(CMFEM_BilinearForm *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->SetDiagonalPolicy(
         mfem::Operator::DiagonalPolicy::DIAG_ONE);
   }

   void CMFEM_BilinearForm_EnableSparseMatrixSorting(CMFEM_BilinearForm
                                                     *bilinear_form,
                                                     int enable)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->EnableSparseMatrixSorting(
         enable != 0);
   }

   void CMFEM_BilinearForm_AddDomainIntegrator_Diffusion(CMFEM_BilinearForm
                                                         *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddDomainIntegrator(
         new mfem::DiffusionIntegrator());
   }

   void CMFEM_BilinearForm_AddDomainIntegrator_DiffusionCoefficient(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddDomainIntegrator(
         new mfem::DiffusionIntegrator(coef));
   }

   void CMFEM_BilinearForm_AddDomainIntegrator_MassCoefficient(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddDomainIntegrator(
         new mfem::MassIntegrator(coef));
   }

   void CMFEM_BilinearForm_AddDomainIntegrator_DiffusionCoefficientMarker(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient,
      const CMFEM_ArrayInt *marker)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddDomainIntegrator(
         new mfem::DiffusionIntegrator(coef),
         const_cast<mfem::Array<int> &>(cmfem::ArrayIntRef(marker)));
   }

   void CMFEM_BilinearForm_AddDomainIntegrator_Elasticity(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_PWConstCoefficient *lambda,
      const CMFEM_PWConstCoefficient *mu)
   {
      auto &lambda_ref = const_cast<mfem::PWConstCoefficient &>(
                            *cmfem::As<const mfem::PWConstCoefficient>(lambda));
      auto &mu_ref = const_cast<mfem::PWConstCoefficient &>(
                        *cmfem::As<const mfem::PWConstCoefficient>(mu));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddDomainIntegrator(
         new mfem::ElasticityIntegrator(lambda_ref, mu_ref));
   }

   void CMFEM_BilinearForm_AddDomainIntegrator_CurlCurl(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddDomainIntegrator(
         new mfem::CurlCurlIntegrator(coef));
   }

   void CMFEM_BilinearForm_AddDomainIntegrator_VectorFEMass(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddDomainIntegrator(
         new mfem::VectorFEMassIntegrator(coef));
   }

   void CMFEM_BilinearForm_AddDomainIntegrator_DivDiv(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddDomainIntegrator(
         new mfem::DivDivIntegrator(coef));
   }

   void CMFEM_BilinearForm_EnableStaticCondensation(CMFEM_BilinearForm
                                                    *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->EnableStaticCondensation();
   }

   void CMFEM_BilinearForm_EnableHybridization(CMFEM_BilinearForm *bilinear_form,
                                               CMFEM_FiniteElementSpace *constraint_space,
                                               const CMFEM_ArrayInt *ess_tdof_list)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->EnableHybridization(
         cmfem::As<mfem::FiniteElementSpace>(constraint_space),
         new mfem::NormalTraceJumpIntegrator(),
         cmfem::ArrayIntRef(ess_tdof_list));
   }

   void CMFEM_BilinearForm_AddInteriorFaceIntegrator_DGDiffusionIntegrator(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient,
      double sigma,
      double kappa)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddInteriorFaceIntegrator(
         new mfem::DGDiffusionIntegrator(
            coef,
            static_cast<mfem::real_t>(sigma),
            static_cast<mfem::real_t>(kappa)));
   }

   void CMFEM_BilinearForm_AddBdrFaceIntegrator_DGDiffusionIntegrator(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ConstantCoefficient *coefficient,
      double sigma,
      double kappa)
   {
      auto &coef = const_cast<mfem::ConstantCoefficient &>(
                      *cmfem::As<const mfem::ConstantCoefficient>(coefficient));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddBdrFaceIntegrator(
         new mfem::DGDiffusionIntegrator(
            coef,
            static_cast<mfem::real_t>(sigma),
            static_cast<mfem::real_t>(kappa)));
   }

   void CMFEM_BilinearForm_AddInteriorFaceIntegrator_DGDiffusionBR2Integrator(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_FiniteElementSpace *fespace,
      double eta)
   {
      auto &fespace_ref = const_cast<mfem::FiniteElementSpace &>(
                             *cmfem::As<const mfem::FiniteElementSpace>(fespace));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddInteriorFaceIntegrator(
         new mfem::DGDiffusionBR2Integrator(
            fespace_ref,
            static_cast<mfem::real_t>(eta)));
   }

   void CMFEM_BilinearForm_AddBdrFaceIntegrator_DGDiffusionBR2Integrator(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_FiniteElementSpace *fespace,
      double eta)
   {
      auto &fespace_ref = const_cast<mfem::FiniteElementSpace &>(
                             *cmfem::As<const mfem::FiniteElementSpace>(fespace));
      cmfem::As<mfem::BilinearForm>(bilinear_form)->AddBdrFaceIntegrator(
         new mfem::DGDiffusionBR2Integrator(
            fespace_ref,
            static_cast<mfem::real_t>(eta)));
   }

   void CMFEM_BilinearForm_Assemble(CMFEM_BilinearForm *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->Assemble();
   }

   void CMFEM_BilinearForm_Finalize(CMFEM_BilinearForm *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->Finalize();
   }

   void CMFEM_BilinearForm_FormLinearSystemSparseMatrix(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ArrayInt *ess_tdof_list,
      CMFEM_GridFunction *x,
      CMFEM_LinearForm *b,
      CMFEM_SparseMatrix *A,
      CMFEM_Vector *X,
      CMFEM_Vector *B)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->FormLinearSystem(
         cmfem::ArrayIntRef(ess_tdof_list),
         *cmfem::As<mfem::GridFunction>(x),
         *cmfem::As<mfem::LinearForm>(b),
         *cmfem::As<mfem::SparseMatrix>(A),
         cmfem::VectorRef(X),
         cmfem::VectorRef(B));
   }

   void CMFEM_BilinearForm_FormLinearSystemSparseMatrixCopyInterior(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ArrayInt *ess_tdof_list,
      CMFEM_GridFunction *x,
      CMFEM_LinearForm *b,
      CMFEM_SparseMatrix *A,
      CMFEM_Vector *X,
      CMFEM_Vector *B,
      int copy_interior)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->FormLinearSystem(
         cmfem::ArrayIntRef(ess_tdof_list),
         *cmfem::As<mfem::GridFunction>(x),
         *cmfem::As<mfem::LinearForm>(b),
         *cmfem::As<mfem::SparseMatrix>(A),
         cmfem::VectorRef(X),
         cmfem::VectorRef(B),
         copy_interior);
   }

   void CMFEM_BilinearForm_FormLinearSystemOperator(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ArrayInt *ess_tdof_list,
      CMFEM_GridFunction *x,
      CMFEM_LinearForm *b,
      CMFEM_OperatorPtr *A,
      CMFEM_Vector *X,
      CMFEM_Vector *B)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->FormLinearSystem(
         cmfem::ArrayIntRef(ess_tdof_list),
         *cmfem::As<mfem::GridFunction>(x),
         *cmfem::As<mfem::LinearForm>(b),
         cmfem::OperatorPtrRef(A),
         cmfem::VectorRef(X),
         cmfem::VectorRef(B));
   }

   void CMFEM_BilinearForm_FormLinearSystemOperatorCopyInterior(
      CMFEM_BilinearForm *bilinear_form,
      const CMFEM_ArrayInt *ess_tdof_list,
      CMFEM_GridFunction *x,
      CMFEM_LinearForm *b,
      CMFEM_OperatorPtr *A,
      CMFEM_Vector *X,
      CMFEM_Vector *B,
      int copy_interior)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->FormLinearSystem(
         cmfem::ArrayIntRef(ess_tdof_list),
         *cmfem::As<mfem::GridFunction>(x),
         *cmfem::As<mfem::LinearForm>(b),
         cmfem::OperatorPtrRef(A),
         cmfem::VectorRef(X),
         cmfem::VectorRef(B),
         copy_interior);
   }

   void CMFEM_BilinearForm_RecoverFEMSolution(const CMFEM_BilinearForm
                                              *bilinear_form,
                                              const CMFEM_Vector *X,
                                              const CMFEM_LinearForm *b,
                                              CMFEM_GridFunction *x)
   {
      const_cast<mfem::BilinearForm *>(
         cmfem::As<const mfem::BilinearForm>(bilinear_form))->RecoverFEMSolution(
            cmfem::VectorRef(X),
            *cmfem::As<const mfem::LinearForm>(b),
            *cmfem::As<mfem::GridFunction>(x));
   }

   void CMFEM_BilinearForm_Update(CMFEM_BilinearForm *bilinear_form)
   {
      cmfem::As<mfem::BilinearForm>(bilinear_form)->Update();
   }

} // extern "C"
