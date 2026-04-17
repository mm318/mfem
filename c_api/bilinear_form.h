#ifndef CMFEM_BILINEAR_FORM_H
#define CMFEM_BILINEAR_FORM_H

#include "common.h"

CMFEM_STORAGE(BilinearForm, 416);

CMFEM_BEGIN_EXTERN_C

CMFEM_BilinearForm *CMFEM_BilinearForm_New(CMFEM_FiniteElementSpace *fespace);
void CMFEM_BilinearForm_Delete(CMFEM_BilinearForm *bilinear_form);
void CMFEM_BilinearForm_SetAssemblyLevelPartial(CMFEM_BilinearForm
                                                *bilinear_form);
void CMFEM_BilinearForm_SetAssemblyLevelFull(CMFEM_BilinearForm *bilinear_form);
void CMFEM_BilinearForm_SetAssemblyLevelElement(CMFEM_BilinearForm
                                                *bilinear_form);
void CMFEM_BilinearForm_SetDiagonalPolicyOne(CMFEM_BilinearForm *bilinear_form);
void CMFEM_BilinearForm_EnableSparseMatrixSorting(CMFEM_BilinearForm
                                                  *bilinear_form, int enable);
void CMFEM_BilinearForm_AddDomainIntegrator_Diffusion(CMFEM_BilinearForm
                                                      *bilinear_form);
void CMFEM_BilinearForm_AddDomainIntegrator_DiffusionCoefficient(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_AddDomainIntegrator_MassCoefficient(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_AddDomainIntegrator_DiffusionCoefficientMarker(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient,
   const CMFEM_ArrayInt *marker);
void CMFEM_BilinearForm_AddDomainIntegrator_Elasticity(CMFEM_BilinearForm
                                                       *bilinear_form, const CMFEM_PWConstCoefficient *lambda,
                                                       const CMFEM_PWConstCoefficient *mu);
void CMFEM_BilinearForm_AddDomainIntegrator_CurlCurl(CMFEM_BilinearForm
                                                     *bilinear_form, const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_AddDomainIntegrator_VectorFEMass(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_AddDomainIntegrator_DivDiv(CMFEM_BilinearForm
                                                   *bilinear_form, const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_EnableStaticCondensation(CMFEM_BilinearForm
                                                 *bilinear_form);
void CMFEM_BilinearForm_EnableHybridization(CMFEM_BilinearForm *bilinear_form,
                                            CMFEM_FiniteElementSpace *constraint_space,
                                            const CMFEM_ArrayInt *ess_tdof_list);
void CMFEM_BilinearForm_AddInteriorFaceIntegrator_DGDiffusionIntegrator(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient,
   double sigma,
   double kappa);
void CMFEM_BilinearForm_AddBdrFaceIntegrator_DGDiffusionIntegrator(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient,
   double sigma,
   double kappa);
void CMFEM_BilinearForm_AddInteriorFaceIntegrator_DGDiffusionBR2Integrator(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_FiniteElementSpace *fespace,
   double eta);
void CMFEM_BilinearForm_AddBdrFaceIntegrator_DGDiffusionBR2Integrator(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_FiniteElementSpace *fespace,
   double eta);
void CMFEM_BilinearForm_Assemble(CMFEM_BilinearForm *bilinear_form);
void CMFEM_BilinearForm_Finalize(CMFEM_BilinearForm *bilinear_form);
void CMFEM_BilinearForm_FormLinearSystemSparseMatrix(CMFEM_BilinearForm
                                                     *bilinear_form, const CMFEM_ArrayInt *ess_tdof_list, CMFEM_GridFunction *x,
                                                     CMFEM_LinearForm *b, CMFEM_SparseMatrix *A, CMFEM_Vector *X, CMFEM_Vector *B);
void CMFEM_BilinearForm_FormLinearSystemSparseMatrixCopyInterior(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ArrayInt *ess_tdof_list,
   CMFEM_GridFunction *x,
   CMFEM_LinearForm *b,
   CMFEM_SparseMatrix *A,
   CMFEM_Vector *X,
   CMFEM_Vector *B,
   int copy_interior);
void CMFEM_BilinearForm_FormLinearSystemOperator(CMFEM_BilinearForm
                                                 *bilinear_form, const CMFEM_ArrayInt *ess_tdof_list, CMFEM_GridFunction *x,
                                                 CMFEM_LinearForm *b, CMFEM_OperatorPtr *A, CMFEM_Vector *X, CMFEM_Vector *B);
void CMFEM_BilinearForm_FormLinearSystemOperatorCopyInterior(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ArrayInt *ess_tdof_list,
   CMFEM_GridFunction *x,
   CMFEM_LinearForm *b,
   CMFEM_OperatorPtr *A,
   CMFEM_Vector *X,
   CMFEM_Vector *B,
   int copy_interior);
void CMFEM_BilinearForm_RecoverFEMSolution(const CMFEM_BilinearForm
                                           *bilinear_form, const CMFEM_Vector *X, const CMFEM_LinearForm *b,
                                           CMFEM_GridFunction *x);
void CMFEM_BilinearForm_Update(CMFEM_BilinearForm *bilinear_form);

CMFEM_END_EXTERN_C

#endif
