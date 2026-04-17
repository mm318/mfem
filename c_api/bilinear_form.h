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
void CMFEM_BilinearForm_AddDomainIntegratorDi(CMFEM_BilinearForm
                                              *bilinear_form);
void CMFEM_BilinearForm_AddDomainIntegratorDiCc(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_AddDomainIntegratorMiCc(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_AddDomainIntegratorDiCcMarker(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient,
   const CMFEM_ArrayInt *marker);
void CMFEM_BilinearForm_AddDomainIntegratorEi(CMFEM_BilinearForm
                                              *bilinear_form, const CMFEM_PWConstCoefficient *lambda,
                                              const CMFEM_PWConstCoefficient *mu);
void CMFEM_BilinearForm_AddDomainIntegratorCci(CMFEM_BilinearForm
                                               *bilinear_form, const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_AddDomainIntegratorVmi(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_AddDomainIntegratorDdi(CMFEM_BilinearForm
                                               *bilinear_form, const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_BilinearForm_EnableStaticCondensation(CMFEM_BilinearForm
                                                 *bilinear_form);
void CMFEM_BilinearForm_EnableHybridization(CMFEM_BilinearForm *bilinear_form,
                                            CMFEM_FiniteElementSpace *constraint_space,
                                            const CMFEM_ArrayInt *ess_tdof_list);
void CMFEM_BilinearForm_AddInteriorFaceIntegratorDgd(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient,
   double sigma,
   double kappa);
void CMFEM_BilinearForm_AddBdrFaceIntegratorDgd(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient,
   double sigma,
   double kappa);
void CMFEM_BilinearForm_AddInteriorFaceIntegratorDgb(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_FiniteElementSpace *fespace,
   double eta);
void CMFEM_BilinearForm_AddBdrFaceIntegratorDgb(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_FiniteElementSpace *fespace,
   double eta);
void CMFEM_BilinearForm_Assemble(CMFEM_BilinearForm *bilinear_form);
void CMFEM_BilinearForm_Finalize(CMFEM_BilinearForm *bilinear_form);
void CMFEM_BilinearForm_FormLinearSystemSm(CMFEM_BilinearForm
                                           *bilinear_form, const CMFEM_ArrayInt *ess_tdof_list, CMFEM_GridFunction *x,
                                           CMFEM_LinearForm *b, CMFEM_SparseMatrix *A, CMFEM_Vector *X, CMFEM_Vector *B);
void CMFEM_BilinearForm_FormLinearSystemSmCopyInterior(
   CMFEM_BilinearForm *bilinear_form,
   const CMFEM_ArrayInt *ess_tdof_list,
   CMFEM_GridFunction *x,
   CMFEM_LinearForm *b,
   CMFEM_SparseMatrix *A,
   CMFEM_Vector *X,
   CMFEM_Vector *B,
   int copy_interior);
void CMFEM_BilinearForm_FormLinearSystemOp(CMFEM_BilinearForm
                                           *bilinear_form, const CMFEM_ArrayInt *ess_tdof_list, CMFEM_GridFunction *x,
                                           CMFEM_LinearForm *b, CMFEM_OperatorPtr *A, CMFEM_Vector *X, CMFEM_Vector *B);
void CMFEM_BilinearForm_FormLinearSystemOpCopyInterior(
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
