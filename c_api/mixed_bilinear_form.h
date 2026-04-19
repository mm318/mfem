#ifndef CMFEM_MIXED_BILINEAR_FORM_H
#define CMFEM_MIXED_BILINEAR_FORM_H

#include "common.h"

CMFEM_STORAGE(MixedBilinearForm, 504);

CMFEM_BEGIN_EXTERN_C

CMFEM_MixedBilinearForm *CMFEM_MixedBilinearForm_New(
   CMFEM_FiniteElementSpace *trial_fespace,
   CMFEM_FiniteElementSpace *test_fespace);
void CMFEM_MixedBilinearForm_Delete(CMFEM_MixedBilinearForm
                                    *mixed_bilinear_form);
void CMFEM_MixedBilinearForm_SetAssemblyLevelPartial(
   CMFEM_MixedBilinearForm *mixed_bilinear_form);
void CMFEM_MixedBilinearForm_AddDomainIntegratorMvgiCc(
   CMFEM_MixedBilinearForm *mixed_bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_MixedBilinearForm_AddDomainIntegratorDiCc(
   CMFEM_MixedBilinearForm *mixed_bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_MixedBilinearForm_AddDomainIntegratorMvciCc(
   CMFEM_MixedBilinearForm *mixed_bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_MixedBilinearForm_AddDomainIntegratorVfdiCc(
   CMFEM_MixedBilinearForm *mixed_bilinear_form,
   const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_MixedBilinearForm_AddDomainIntegratorMsi(
   CMFEM_MixedBilinearForm *mixed_bilinear_form);
void CMFEM_MixedBilinearForm_AddTraceFaceIntegratorTji(
   CMFEM_MixedBilinearForm *mixed_bilinear_form);
void CMFEM_MixedBilinearForm_Assemble(CMFEM_MixedBilinearForm
                                      *mixed_bilinear_form);
void CMFEM_MixedBilinearForm_Finalize(CMFEM_MixedBilinearForm
                                      *mixed_bilinear_form);
void CMFEM_MixedBilinearForm_EliminateTrialEssentialBCAiVecLf(
   CMFEM_MixedBilinearForm *mixed_bilinear_form,
   const CMFEM_ArrayInt *essential_bdr,
   const CMFEM_Vector *solution,
   CMFEM_LinearForm *rhs);
void CMFEM_MixedBilinearForm_EliminateTrialEssentialBCAiVecVec(
   CMFEM_MixedBilinearForm *mixed_bilinear_form,
   const CMFEM_ArrayInt *essential_bdr,
   const CMFEM_Vector *solution,
   CMFEM_Vector *rhs);
CMFEM_SparseMatrix *CMFEM_MixedBilinearForm_SpMat(
   CMFEM_MixedBilinearForm *mixed_bilinear_form);
void CMFEM_MixedBilinearForm_Mult(
   const CMFEM_MixedBilinearForm *mixed_bilinear_form,
   const CMFEM_Vector *x,
   CMFEM_Vector *y);

CMFEM_END_EXTERN_C

#endif
