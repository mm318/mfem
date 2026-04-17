#ifndef CMFEM_LINEAR_FORM_H
#define CMFEM_LINEAR_FORM_H

#include "common.h"

CMFEM_STORAGE(LinearForm, 384);

CMFEM_BEGIN_EXTERN_C

CMFEM_LinearForm *CMFEM_LinearForm_New(CMFEM_FiniteElementSpace *fespace);
void CMFEM_LinearForm_Delete(CMFEM_LinearForm *linear_form);
void CMFEM_LinearForm_AddDomainIntegratorDliCc(
   CMFEM_LinearForm *linear_form, const CMFEM_ConstantCoefficient *coefficient);
void CMFEM_LinearForm_AddDomainIntegratorDliFc(
   CMFEM_LinearForm *linear_form,
   const CMFEM_FunctionCoefficient *coefficient);
void
CMFEM_LinearForm_AddDomainIntegratorDliCcMarker(
   CMFEM_LinearForm *linear_form,
   const CMFEM_ConstantCoefficient *coefficient,
   const CMFEM_ArrayInt *marker);
void CMFEM_LinearForm_AddDomainIntegratorVfd(
   CMFEM_LinearForm *linear_form,
   const CMFEM_VectorFunctionCoefficient *coefficient);
void CMFEM_LinearForm_AddBdrFaceIntegratorDgl(
   CMFEM_LinearForm *linear_form,
   const CMFEM_ConstantCoefficient *u_coeff,
   const CMFEM_ConstantCoefficient *q_coeff,
   double sigma,
   double kappa);
void CMFEM_LinearForm_AddBoundaryIntegratorVbl(
   CMFEM_LinearForm *linear_form, const CMFEM_VectorArrayCoefficient *coefficient);
void CMFEM_LinearForm_Assemble(CMFEM_LinearForm *linear_form);
void CMFEM_LinearForm_Update(CMFEM_LinearForm *linear_form);

CMFEM_END_EXTERN_C

#endif
