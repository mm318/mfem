#ifndef CMFEM_SESQUILINEAR_FORM_H
#define CMFEM_SESQUILINEAR_FORM_H

#include "../common.h"

CMFEM_STORAGE(SesquilinearForm, 32);

CMFEM_BEGIN_EXTERN_C

CMFEM_SesquilinearForm *CMFEM_SesquilinearForm_NewFesConv(
   CMFEM_FiniteElementSpace *fespace,
   enum CMFEM_ComplexConvention convention);
void CMFEM_SesquilinearForm_Delete(CMFEM_SesquilinearForm *sesquilinear_form);
void CMFEM_SesquilinearForm_SetAssemblyLevelPartial(
   CMFEM_SesquilinearForm *sesquilinear_form);
void CMFEM_SesquilinearForm_AddDomainIntegratorDiCcNull(
   CMFEM_SesquilinearForm *sesquilinear_form,
   const CMFEM_ConstantCoefficient *real_coefficient);
void CMFEM_SesquilinearForm_AddDomainIntegratorMiCcMiCc(
   CMFEM_SesquilinearForm *sesquilinear_form,
   const CMFEM_ConstantCoefficient *real_coefficient,
   const CMFEM_ConstantCoefficient *imag_coefficient);
void CMFEM_SesquilinearForm_AddDomainIntegratorCciCcNull(
   CMFEM_SesquilinearForm *sesquilinear_form,
   const CMFEM_ConstantCoefficient *real_coefficient);
void CMFEM_SesquilinearForm_AddDomainIntegratorVmiCcVmiCc(
   CMFEM_SesquilinearForm *sesquilinear_form,
   const CMFEM_ConstantCoefficient *real_coefficient,
   const CMFEM_ConstantCoefficient *imag_coefficient);
void CMFEM_SesquilinearForm_AddDomainIntegratorDdiCcNull(
   CMFEM_SesquilinearForm *sesquilinear_form,
   const CMFEM_ConstantCoefficient *real_coefficient);
void CMFEM_SesquilinearForm_Assemble(CMFEM_SesquilinearForm *sesquilinear_form);
void CMFEM_SesquilinearForm_AssembleSkipZeros(
   CMFEM_SesquilinearForm *sesquilinear_form,
   int skip_zeros);
void CMFEM_SesquilinearForm_FormLinearSystemOp(
   CMFEM_SesquilinearForm *sesquilinear_form,
   const CMFEM_ArrayInt *ess_tdof_list,
   CMFEM_ComplexGridFunction *x,
   CMFEM_ComplexLinearForm *b,
   CMFEM_OperatorPtr *A,
   CMFEM_Vector *X,
   CMFEM_Vector *B);
void CMFEM_SesquilinearForm_RecoverFEMSolution(
   const CMFEM_SesquilinearForm *sesquilinear_form,
   const CMFEM_Vector *X,
   const CMFEM_ComplexLinearForm *b,
   CMFEM_ComplexGridFunction *x);

CMFEM_END_EXTERN_C

#endif
