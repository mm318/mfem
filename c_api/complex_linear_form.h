#ifndef CMFEM_COMPLEX_LINEAR_FORM_H
#define CMFEM_COMPLEX_LINEAR_FORM_H

#include "common.h"

CMFEM_STORAGE(ComplexLinearForm, 56);

CMFEM_BEGIN_EXTERN_C

CMFEM_ComplexLinearForm *CMFEM_ComplexLinearForm_NewFesConv(
   CMFEM_FiniteElementSpace *fespace,
   enum CMFEM_ComplexConvention convention);
void CMFEM_ComplexLinearForm_Delete(CMFEM_ComplexLinearForm *linear_form);
void CMFEM_ComplexLinearForm_Assign(CMFEM_ComplexLinearForm *linear_form,
                                    double real_value,
                                    double imag_value);
void CMFEM_ComplexLinearForm_AddDomainIntegratorNullVfdVfc(
   CMFEM_ComplexLinearForm *linear_form,
   const CMFEM_VectorFunctionCoefficient *coefficient);
void CMFEM_ComplexLinearForm_Assemble(CMFEM_ComplexLinearForm *linear_form);

CMFEM_END_EXTERN_C

#endif
