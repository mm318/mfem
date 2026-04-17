#ifndef CMFEM_PW_CONST_COEFFICIENT_H
#define CMFEM_PW_CONST_COEFFICIENT_H

#include "common.h"

CMFEM_STORAGE(PWConstCoefficient, 56);

CMFEM_BEGIN_EXTERN_C

CMFEM_PWConstCoefficient *CMFEM_PWConstCoefficient_New(const CMFEM_Vector *values);
void CMFEM_PWConstCoefficient_Delete(CMFEM_PWConstCoefficient *coefficient);

CMFEM_END_EXTERN_C

#endif
