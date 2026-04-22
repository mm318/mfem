#ifndef CMFEM_CONSTANT_COEFFICIENT_H
#define CMFEM_CONSTANT_COEFFICIENT_H

#include "../common.h"

CMFEM_STORAGE(ConstantCoefficient, 24);

CMFEM_BEGIN_EXTERN_C

CMFEM_ConstantCoefficient *CMFEM_ConstantCoefficient_New(double value);
void CMFEM_ConstantCoefficient_Delete(CMFEM_ConstantCoefficient *coefficient);

CMFEM_END_EXTERN_C

#endif
