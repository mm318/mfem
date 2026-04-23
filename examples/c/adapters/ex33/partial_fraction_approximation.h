#ifndef CMFEM_PARTIAL_FRACTION_APPROXIMATION_H
#define CMFEM_PARTIAL_FRACTION_APPROXIMATION_H

#include "../../../../c_api/common.h"

CMFEM_BEGIN_EXTERN_C

void CMFEM_ComputePartialFractionApproximation(double *alpha,
                                               int *size,
                                               double **coeffs,
                                               double **poles);
void CMFEM_FreeDoubles(double *values);

CMFEM_END_EXTERN_C

#endif
