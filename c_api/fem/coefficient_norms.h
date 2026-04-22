#ifndef CMFEM_COEFFICIENT_NORMS_H
#define CMFEM_COEFFICIENT_NORMS_H

#include "../common.h"

CMFEM_BEGIN_EXTERN_C

double CMFEM_ComputeLpNormFcMeshOrder(double p,
                                      const CMFEM_FunctionCoefficient *coefficient,
                                      const CMFEM_Mesh *mesh,
                                      int quadrature_order);
double CMFEM_ComputeLpNormVfcMeshOrder(
   double p,
   const CMFEM_VectorFunctionCoefficient *coefficient,
   const CMFEM_Mesh *mesh,
   int quadrature_order);

CMFEM_END_EXTERN_C

#endif
