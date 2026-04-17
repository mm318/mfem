#ifndef CMFEM_COEFFICIENT_REFINER_H
#define CMFEM_COEFFICIENT_REFINER_H

#include "common.h"

CMFEM_STORAGE(CoefficientRefiner, 208);

CMFEM_BEGIN_EXTERN_C

CMFEM_CoefficientRefiner *CMFEM_CoefficientRefiner_NewFc(
   const CMFEM_FunctionCoefficient *coefficient,
   int order);
void CMFEM_CoefficientRefiner_Delete(CMFEM_CoefficientRefiner *refiner);
void CMFEM_CoefficientRefiner_SetThreshold(CMFEM_CoefficientRefiner *refiner,
                                           double threshold);
void CMFEM_CoefficientRefiner_SetMaxElements(CMFEM_CoefficientRefiner *refiner,
                                             long long max_elements);
void CMFEM_CoefficientRefiner_ResetCoefficientFc(
   CMFEM_CoefficientRefiner *refiner,
   const CMFEM_FunctionCoefficient *coefficient);
void CMFEM_CoefficientRefiner_SetOrder(CMFEM_CoefficientRefiner *refiner,
                                       int order);
void CMFEM_CoefficientRefiner_SetNCLimit(CMFEM_CoefficientRefiner *refiner,
                                         int nc_limit);
void CMFEM_CoefficientRefiner_SetIntRuleOrder(
   CMFEM_CoefficientRefiner *refiner,
   int quadrature_order);
void CMFEM_CoefficientRefiner_PrintWarnings(CMFEM_CoefficientRefiner *refiner);
int CMFEM_CoefficientRefiner_PreprocessMesh(CMFEM_CoefficientRefiner *refiner,
                                            CMFEM_Mesh *mesh);
double CMFEM_CoefficientRefiner_GetOsc(
   const CMFEM_CoefficientRefiner *refiner);

CMFEM_END_EXTERN_C

#endif
