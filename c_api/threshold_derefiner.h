#ifndef CMFEM_THRESHOLD_DEREFINER_H
#define CMFEM_THRESHOLD_DEREFINER_H

#include "common.h"

CMFEM_STORAGE(ThresholdDerefiner, 40);

CMFEM_BEGIN_EXTERN_C

CMFEM_ThresholdDerefiner *
CMFEM_ThresholdDerefiner_NewZienkiewiczZhuEstimator(
   CMFEM_ZienkiewiczZhuEstimator *estimator);
CMFEM_ThresholdDerefiner *
CMFEM_ThresholdDerefiner_NewKellyErrorEstimator(
   CMFEM_KellyErrorEstimator *estimator);
void CMFEM_ThresholdDerefiner_Delete(CMFEM_ThresholdDerefiner *derefiner);
void CMFEM_ThresholdDerefiner_SetThreshold(CMFEM_ThresholdDerefiner *derefiner,
                                           double threshold);
void CMFEM_ThresholdDerefiner_SetNCLimit(CMFEM_ThresholdDerefiner *derefiner,
                                         int nc_limit);
void CMFEM_ThresholdDerefiner_Reset(CMFEM_ThresholdDerefiner *derefiner);
int CMFEM_ThresholdDerefiner_Apply(CMFEM_ThresholdDerefiner *derefiner,
                                   CMFEM_Mesh *mesh);

CMFEM_END_EXTERN_C

#endif
