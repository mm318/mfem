#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ThresholdDerefiner, mfem::ThresholdDerefiner);

} // namespace

extern "C" {

   CMFEM_ThresholdDerefiner *
   CMFEM_ThresholdDerefiner_NewZze(
      CMFEM_ZienkiewiczZhuEstimator *estimator)
   {
      return reinterpret_cast<CMFEM_ThresholdDerefiner *>(
                new mfem::ThresholdDerefiner(
                   *cmfem::As<mfem::ZienkiewiczZhuEstimator>(estimator)));
   }

   CMFEM_ThresholdDerefiner *
   CMFEM_ThresholdDerefiner_NewKee(
      CMFEM_KellyErrorEstimator *estimator)
   {
      return reinterpret_cast<CMFEM_ThresholdDerefiner *>(
                new mfem::ThresholdDerefiner(
                   *cmfem::As<mfem::KellyErrorEstimator>(estimator)));
   }

   void CMFEM_ThresholdDerefiner_Delete(CMFEM_ThresholdDerefiner *derefiner)
   {
      delete cmfem::As<mfem::ThresholdDerefiner>(derefiner);
   }

   void CMFEM_ThresholdDerefiner_SetThreshold(CMFEM_ThresholdDerefiner *derefiner,
                                              double threshold)
   {
      cmfem::As<mfem::ThresholdDerefiner>(derefiner)->SetThreshold(
         static_cast<mfem::real_t>(threshold));
   }

   void CMFEM_ThresholdDerefiner_SetNCLimit(CMFEM_ThresholdDerefiner *derefiner,
                                            int nc_limit)
   {
      cmfem::As<mfem::ThresholdDerefiner>(derefiner)->SetNCLimit(nc_limit);
   }

   void CMFEM_ThresholdDerefiner_Reset(CMFEM_ThresholdDerefiner *derefiner)
   {
      cmfem::As<mfem::ThresholdDerefiner>(derefiner)->Reset();
   }

   int CMFEM_ThresholdDerefiner_Apply(CMFEM_ThresholdDerefiner *derefiner,
                                      CMFEM_Mesh *mesh)
   {
      return cmfem::As<mfem::ThresholdDerefiner>(derefiner)->Apply(
                *cmfem::As<mfem::Mesh>(mesh)) ? 1 : 0;
   }

} // extern "C"
