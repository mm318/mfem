#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ThresholdRefiner, mfem::ThresholdRefiner);

} // namespace

extern "C" {

   CMFEM_ThresholdRefiner *
   CMFEM_ThresholdRefiner_NewZze(
      CMFEM_ZienkiewiczZhuEstimator *estimator)
   {
      return reinterpret_cast<CMFEM_ThresholdRefiner *>(
                new mfem::ThresholdRefiner(
                   *cmfem::As<mfem::ZienkiewiczZhuEstimator>(estimator)));
   }

   CMFEM_ThresholdRefiner *
   CMFEM_ThresholdRefiner_NewLzz(
      CMFEM_LSZienkiewiczZhuEstimator *estimator)
   {
      return reinterpret_cast<CMFEM_ThresholdRefiner *>(
                new mfem::ThresholdRefiner(
                   *cmfem::As<mfem::LSZienkiewiczZhuEstimator>(estimator)));
   }

   CMFEM_ThresholdRefiner *
   CMFEM_ThresholdRefiner_NewKee(
      CMFEM_KellyErrorEstimator *estimator)
   {
      return reinterpret_cast<CMFEM_ThresholdRefiner *>(
                new mfem::ThresholdRefiner(
                   *cmfem::As<mfem::KellyErrorEstimator>(estimator)));
   }

   void CMFEM_ThresholdRefiner_Delete(CMFEM_ThresholdRefiner *refiner)
   {
      delete cmfem::As<mfem::ThresholdRefiner>(refiner);
   }

   void CMFEM_ThresholdRefiner_SetTotalErrorFraction(CMFEM_ThresholdRefiner
                                                     *refiner,
                                                     double fraction)
   {
      cmfem::As<mfem::ThresholdRefiner>(refiner)->SetTotalErrorFraction(
         static_cast<mfem::real_t>(fraction));
   }

   void CMFEM_ThresholdRefiner_SetLocalErrorGoal(CMFEM_ThresholdRefiner *refiner,
                                                 double goal)
   {
      cmfem::As<mfem::ThresholdRefiner>(refiner)->SetLocalErrorGoal(
         static_cast<mfem::real_t>(goal));
   }

   void CMFEM_ThresholdRefiner_PreferConformingRefinement(
      CMFEM_ThresholdRefiner *refiner)
   {
      cmfem::As<mfem::ThresholdRefiner>(refiner)->PreferConformingRefinement();
   }

   void CMFEM_ThresholdRefiner_SetNCLimit(CMFEM_ThresholdRefiner *refiner,
                                          int nc_limit)
   {
      cmfem::As<mfem::ThresholdRefiner>(refiner)->SetNCLimit(nc_limit);
   }

   void CMFEM_ThresholdRefiner_Reset(CMFEM_ThresholdRefiner *refiner)
   {
      cmfem::As<mfem::ThresholdRefiner>(refiner)->Reset();
   }

   int CMFEM_ThresholdRefiner_Apply(CMFEM_ThresholdRefiner *refiner,
                                    CMFEM_Mesh *mesh)
   {
      return cmfem::As<mfem::ThresholdRefiner>(refiner)->Apply(
                *cmfem::As<mfem::Mesh>(mesh)) ? 1 : 0;
   }

   int CMFEM_ThresholdRefiner_Stop(const CMFEM_ThresholdRefiner *refiner)
   {
      return cmfem::As<const mfem::ThresholdRefiner>(refiner)->Stop() ? 1 : 0;
   }

} // extern "C"
