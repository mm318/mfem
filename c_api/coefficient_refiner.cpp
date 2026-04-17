#include "common.hpp"

#include <array>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_CoefficientRefiner, mfem::CoefficientRefiner);

using IntRuleArray = std::array<const mfem::IntegrationRule *,
      mfem::Geometry::NumGeom>;

std::mutex &IntRuleMutex()
{
   static std::mutex mutex;
   return mutex;
}

std::unordered_map<const CMFEM_CoefficientRefiner *,
    std::unique_ptr<IntRuleArray>> &IntRuleCache()
{
   static std::unordered_map<const CMFEM_CoefficientRefiner *,
          std::unique_ptr<IntRuleArray>> cache;
   return cache;
}

void ClearCachedIntRules(const CMFEM_CoefficientRefiner *refiner)
{
   std::lock_guard<std::mutex> lock(IntRuleMutex());
   IntRuleCache().erase(refiner);
}

} // namespace

extern "C" {

   CMFEM_CoefficientRefiner *CMFEM_CoefficientRefiner_NewFc(
      const CMFEM_FunctionCoefficient *coefficient,
      int order)
   {
      auto &coeff = const_cast<mfem::FunctionCoefficient &>(
                       *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      return reinterpret_cast<CMFEM_CoefficientRefiner *>(
                new mfem::CoefficientRefiner(coeff, order));
   }

   void CMFEM_CoefficientRefiner_Delete(CMFEM_CoefficientRefiner *refiner)
   {
      ClearCachedIntRules(refiner);
      delete cmfem::As<mfem::CoefficientRefiner>(refiner);
   }

   void CMFEM_CoefficientRefiner_SetThreshold(CMFEM_CoefficientRefiner *refiner,
                                              double threshold)
   {
      cmfem::As<mfem::CoefficientRefiner>(refiner)->SetThreshold(
         static_cast<mfem::real_t>(threshold));
   }

   void CMFEM_CoefficientRefiner_SetMaxElements(CMFEM_CoefficientRefiner
                                                *refiner,
                                                long long max_elements)
   {
      cmfem::As<mfem::CoefficientRefiner>(refiner)->SetMaxElements(max_elements);
   }

   void CMFEM_CoefficientRefiner_ResetCoefficientFc(
      CMFEM_CoefficientRefiner *refiner,
      const CMFEM_FunctionCoefficient *coefficient)
   {
      auto &coeff = const_cast<mfem::FunctionCoefficient &>(
                       *cmfem::As<const mfem::FunctionCoefficient>(coefficient));
      cmfem::As<mfem::CoefficientRefiner>(refiner)->ResetCoefficient(coeff);
   }

   void CMFEM_CoefficientRefiner_SetOrder(CMFEM_CoefficientRefiner *refiner,
                                          int order)
   {
      cmfem::As<mfem::CoefficientRefiner>(refiner)->SetOrder(order);
   }

   void CMFEM_CoefficientRefiner_SetNCLimit(CMFEM_CoefficientRefiner *refiner,
                                            int nc_limit)
   {
      cmfem::As<mfem::CoefficientRefiner>(refiner)->SetNCLimit(nc_limit);
   }

   void CMFEM_CoefficientRefiner_SetIntRuleOrder(
      CMFEM_CoefficientRefiner *refiner,
      int quadrature_order)
   {
      auto owned_rules = std::make_unique<IntRuleArray>();
      for (int i = 0; i < mfem::Geometry::NumGeom; i++)
      {
         (*owned_rules)[i] = &(mfem::IntRules.Get(i, quadrature_order));
      }

      auto *rules = owned_rules->data();
      {
         std::lock_guard<std::mutex> lock(IntRuleMutex());
         IntRuleCache()[refiner] = std::move(owned_rules);
      }
      cmfem::As<mfem::CoefficientRefiner>(refiner)->SetIntRule(rules);
   }

   void CMFEM_CoefficientRefiner_PrintWarnings(CMFEM_CoefficientRefiner *refiner)
   {
      cmfem::As<mfem::CoefficientRefiner>(refiner)->PrintWarnings();
   }

   int CMFEM_CoefficientRefiner_PreprocessMesh(CMFEM_CoefficientRefiner *refiner,
                                               CMFEM_Mesh *mesh)
   {
      return cmfem::As<mfem::CoefficientRefiner>(refiner)->PreprocessMesh(
                *cmfem::As<mfem::Mesh>(mesh));
   }

   double CMFEM_CoefficientRefiner_GetOsc(
      const CMFEM_CoefficientRefiner *refiner)
   {
      return cmfem::As<const mfem::CoefficientRefiner>(refiner)->GetOsc();
   }

} // extern "C"
