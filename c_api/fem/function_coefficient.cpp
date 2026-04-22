#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_FunctionCoefficient, mfem::FunctionCoefficient);

} // namespace

extern "C" {

   CMFEM_FunctionCoefficient *CMFEM_FunctionCoefficient_New(
      CMFEM_FunctionCallback callback,
      void *context)
   {
      auto function = [callback, context](const mfem::Vector &x)
      {
         return static_cast<mfem::real_t>(
                   callback(reinterpret_cast<const CMFEM_Vector *>(&x), context));
      };
      return reinterpret_cast<CMFEM_FunctionCoefficient *>(
                new mfem::FunctionCoefficient(std::move(function)));
   }

   CMFEM_FunctionCoefficient *CMFEM_FunctionCoefficient_NewTimeDependent(
      CMFEM_TimeDependentFunctionCallback callback,
      void *context)
   {
      auto function = [callback, context](const mfem::Vector &x, mfem::real_t time)
      {
         return static_cast<mfem::real_t>(
                   callback(reinterpret_cast<const CMFEM_Vector *>(&x), time,
                            context));
      };
      return reinterpret_cast<CMFEM_FunctionCoefficient *>(
                new mfem::FunctionCoefficient(std::move(function)));
   }

   void CMFEM_FunctionCoefficient_SetTime(CMFEM_FunctionCoefficient *coefficient,
                                          double time)
   {
      cmfem::As<mfem::FunctionCoefficient>(coefficient)->SetTime(
         static_cast<mfem::real_t>(time));
   }

   void CMFEM_FunctionCoefficient_Delete(CMFEM_FunctionCoefficient *coefficient)
   {
      delete cmfem::As<mfem::FunctionCoefficient>(coefficient);
   }

} // extern "C"
