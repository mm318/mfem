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

   void CMFEM_FunctionCoefficient_Delete(CMFEM_FunctionCoefficient *coefficient)
   {
      delete cmfem::As<mfem::FunctionCoefficient>(coefficient);
   }

} // extern "C"
