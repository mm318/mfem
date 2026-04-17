#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_VectorFunctionCoefficient,
                  mfem::VectorFunctionCoefficient);

} // namespace

extern "C" {

   CMFEM_VectorFunctionCoefficient *CMFEM_VectorFunctionCoefficient_New(
      int dim,
      CMFEM_VectorFunctionCallback callback,
      void *context)
   {
      auto function = [callback, context](const mfem::Vector &x, mfem::Vector &value)
      {
         callback(reinterpret_cast<const CMFEM_Vector *>(&x),
                  reinterpret_cast<CMFEM_Vector *>(&value),
                  context);
      };
      return reinterpret_cast<CMFEM_VectorFunctionCoefficient *>(
                new mfem::VectorFunctionCoefficient(dim, std::move(function)));
   }

   void CMFEM_VectorFunctionCoefficient_Delete(CMFEM_VectorFunctionCoefficient
                                               *coefficient)
   {
      delete cmfem::As<mfem::VectorFunctionCoefficient>(coefficient);
   }

} // extern "C"
