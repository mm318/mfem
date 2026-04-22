#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_GridFunctionCoefficient, mfem::GridFunctionCoefficient);

} // namespace

extern "C" {

   CMFEM_GridFunctionCoefficient *CMFEM_GridFunctionCoefficient_NewGf(
      const CMFEM_GridFunction *grid_function)
   {
      return reinterpret_cast<CMFEM_GridFunctionCoefficient *>(
                new mfem::GridFunctionCoefficient(cmfem::As<const mfem::GridFunction>(
                                                     grid_function)));
   }

   void CMFEM_GridFunctionCoefficient_Delete(
      CMFEM_GridFunctionCoefficient *coefficient)
   {
      delete cmfem::As<mfem::GridFunctionCoefficient>(coefficient);
   }

} // extern "C"
