#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_VectorConstantCoefficient,
                  mfem::VectorConstantCoefficient);

} // namespace

extern "C" {

   CMFEM_VectorConstantCoefficient *CMFEM_VectorConstantCoefficient_New(
      const CMFEM_Vector *value)
   {
      return reinterpret_cast<CMFEM_VectorConstantCoefficient *>(
                new mfem::VectorConstantCoefficient(cmfem::VectorRef(value)));
   }

   void CMFEM_VectorConstantCoefficient_Delete(
      CMFEM_VectorConstantCoefficient *coefficient)
   {
      delete cmfem::As<mfem::VectorConstantCoefficient>(coefficient);
   }

} // extern "C"
