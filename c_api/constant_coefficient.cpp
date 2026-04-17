#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ConstantCoefficient, mfem::ConstantCoefficient);

} // namespace

extern "C" {

   CMFEM_ConstantCoefficient *CMFEM_ConstantCoefficient_New(double value)
   {
      return reinterpret_cast<CMFEM_ConstantCoefficient *>(
                new mfem::ConstantCoefficient(static_cast<mfem::real_t>(value)));
   }

   void CMFEM_ConstantCoefficient_Delete(CMFEM_ConstantCoefficient *coefficient)
   {
      delete cmfem::As<mfem::ConstantCoefficient>(coefficient);
   }

} // extern "C"
