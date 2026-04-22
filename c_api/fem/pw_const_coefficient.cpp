#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_PWConstCoefficient, mfem::PWConstCoefficient);

} // namespace

extern "C" {

   CMFEM_PWConstCoefficient *CMFEM_PWConstCoefficient_New(const CMFEM_Vector
                                                          *values)
   {
      return reinterpret_cast<CMFEM_PWConstCoefficient *>(
                new mfem::PWConstCoefficient(cmfem::VectorRef(values)));
   }

   void CMFEM_PWConstCoefficient_Delete(CMFEM_PWConstCoefficient *coefficient)
   {
      delete cmfem::As<mfem::PWConstCoefficient>(coefficient);
   }

} // extern "C"
