#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_VectorArrayCoefficient, mfem::VectorArrayCoefficient);

} // namespace

extern "C" {

CMFEM_VectorArrayCoefficient *CMFEM_VectorArrayCoefficient_New(int dim)
{
   return reinterpret_cast<CMFEM_VectorArrayCoefficient *>(
      new mfem::VectorArrayCoefficient(dim));
}

void CMFEM_VectorArrayCoefficient_Delete(CMFEM_VectorArrayCoefficient *coefficient)
{
   delete cmfem::As<mfem::VectorArrayCoefficient>(coefficient);
}

void CMFEM_VectorArrayCoefficient_SetConstantCoefficient(
   CMFEM_VectorArrayCoefficient *coefficient,
   int component,
   CMFEM_ConstantCoefficient *entry)
{
   cmfem::As<mfem::VectorArrayCoefficient>(coefficient)->Set(
      component,
      cmfem::As<mfem::ConstantCoefficient>(entry));
}

void CMFEM_VectorArrayCoefficient_SetPWConstCoefficient(
   CMFEM_VectorArrayCoefficient *coefficient,
   int component,
   CMFEM_PWConstCoefficient *entry)
{
   cmfem::As<mfem::VectorArrayCoefficient>(coefficient)->Set(
      component,
      cmfem::As<mfem::PWConstCoefficient>(entry));
}

} // extern "C"
