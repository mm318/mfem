#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_RAPOperator, mfem::RAPOperator);

} // namespace

extern "C" {

   CMFEM_RAPOperator *CMFEM_RAPOperator_NewBopSmBop(
      const CMFEM_BlockOperator *Rt,
      const CMFEM_SparseMatrix *A,
      const CMFEM_BlockOperator *P)
   {
      return reinterpret_cast<CMFEM_RAPOperator *>(
                new mfem::RAPOperator(
                   *cmfem::As<const mfem::BlockOperator>(Rt),
                   *cmfem::As<const mfem::SparseMatrix>(A),
                   *cmfem::As<const mfem::BlockOperator>(P)));
   }

   void CMFEM_RAPOperator_Delete(CMFEM_RAPOperator *rap_operator)
   {
      delete cmfem::As<mfem::RAPOperator>(rap_operator);
   }

} // extern "C"
