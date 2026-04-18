#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_BlockOperator, mfem::BlockOperator);

} // namespace

extern "C" {

   CMFEM_BlockOperator *CMFEM_BlockOperator_NewAiAi(
      const CMFEM_ArrayInt *row_offsets,
      const CMFEM_ArrayInt *col_offsets)
   {
      return reinterpret_cast<CMFEM_BlockOperator *>(
                new mfem::BlockOperator(cmfem::ArrayIntRef(row_offsets),
                                        cmfem::ArrayIntRef(col_offsets)));
   }

   void CMFEM_BlockOperator_Delete(CMFEM_BlockOperator *block_operator)
   {
      delete cmfem::As<mfem::BlockOperator>(block_operator);
   }

   void CMFEM_BlockOperator_SetBlockSm(CMFEM_BlockOperator *block_operator,
                                       int row_block,
                                       int col_block,
                                       CMFEM_SparseMatrix *op)
   {
      cmfem::As<mfem::BlockOperator>(block_operator)->SetBlock(
         row_block, col_block, cmfem::As<mfem::SparseMatrix>(op));
   }

   void CMFEM_BlockOperator_Mult(const CMFEM_BlockOperator *block_operator,
                                 const CMFEM_Vector *x,
                                 CMFEM_Vector *y)
   {
      cmfem::As<const mfem::BlockOperator>(block_operator)->Mult(
         cmfem::VectorRef(x),
         cmfem::VectorRef(y));
   }

   void CMFEM_BlockOperator_MultTranspose(
      const CMFEM_BlockOperator *block_operator,
      const CMFEM_Vector *x,
      CMFEM_Vector *y)
   {
      cmfem::As<const mfem::BlockOperator>(block_operator)->MultTranspose(
         cmfem::VectorRef(x),
         cmfem::VectorRef(y));
   }

} // extern "C"
