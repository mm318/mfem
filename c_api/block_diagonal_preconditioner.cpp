#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_BlockDiagonalPreconditioner,
                  mfem::BlockDiagonalPreconditioner);

} // namespace

extern "C" {

   CMFEM_BlockDiagonalPreconditioner *CMFEM_BlockDiagonalPreconditioner_NewAi(
      const CMFEM_ArrayInt *offsets)
   {
      return reinterpret_cast<CMFEM_BlockDiagonalPreconditioner *>(
                new mfem::BlockDiagonalPreconditioner(cmfem::ArrayIntRef(offsets)));
   }

   void CMFEM_BlockDiagonalPreconditioner_Delete(
      CMFEM_BlockDiagonalPreconditioner *preconditioner)
   {
      delete cmfem::As<mfem::BlockDiagonalPreconditioner>(preconditioner);
   }

   void CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockCgs(
      CMFEM_BlockDiagonalPreconditioner *preconditioner,
      int block,
      CMFEM_CGSolver *solver)
   {
      cmfem::As<mfem::BlockDiagonalPreconditioner>(preconditioner)
      ->SetDiagonalBlock(block, cmfem::As<mfem::CGSolver>(solver));
   }

} // extern "C"
