#ifndef CMFEM_BLOCK_DIAGONAL_PRECONDITIONER_H
#define CMFEM_BLOCK_DIAGONAL_PRECONDITIONER_H

#include "common.h"

CMFEM_STORAGE(BlockDiagonalPreconditioner, 208);

CMFEM_BEGIN_EXTERN_C

CMFEM_BlockDiagonalPreconditioner *CMFEM_BlockDiagonalPreconditioner_NewAi(
   const CMFEM_ArrayInt *offsets);
void CMFEM_BlockDiagonalPreconditioner_Delete(
   CMFEM_BlockDiagonalPreconditioner *preconditioner);
void CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockCgs(
   CMFEM_BlockDiagonalPreconditioner *preconditioner,
   int block,
   CMFEM_CGSolver *solver);

CMFEM_END_EXTERN_C

#endif
