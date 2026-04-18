#ifndef CMFEM_BLOCK_OPERATOR_H
#define CMFEM_BLOCK_OPERATOR_H

#include "common.h"

CMFEM_STORAGE(BlockOperator, 328);

CMFEM_BEGIN_EXTERN_C

CMFEM_BlockOperator *CMFEM_BlockOperator_NewAiAi(
   const CMFEM_ArrayInt *row_offsets,
   const CMFEM_ArrayInt *col_offsets);
void CMFEM_BlockOperator_Delete(CMFEM_BlockOperator *block_operator);
void CMFEM_BlockOperator_SetBlockSm(CMFEM_BlockOperator *block_operator,
                                    int row_block,
                                    int col_block,
                                    CMFEM_SparseMatrix *op);
void CMFEM_BlockOperator_Mult(const CMFEM_BlockOperator *block_operator,
                              const CMFEM_Vector *x,
                              CMFEM_Vector *y);
void CMFEM_BlockOperator_MultTranspose(
   const CMFEM_BlockOperator *block_operator,
   const CMFEM_Vector *x,
   CMFEM_Vector *y);

CMFEM_END_EXTERN_C

#endif
