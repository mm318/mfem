#ifndef CMFEM_SPARSE_MATRIX_H
#define CMFEM_SPARSE_MATRIX_H

#include "common.h"

CMFEM_STORAGE(SparseMatrix, 144);

CMFEM_BEGIN_EXTERN_C

CMFEM_SparseMatrix CMFEM_SparseMatrix_Construct(void);
CMFEM_SparseMatrix *CMFEM_SparseMatrix_New(void);
CMFEM_SparseMatrix CMFEM_SparseMatrix_Copy(const CMFEM_SparseMatrix *matrix);
CMFEM_SparseMatrix *CMFEM_SparseMatrix_NewCopy(const CMFEM_SparseMatrix
                                               *matrix);
void CMFEM_SparseMatrix_Destroy(CMFEM_SparseMatrix *matrix);
void CMFEM_SparseMatrix_Delete(CMFEM_SparseMatrix *matrix);
int CMFEM_SparseMatrix_Height(const CMFEM_SparseMatrix *matrix);
void CMFEM_SparseMatrix_Mult(const CMFEM_SparseMatrix *matrix,
                             const CMFEM_Vector *x,
                             CMFEM_Vector *y);
CMFEM_SparseMatrix *CMFEM_AddSmSm(double a,
                                  const CMFEM_SparseMatrix *A,
                                  double b,
                                  const CMFEM_SparseMatrix *B);

CMFEM_END_EXTERN_C

#endif
