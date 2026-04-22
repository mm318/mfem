#ifndef CMFEM_SPARSE_MATRIX_H
#define CMFEM_SPARSE_MATRIX_H

#include "../common.h"

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
void CMFEM_SparseMatrix_MultTranspose(const CMFEM_SparseMatrix *matrix,
                                      const CMFEM_Vector *x,
                                      CMFEM_Vector *y);
double CMFEM_SparseMatrix_InnerProduct(const CMFEM_SparseMatrix *matrix,
                                       const CMFEM_Vector *x,
                                       const CMFEM_Vector *y);
void CMFEM_SparseMatrix_GetDiag(const CMFEM_SparseMatrix *matrix,
                                CMFEM_Vector *diag);
void CMFEM_SparseMatrix_Scale(CMFEM_SparseMatrix *matrix, double scale);
void CMFEM_SparseMatrix_ScaleRow(CMFEM_SparseMatrix *matrix,
                                 int row,
                                 double scale);
CMFEM_SparseMatrix *CMFEM_TransposeSm(const CMFEM_SparseMatrix *matrix);
CMFEM_SparseMatrix *CMFEM_MultSmSm(const CMFEM_SparseMatrix *A,
                                   const CMFEM_SparseMatrix *B);
CMFEM_SparseMatrix *CMFEM_MultAtDASmVec(const CMFEM_SparseMatrix *A,
                                        const CMFEM_Vector *diag);
CMFEM_SparseMatrix *CMFEM_RAPSmSmSm(const CMFEM_SparseMatrix *Rt,
                                    const CMFEM_SparseMatrix *A,
                                    const CMFEM_SparseMatrix *P);
CMFEM_SparseMatrix *CMFEM_AddSmSm(double a,
                                  const CMFEM_SparseMatrix *A,
                                  double b,
                                  const CMFEM_SparseMatrix *B);

CMFEM_END_EXTERN_C

#endif
