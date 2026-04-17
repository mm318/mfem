#ifndef CMFEM_SOLVERS_H
#define CMFEM_SOLVERS_H

#include "common.h"

CMFEM_BEGIN_EXTERN_C

void CMFEM_PCGSmGs(const CMFEM_SparseMatrix *A,
                   CMFEM_GSSmoother *M,
                   const CMFEM_Vector *B,
                   CMFEM_Vector *X,
                   int print_iter,
                   int max_iter,
                   double rtol,
                   double atol);
void CMFEM_PCGOpGs(const CMFEM_OperatorPtr *A,
                   CMFEM_GSSmoother *M,
                   const CMFEM_Vector *B,
                   CMFEM_Vector *X,
                   int print_iter,
                   int max_iter,
                   double rtol,
                   double atol);
void CMFEM_PCGOpOjs(const CMFEM_OperatorPtr *A,
                    CMFEM_OperatorJacobiSmoother *M,
                    const CMFEM_Vector *B,
                    CMFEM_Vector *X,
                    int print_iter,
                    int max_iter,
                    double rtol,
                    double atol);
void CMFEM_CGOp(const CMFEM_OperatorPtr *A,
                const CMFEM_Vector *B,
                CMFEM_Vector *X,
                int print_iter,
                int max_iter,
                double rtol,
                double atol);
void CMFEM_GMRESOpGs(const CMFEM_OperatorPtr *A,
                     CMFEM_GSSmoother *M,
                     const CMFEM_Vector *B,
                     CMFEM_Vector *X,
                     int print_iter,
                     int max_iter,
                     int restart,
                     double rtol,
                     double atol);
int CMFEM_UsesTensorBasis(const CMFEM_FiniteElementSpace *fespace);

CMFEM_END_EXTERN_C

#endif
