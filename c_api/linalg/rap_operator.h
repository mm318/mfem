#ifndef CMFEM_RAP_OPERATOR_H
#define CMFEM_RAP_OPERATOR_H

#include "../common.h"

CMFEM_STORAGE(RAPOperator, 128);

CMFEM_BEGIN_EXTERN_C

CMFEM_RAPOperator *CMFEM_RAPOperator_NewBopSmBop(
   const CMFEM_BlockOperator *Rt,
   const CMFEM_SparseMatrix *A,
   const CMFEM_BlockOperator *P);
void CMFEM_RAPOperator_Delete(CMFEM_RAPOperator *rap_operator);

CMFEM_END_EXTERN_C

#endif
