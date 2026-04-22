#ifndef CMFEM_OPERATOR_PTR_H
#define CMFEM_OPERATOR_PTR_H

#include "../common.h"

CMFEM_STORAGE(OperatorPtr, 16);

CMFEM_BEGIN_EXTERN_C

CMFEM_OperatorPtr CMFEM_OperatorPtr_Construct(void);
CMFEM_OperatorPtr *CMFEM_OperatorPtr_New(void);
CMFEM_OperatorPtr CMFEM_OperatorPtr_Copy(const CMFEM_OperatorPtr *op);
CMFEM_OperatorPtr *CMFEM_OperatorPtr_NewCopy(const CMFEM_OperatorPtr *op);
void CMFEM_OperatorPtr_Destroy(CMFEM_OperatorPtr *op);
void CMFEM_OperatorPtr_Delete(CMFEM_OperatorPtr *op);
int CMFEM_OperatorPtr_Height(const CMFEM_OperatorPtr *op);

CMFEM_END_EXTERN_C

#endif
