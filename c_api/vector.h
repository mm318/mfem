#ifndef CMFEM_VECTOR_H
#define CMFEM_VECTOR_H

#include "common.h"

CMFEM_STORAGE(Vector, 40);

CMFEM_BEGIN_EXTERN_C

CMFEM_Vector CMFEM_Vector_Construct(void);
CMFEM_Vector CMFEM_Vector_ConstructSize(int size);
CMFEM_Vector *CMFEM_Vector_New(void);
CMFEM_Vector *CMFEM_Vector_NewSize(int size);
CMFEM_Vector CMFEM_Vector_Copy(const CMFEM_Vector *vector);
CMFEM_Vector *CMFEM_Vector_NewCopy(const CMFEM_Vector *vector);
void CMFEM_Vector_Destroy(CMFEM_Vector *vector);
void CMFEM_Vector_Delete(CMFEM_Vector *vector);
void CMFEM_Vector_Assign(CMFEM_Vector *vector, double value);
void CMFEM_Vector_Set(CMFEM_Vector *vector, int index, double value);
double CMFEM_Vector_Get(const CMFEM_Vector *vector, int index);
int CMFEM_Vector_Size(const CMFEM_Vector *vector);
void CMFEM_Vector_Neg(CMFEM_Vector *vector);
void CMFEM_Vector_Add(CMFEM_Vector *vector, const CMFEM_Vector *other);
void CMFEM_Vector_SetSubVectorAi(CMFEM_Vector *vector,
                                 const CMFEM_ArrayInt *indices,
                                 double value);

CMFEM_END_EXTERN_C

#endif
