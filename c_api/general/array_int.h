#ifndef CMFEM_ARRAY_INT_H
#define CMFEM_ARRAY_INT_H

#include "../common.h"

CMFEM_STORAGE(ArrayInt, 32);

CMFEM_BEGIN_EXTERN_C

CMFEM_ArrayInt CMFEM_ArrayInt_Construct(void);
CMFEM_ArrayInt CMFEM_ArrayInt_ConstructSize(int size);
CMFEM_ArrayInt *CMFEM_ArrayInt_New(void);
CMFEM_ArrayInt *CMFEM_ArrayInt_NewSize(int size);
CMFEM_ArrayInt CMFEM_ArrayInt_Copy(const CMFEM_ArrayInt *array);
CMFEM_ArrayInt *CMFEM_ArrayInt_NewCopy(const CMFEM_ArrayInt *array);
void CMFEM_ArrayInt_Destroy(CMFEM_ArrayInt *array);
void CMFEM_ArrayInt_Delete(CMFEM_ArrayInt *array);
void CMFEM_ArrayInt_Assign(CMFEM_ArrayInt *array, int value);
void CMFEM_ArrayInt_Set(CMFEM_ArrayInt *array, int index, int value);
int CMFEM_ArrayInt_Get(const CMFEM_ArrayInt *array, int index);
int CMFEM_ArrayInt_Size(const CMFEM_ArrayInt *array);
int CMFEM_ArrayInt_Max(const CMFEM_ArrayInt *array);

CMFEM_END_EXTERN_C

#endif
