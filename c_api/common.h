#ifndef CMFEM_COMMON_H
#define CMFEM_COMMON_H

#include <stddef.h>

#if defined(__cplusplus)
#define CMFEM_BEGIN_EXTERN_C extern "C" {
#define CMFEM_END_EXTERN_C }
#else
#define CMFEM_BEGIN_EXTERN_C
#define CMFEM_END_EXTERN_C
#endif

typedef struct CMFEM_ArrayInt CMFEM_ArrayInt;
typedef struct CMFEM_Mesh CMFEM_Mesh;
typedef struct CMFEM_H1_FECollection CMFEM_H1_FECollection;
typedef struct CMFEM_ND_FECollection CMFEM_ND_FECollection;
typedef struct CMFEM_RT_FECollection CMFEM_RT_FECollection;
typedef struct CMFEM_DG_FECollection CMFEM_DG_FECollection;
typedef struct CMFEM_DG_Interface_FECollection CMFEM_DG_Interface_FECollection;
typedef struct CMFEM_FiniteElementSpace CMFEM_FiniteElementSpace;
typedef struct CMFEM_GridFunction CMFEM_GridFunction;
typedef struct CMFEM_LinearForm CMFEM_LinearForm;
typedef struct CMFEM_BilinearForm CMFEM_BilinearForm;
typedef struct CMFEM_SparseMatrix CMFEM_SparseMatrix;
typedef struct CMFEM_Vector CMFEM_Vector;
typedef struct CMFEM_GSSmoother CMFEM_GSSmoother;
typedef struct CMFEM_Device CMFEM_Device;
typedef struct CMFEM_OperatorPtr CMFEM_OperatorPtr;
typedef struct CMFEM_OperatorJacobiSmoother CMFEM_OperatorJacobiSmoother;
typedef struct CMFEM_ConstantCoefficient CMFEM_ConstantCoefficient;
typedef struct CMFEM_FunctionCoefficient CMFEM_FunctionCoefficient;
typedef struct CMFEM_PWConstCoefficient CMFEM_PWConstCoefficient;
typedef struct CMFEM_VectorArrayCoefficient CMFEM_VectorArrayCoefficient;
typedef struct CMFEM_VectorFunctionCoefficient CMFEM_VectorFunctionCoefficient;

#define CMFEM_STORAGE(Name, Size) struct CMFEM_##Name { char data[Size]; }

#endif
