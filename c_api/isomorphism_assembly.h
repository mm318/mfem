#ifndef CMFEM_ISOMORPHISM_ASSEMBLY_H
#define CMFEM_ISOMORPHISM_ASSEMBLY_H

#include "common.h"

CMFEM_BEGIN_EXTERN_C

void CMFEM_AssembleIsomorphismLinearFormRtGf(
   const CMFEM_FiniteElementSpace *fespace,
   const CMFEM_GridFunction *psi,
   CMFEM_Vector *out);
CMFEM_SparseMatrix *CMFEM_AssembleDIsomorphismMassMatrixRtGf(
   const CMFEM_FiniteElementSpace *fespace,
   const CMFEM_GridFunction *psi,
   double eps);

CMFEM_END_EXTERN_C

#endif
