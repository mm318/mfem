#ifndef CMFEM_CONSTRAINTS_H
#define CMFEM_CONSTRAINTS_H

#include "common.h"

CMFEM_BEGIN_EXTERN_C

CMFEM_SparseMatrix *CMFEM_BuildNormalConstraints(
   CMFEM_FiniteElementSpace *fespace,
   CMFEM_ArrayInt *constrained_att,
   CMFEM_ArrayInt *constraint_rowstarts);

CMFEM_END_EXTERN_C

#endif
