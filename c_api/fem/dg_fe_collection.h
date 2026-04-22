#ifndef CMFEM_DG_FE_COLLECTION_H
#define CMFEM_DG_FE_COLLECTION_H

#include "../common.h"

CMFEM_STORAGE(DgFeCollection, 488);

enum CMFEM_BasisType
{
   CMFEM_BASIS_GAUSS_LEGENDRE = 0,
   CMFEM_BASIS_GAUSS_LOBATTO = 1,
};

CMFEM_BEGIN_EXTERN_C

CMFEM_DgFeCollection *CMFEM_DgFeCollection_NewOrderDimBasis(int order,
                                                            int dim,
                                                            int basis_type);
void CMFEM_DgFeCollection_Delete(CMFEM_DgFeCollection *fec);

CMFEM_END_EXTERN_C

#endif
