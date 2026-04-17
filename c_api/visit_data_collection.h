#ifndef CMFEM_VISIT_DATA_COLLECTION_H
#define CMFEM_VISIT_DATA_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(VisItDataCollection, 216);

CMFEM_BEGIN_EXTERN_C

CMFEM_VisItDataCollection *CMFEM_VisItDataCollection_New(
   const char *collection_name,
   CMFEM_Mesh *mesh);
void CMFEM_VisItDataCollection_Delete(CMFEM_VisItDataCollection *collection);
void CMFEM_VisItDataCollection_RegisterFieldGf(
   CMFEM_VisItDataCollection *collection,
   const char *field_name,
   CMFEM_GridFunction *grid_function);
void CMFEM_VisItDataCollection_SetCycle(CMFEM_VisItDataCollection *collection,
                                        int cycle);
void CMFEM_VisItDataCollection_SetTime(CMFEM_VisItDataCollection *collection,
                                       double time);
void CMFEM_VisItDataCollection_Save(CMFEM_VisItDataCollection *collection);

CMFEM_END_EXTERN_C

#endif
