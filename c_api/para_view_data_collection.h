#ifndef CMFEM_PARA_VIEW_DATA_COLLECTION_H
#define CMFEM_PARA_VIEW_DATA_COLLECTION_H

#include "common.h"

CMFEM_STORAGE(ParaViewDataCollection, 584);

CMFEM_BEGIN_EXTERN_C

CMFEM_ParaViewDataCollection *CMFEM_ParaViewDataCollection_New(
   const char *collection_name,
   CMFEM_Mesh *mesh);
void CMFEM_ParaViewDataCollection_Delete(
   CMFEM_ParaViewDataCollection *collection);
void CMFEM_ParaViewDataCollection_SetPrefixPath(
   CMFEM_ParaViewDataCollection *collection,
   const char *prefix_path);
void CMFEM_ParaViewDataCollection_SetLevelsOfDetail(
   CMFEM_ParaViewDataCollection *collection,
   int levels_of_detail);
void CMFEM_ParaViewDataCollection_SetCycle(
   CMFEM_ParaViewDataCollection *collection,
   int cycle);
void CMFEM_ParaViewDataCollection_SetDataFormatBinary(
   CMFEM_ParaViewDataCollection *collection);
void CMFEM_ParaViewDataCollection_SetHighOrderOutput(
   CMFEM_ParaViewDataCollection *collection,
   int high_order_output);
void CMFEM_ParaViewDataCollection_SetTime(
   CMFEM_ParaViewDataCollection *collection,
   double time);
void CMFEM_ParaViewDataCollection_RegisterFieldGf(
   CMFEM_ParaViewDataCollection *collection,
   const char *field_name,
   CMFEM_GridFunction *grid_function);
void CMFEM_ParaViewDataCollection_Save(
   CMFEM_ParaViewDataCollection *collection);

CMFEM_END_EXTERN_C

#endif
