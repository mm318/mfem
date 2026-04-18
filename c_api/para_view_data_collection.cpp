#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ParaViewDataCollection, mfem::ParaViewDataCollection);

} // namespace

extern "C" {

   CMFEM_ParaViewDataCollection *CMFEM_ParaViewDataCollection_New(
      const char *collection_name,
      CMFEM_Mesh *mesh)
   {
      return reinterpret_cast<CMFEM_ParaViewDataCollection *>(
                new mfem::ParaViewDataCollection(collection_name,
                                                 cmfem::As<mfem::Mesh>(mesh)));
   }

   void CMFEM_ParaViewDataCollection_Delete(
      CMFEM_ParaViewDataCollection *collection)
   {
      delete cmfem::As<mfem::ParaViewDataCollection>(collection);
   }

   void CMFEM_ParaViewDataCollection_SetPrefixPath(
      CMFEM_ParaViewDataCollection *collection,
      const char *prefix_path)
   {
      cmfem::As<mfem::ParaViewDataCollection>(collection)->SetPrefixPath(
         prefix_path);
   }

   void CMFEM_ParaViewDataCollection_SetLevelsOfDetail(
      CMFEM_ParaViewDataCollection *collection,
      int levels_of_detail)
   {
      cmfem::As<mfem::ParaViewDataCollection>(collection)->SetLevelsOfDetail(
         levels_of_detail);
   }

   void CMFEM_ParaViewDataCollection_SetCycle(
      CMFEM_ParaViewDataCollection *collection,
      int cycle)
   {
      cmfem::As<mfem::ParaViewDataCollection>(collection)->SetCycle(cycle);
   }

   void CMFEM_ParaViewDataCollection_SetDataFormatBinary(
      CMFEM_ParaViewDataCollection *collection)
   {
      cmfem::As<mfem::ParaViewDataCollection>(collection)->SetDataFormat(
         mfem::VTKFormat::BINARY);
   }

   void CMFEM_ParaViewDataCollection_SetHighOrderOutput(
      CMFEM_ParaViewDataCollection *collection,
      int high_order_output)
   {
      cmfem::As<mfem::ParaViewDataCollection>(collection)->SetHighOrderOutput(
         high_order_output != 0);
   }

   void CMFEM_ParaViewDataCollection_SetTime(
      CMFEM_ParaViewDataCollection *collection,
      double time)
   {
      cmfem::As<mfem::ParaViewDataCollection>(collection)->SetTime(
         static_cast<mfem::real_t>(time));
   }

   void CMFEM_ParaViewDataCollection_RegisterFieldGf(
      CMFEM_ParaViewDataCollection *collection,
      const char *field_name,
      CMFEM_GridFunction *grid_function)
   {
      cmfem::As<mfem::ParaViewDataCollection>(collection)->RegisterField(
         field_name,
         cmfem::As<mfem::GridFunction>(grid_function));
   }

   void CMFEM_ParaViewDataCollection_Save(
      CMFEM_ParaViewDataCollection *collection)
   {
      cmfem::As<mfem::ParaViewDataCollection>(collection)->Save();
   }

} // extern "C"
