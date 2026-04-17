#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_VisItDataCollection, mfem::VisItDataCollection);

} // namespace

extern "C" {

   CMFEM_VisItDataCollection *CMFEM_VisItDataCollection_New(
      const char *collection_name,
      CMFEM_Mesh *mesh)
   {
      return reinterpret_cast<CMFEM_VisItDataCollection *>(
                new mfem::VisItDataCollection(collection_name,
                                              cmfem::As<mfem::Mesh>(mesh)));
   }

   void CMFEM_VisItDataCollection_Delete(CMFEM_VisItDataCollection *collection)
   {
      delete cmfem::As<mfem::VisItDataCollection>(collection);
   }

   void CMFEM_VisItDataCollection_RegisterFieldGf(
      CMFEM_VisItDataCollection *collection,
      const char *field_name,
      CMFEM_GridFunction *grid_function)
   {
      cmfem::As<mfem::VisItDataCollection>(collection)->RegisterField(
         field_name,
         cmfem::As<mfem::GridFunction>(grid_function));
   }

   void CMFEM_VisItDataCollection_SetCycle(CMFEM_VisItDataCollection *collection,
                                           int cycle)
   {
      cmfem::As<mfem::VisItDataCollection>(collection)->SetCycle(cycle);
   }

   void CMFEM_VisItDataCollection_SetTime(CMFEM_VisItDataCollection *collection,
                                          double time)
   {
      cmfem::As<mfem::VisItDataCollection>(collection)->SetTime(
         static_cast<mfem::real_t>(time));
   }

   void CMFEM_VisItDataCollection_Save(CMFEM_VisItDataCollection *collection)
   {
      cmfem::As<mfem::VisItDataCollection>(collection)->Save();
   }

} // extern "C"
