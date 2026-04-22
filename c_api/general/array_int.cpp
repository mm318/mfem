#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_ArrayInt, mfem::Array<int>);

} // namespace

extern "C" {

   CMFEM_ArrayInt CMFEM_ArrayInt_Construct(void)
   {
      alignas(mfem::Array<int>) CMFEM_ArrayInt array;
      new (cmfem::As<mfem::Array<int> >(&array)) mfem::Array<int>();
      return array;
   }

   CMFEM_ArrayInt CMFEM_ArrayInt_ConstructSize(int size)
   {
      alignas(mfem::Array<int>) CMFEM_ArrayInt array;
      new (cmfem::As<mfem::Array<int> >(&array)) mfem::Array<int>(size);
      return array;
   }

   CMFEM_ArrayInt *CMFEM_ArrayInt_New(void)
   {
      return reinterpret_cast<CMFEM_ArrayInt *>(new mfem::Array<int>());
   }

   CMFEM_ArrayInt *CMFEM_ArrayInt_NewSize(int size)
   {
      return reinterpret_cast<CMFEM_ArrayInt *>(new mfem::Array<int>(size));
   }

   CMFEM_ArrayInt CMFEM_ArrayInt_Copy(const CMFEM_ArrayInt *array)
   {
      alignas(mfem::Array<int>) CMFEM_ArrayInt copy;
      new (cmfem::As<mfem::Array<int> >(&copy)) mfem::Array<int>(cmfem::ArrayIntRef(
                                                                    array));
      return copy;
   }

   CMFEM_ArrayInt *CMFEM_ArrayInt_NewCopy(const CMFEM_ArrayInt *array)
   {
      return reinterpret_cast<CMFEM_ArrayInt *>(
                new mfem::Array<int>(cmfem::ArrayIntRef(array)));
   }

   void CMFEM_ArrayInt_Destroy(CMFEM_ArrayInt *array)
   {
      cmfem::As<mfem::Array<int> >(array)->~Array();
   }

   void CMFEM_ArrayInt_Delete(CMFEM_ArrayInt *array)
   {
      delete cmfem::As<mfem::Array<int> >(array);
   }

   void CMFEM_ArrayInt_Assign(CMFEM_ArrayInt *array, int value)
   {
      cmfem::ArrayIntRef(array) = value;
   }

   void CMFEM_ArrayInt_Set(CMFEM_ArrayInt *array, int index, int value)
   {
      cmfem::ArrayIntRef(array)[index] = value;
   }

   int CMFEM_ArrayInt_Get(const CMFEM_ArrayInt *array, int index)
   {
      return cmfem::ArrayIntRef(array)[index];
   }

   int CMFEM_ArrayInt_Size(const CMFEM_ArrayInt *array)
   {
      return cmfem::ArrayIntRef(array).Size();
   }

   int CMFEM_ArrayInt_Max(const CMFEM_ArrayInt *array)
   {
      return cmfem::ArrayIntRef(array).Max();
   }

} // extern "C"
