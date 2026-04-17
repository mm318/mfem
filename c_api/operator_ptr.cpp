#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_OperatorPtr, mfem::OperatorPtr);

} // namespace

extern "C" {

   CMFEM_OperatorPtr CMFEM_OperatorPtr_Construct(void)
   {
      alignas(mfem::OperatorPtr) CMFEM_OperatorPtr op;
      new (cmfem::As<mfem::OperatorPtr>(&op)) mfem::OperatorPtr();
      return op;
   }

   CMFEM_OperatorPtr *CMFEM_OperatorPtr_New(void)
   {
      return reinterpret_cast<CMFEM_OperatorPtr *>(new mfem::OperatorPtr());
   }

   CMFEM_OperatorPtr CMFEM_OperatorPtr_Copy(const CMFEM_OperatorPtr *op)
   {
      alignas(mfem::OperatorPtr) CMFEM_OperatorPtr copy;
      new (cmfem::As<mfem::OperatorPtr>(&copy)) mfem::OperatorPtr(
         cmfem::OperatorPtrRef(op));
      return copy;
   }

   CMFEM_OperatorPtr *CMFEM_OperatorPtr_NewCopy(const CMFEM_OperatorPtr *op)
   {
      return reinterpret_cast<CMFEM_OperatorPtr *>(
                new mfem::OperatorPtr(cmfem::OperatorPtrRef(op)));
   }

   void CMFEM_OperatorPtr_Destroy(CMFEM_OperatorPtr *op)
   {
      cmfem::As<mfem::OperatorPtr>(op)->~OperatorHandle();
   }

   void CMFEM_OperatorPtr_Delete(CMFEM_OperatorPtr *op)
   {
      delete cmfem::As<mfem::OperatorPtr>(op);
   }

   int CMFEM_OperatorPtr_Height(const CMFEM_OperatorPtr *op)
   {
      return cmfem::OperatorPtrRef(op)->Height();
   }

} // extern "C"
