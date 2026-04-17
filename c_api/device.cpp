#include "common.hpp"

namespace
{

CMFEM_ASSERT_TYPE(CMFEM_Device, mfem::Device);

} // namespace

extern "C" {

   CMFEM_Device *CMFEM_Device_New(const char *device_config)
   {
      return reinterpret_cast<CMFEM_Device *>(new mfem::Device(device_config));
   }

   void CMFEM_Device_Delete(CMFEM_Device *device)
   {
      delete cmfem::As<mfem::Device>(device);
   }

   void CMFEM_Device_Print(CMFEM_Device *device)
   {
      cmfem::As<mfem::Device>(device)->Print();
   }

   int CMFEM_Device_IsEnabled(void)
   {
      return mfem::Device::IsEnabled();
   }

} // extern "C"
