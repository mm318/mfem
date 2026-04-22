#ifndef CMFEM_DEVICE_H
#define CMFEM_DEVICE_H

#include "../common.h"

CMFEM_STORAGE(Device, 40);

CMFEM_BEGIN_EXTERN_C

CMFEM_Device *CMFEM_Device_New(const char *device_config);
void CMFEM_Device_Delete(CMFEM_Device *device);
void CMFEM_Device_Print(CMFEM_Device *device);
int CMFEM_Device_IsEnabled(void);

CMFEM_END_EXTERN_C

#endif
