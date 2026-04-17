#ifndef CMFEM_VISUALIZATION_H
#define CMFEM_VISUALIZATION_H

#include "common.h"

CMFEM_BEGIN_EXTERN_C

void CMFEM_SendSolutionToGLVis(const CMFEM_Mesh *mesh,
                               const CMFEM_GridFunction *grid_function,
                               const char *host,
                               int port);

CMFEM_END_EXTERN_C

#endif
