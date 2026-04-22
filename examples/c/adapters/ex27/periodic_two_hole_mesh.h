#ifndef CMFEM_PERIODIC_TWO_HOLE_MESH_H
#define CMFEM_PERIODIC_TWO_HOLE_MESH_H

#include "../../../c_api/common.h"

CMFEM_BEGIN_EXTERN_C

CMFEM_Mesh *CMFEM_NewPeriodicTwoHoleMesh2d(int ref_levels, double radius);

CMFEM_END_EXTERN_C

#endif
