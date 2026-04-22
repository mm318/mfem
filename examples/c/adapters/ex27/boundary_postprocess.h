#ifndef CMFEM_BOUNDARY_POSTPROCESS_H
#define CMFEM_BOUNDARY_POSTPROCESS_H

#include "../../../c_api/common.h"

CMFEM_BEGIN_EXTERN_C

double CMFEM_IntegrateBoundaryConditionGfAi(const CMFEM_GridFunction
                                            *grid_function,
                                            const CMFEM_ArrayInt *marker,
                                            double alpha,
                                            double beta,
                                            double gamma,
                                            double *error);

CMFEM_END_EXTERN_C

#endif
