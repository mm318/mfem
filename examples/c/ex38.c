//                                MFEM Example 38
//
// Compile with: make ex38
//
// Sample runs:
// (since all sample runs require LAPACK or ALGOIM, the * symbol is used to
//  exclude them from the automatically generated internal MFEM tests).
//              * ex38
//              * ex38 -i volumetric1d
//              * ex38 -i surface2d
//              * ex38 -i surface2d -o 4 -r 5 -m 1
//              * ex38 -i volumetric2d
//              * ex38 -i volumetric2d -o 4 -r 5 -m 1
//              * ex38 -i surface3d
//              * ex38 -i surface3d -o 3 -r 4 -m 1
//              * ex38 -i volumetric3d
//              * ex38 -i volumetric3d -o 3 -r 4 -m 1
//
// Description: This example code demonstrates the use of MFEM to integrate
//              functions over implicit interfaces and subdomains bounded by
//              implicit interfaces.
//
//              The quadrature rules are constructed by means of moment-fitting.
//              The interface is given by the zero isoline of a level-set
//              function phi and the subdomain is given as the domain where
//              phi>0 holds.

#include "cmfem.h"

#include <stdio.h>

enum
{
   cmfem_mfem_skip_return_value_ = 242
};

int main(int argc, char *argv[])
{
   (void)argc;
   (void)argv;

#if defined(MFEM_USE_LAPACK) || defined(MFEM_USE_ALGOIM)
   fprintf(stderr,
           "ex38.c is not implemented yet for LAPACK/ALGOIM-enabled builds.\n");
   return 1;
#else
   puts("MFEM must be built with LAPACK or ALGOIM for this example.");
   return cmfem_mfem_skip_return_value_;
#endif
}
