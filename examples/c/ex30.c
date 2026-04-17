//                               MFEM C Example 30
//
// Description: This example demonstrates adaptive mesh refinement
//              preprocessing which lowers the data oscillation to a
//              user-defined relative threshold. There is no PDE being solved.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

// Piecewise-affine function which is sometimes mesh-conforming.
static double affine_function(const CMFEM_Vector *p, void *context)
{
   const double x = CMFEM_Vector_Get(p, 0);
   const double y = CMFEM_Vector_Get(p, 1);
   (void)context;
   if (x < 0.0)
   {
      return 1.0 + x + y;
   }
   return 1.0;
}

// Piecewise-constant function which is never mesh-conforming.
static double jump_function(const CMFEM_Vector *p, void *context)
{
   const double radius = hypot(CMFEM_Vector_Get(p, 0), CMFEM_Vector_Get(p, 1));
   (void)context;
   if (radius > 0.4 && radius < 0.6)
   {
      return 1.0;
   }
   return 5.0;
}

// Singular function derived from the "steep wavefront" problem.
static double singular_function(const CMFEM_Vector *p, void *context)
{
   const double x = CMFEM_Vector_Get(p, 0);
   const double y = CMFEM_Vector_Get(p, 1);
   const double alpha = 1000.0;
   const double xc = 0.75;
   const double yc = 0.5;
   const double r0 = 0.7;
   const double r = hypot(x - xc, y - yc);
   const double num = -(alpha - pow(alpha, 3.0) * (r * r - r0 * r0));
   double denom = pow(r * (alpha * alpha * r0 * r0 +
                           alpha * alpha * r * r -
                           2.0 * alpha * alpha * r0 * r + 1.0), 2.0);
   (void)context;
   if (denom < 1.0e-8)
   {
      denom = 1.0e-8;
   }
   return num / denom;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int nc_limit = 1;
   int max_elems = 100 * 1000;
   double max_elems_double = (double)max_elems;
   int visualization = 1;
   double osc_threshold = 1.0e-3;
   int enriched_order = 5;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-l", "--nc-limit",
                                      &nc_limit);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-me", "--max-elems",
                                         &max_elems_double);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-e", "--error",
                                         &osc_threshold);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-eo",
                                      "--enriched_order", &enriched_order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }
   max_elems = (int)max_elems_double;

   // 2. Read the mesh and convert NURBS meshes to piecewise-polynomial curved
   //    meshes when needed.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   if (CMFEM_Mesh_HasNURBSext(mesh))
   {
      for (i = 0; i < 2; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }
      CMFEM_Mesh_SetCurvature(mesh, 2);
   }

   // 3. Define functions and the coefficient refiner.
   CMFEM_FunctionCoefficient *affine_coeff =
      CMFEM_FunctionCoefficient_New(affine_function, NULL);
   CMFEM_FunctionCoefficient *jump_coeff =
      CMFEM_FunctionCoefficient_New(jump_function, NULL);
   CMFEM_FunctionCoefficient *singular_coeff =
      CMFEM_FunctionCoefficient_New(singular_function, NULL);
   CMFEM_CoefficientRefiner *coeff_refiner =
      CMFEM_CoefficientRefiner_NewFc(affine_coeff, order);

   // 4. Apply custom refiner settings.
   CMFEM_CoefficientRefiner_SetIntRuleOrder(coeff_refiner,
                                            2 * order + enriched_order);
   CMFEM_CoefficientRefiner_SetMaxElements(coeff_refiner, max_elems);
   CMFEM_CoefficientRefiner_SetThreshold(coeff_refiner, osc_threshold);
   CMFEM_CoefficientRefiner_SetNCLimit(coeff_refiner, nc_limit);
   CMFEM_CoefficientRefiner_PrintWarnings(coeff_refiner);

   // 5. Preprocess mesh to control oscillation for the piecewise-affine
   //    function.
   CMFEM_CoefficientRefiner_PreprocessMesh(coeff_refiner, mesh);
   printf("\nFunction 0 (affine)\n");
   printf("Number of Elements %d\n", CMFEM_Mesh_GetNE(mesh));
   printf("Osc error %.14g\n", CMFEM_CoefficientRefiner_GetOsc(coeff_refiner));

   // 6. Preprocess mesh to control oscillation for the jump function.
   CMFEM_CoefficientRefiner_ResetCoefficientFc(coeff_refiner, jump_coeff);
   CMFEM_CoefficientRefiner_PreprocessMesh(coeff_refiner, mesh);
   printf("\nFunction 1 (discontinuous)\n");
   printf("Number of Elements %d\n", CMFEM_Mesh_GetNE(mesh));
   printf("Osc error %.14g\n", CMFEM_CoefficientRefiner_GetOsc(coeff_refiner));

   // 7. Preprocess mesh to control oscillation for the singular function.
   CMFEM_CoefficientRefiner_ResetCoefficientFc(coeff_refiner,
                                               singular_coeff);
   CMFEM_CoefficientRefiner_PreprocessMesh(coeff_refiner, mesh);
   printf("\nFunction 2 (singular)\n");
   printf("Number of Elements %d\n", CMFEM_Mesh_GetNE(mesh));
   printf("Osc error %.14g\n", CMFEM_CoefficientRefiner_GetOsc(coeff_refiner));

   // 8. Send the final mesh to GLVis.
   if (visualization)
   {
      CMFEM_SendMeshToGLVis(mesh, "localhost", 19916);
   }

   CMFEM_CoefficientRefiner_Delete(coeff_refiner);
   CMFEM_FunctionCoefficient_Delete(singular_coeff);
   CMFEM_FunctionCoefficient_Delete(jump_coeff);
   CMFEM_FunctionCoefficient_Delete(affine_coeff);
   CMFEM_Mesh_Delete(mesh);
   return 0;
}
