//                              MFEM C Example 19
//
// Description: This example solves a quasi-static incompressible nonlinear
//              elasticity problem using a focused CMFEM bridge for the block
//              nonlinear solve.

#include "common.h"
#include "cmfem.h"
#include "adapters/ex19/incompressible_rubber.h"

#include <stdio.h>

static void reference_configuration(const CMFEM_Vector *x,
                                    CMFEM_Vector *y,
                                    void *context)
{
   int i;
   (void)context;
   for (i = 0; i < CMFEM_Vector_Size(x); i++)
   {
      CMFEM_Vector_Set(y, i, CMFEM_Vector_Get(x, i));
   }
}

static void initial_deformation(const CMFEM_Vector *x,
                                CMFEM_Vector *y,
                                void *context)
{
   int i;
   (void)context;
   for (i = 0; i < CMFEM_Vector_Size(x); i++)
   {
      CMFEM_Vector_Set(y, i, CMFEM_Vector_Get(x, i));
   }
   if (CMFEM_Vector_Size(x) > 1)
   {
      CMFEM_Vector_Set(y,
                       1,
                       CMFEM_Vector_Get(x, 1) + 0.25 * CMFEM_Vector_Get(x, 0));
   }
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("beam-tet.mesh");
   int ref_levels = 0;
   int order = 2;
   int visualization = 1;
   double newton_rel_tol = 1.0e-4;
   double newton_abs_tol = 1.0e-6;
   int newton_iter = 500;
   double mu = 1.0;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-r", "--refine",
                                      &ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-rel",
                                         "--relative-tolerance",
                                         &newton_rel_tol);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-abs",
                                         "--absolute-tolerance",
                                         &newton_abs_tol);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-it",
                                      "--newton-iterations", &newton_iter);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-mu",
                                         "--shear-modulus", &mu);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Read and refine the mesh.
   {
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      const int dim = CMFEM_Mesh_Dimension(mesh);
      CMFEM_H1FeCollection *quad_coll = NULL;
      CMFEM_H1FeCollection *lin_coll = NULL;
      CMFEM_FiniteElementSpace *R_space = NULL;
      CMFEM_FiniteElementSpace *W_space = NULL;
      CMFEM_GridFunction *x_gf = NULL;
      CMFEM_GridFunction *x_ref = NULL;
      CMFEM_GridFunction *x_def = NULL;
      CMFEM_GridFunction *p_gf = NULL;
      CMFEM_VectorFunctionCoefficient *deform = NULL;
      CMFEM_VectorFunctionCoefficient *refconfig = NULL;
      void *oper = NULL;

      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr_u = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr_p = CMFEM_ArrayInt_Construct();

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      // 3. Define the Taylor-Hood displacement and pressure spaces.
      quad_coll = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      lin_coll = CMFEM_H1FeCollection_NewOrderDim(order - 1, dim);
      R_space = CMFEM_FiniteElementSpace_NewMeshH1VDimByVdim(mesh,
                                                             quad_coll,
                                                             dim);
      W_space = CMFEM_FiniteElementSpace_NewMeshH1(mesh, lin_coll);

      printf("***********************************************************\n");
      printf("dim(u) = %d\n", CMFEM_FiniteElementSpace_GetTrueVSize(R_space));
      printf("dim(p) = %d\n", CMFEM_FiniteElementSpace_GetTrueVSize(W_space));
      printf("dim(u+p) = %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(R_space) +
             CMFEM_FiniteElementSpace_GetTrueVSize(W_space));
      printf("***********************************************************\n");

      // 4. Define the Dirichlet conditions on boundary attributes 1 and 2 for
      //    the displacement and on none for the pressure.
      ess_bdr_u = CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(
                                                  mesh));
      ess_bdr_p = CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(
                                                  mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr_u, 0);
      CMFEM_ArrayInt_Assign(&ess_bdr_p, 0);
      if (CMFEM_ArrayInt_Size(&ess_bdr_u) > 0)
      {
         CMFEM_ArrayInt_Set(&ess_bdr_u, 0, 1);
      }
      if (CMFEM_ArrayInt_Size(&ess_bdr_u) > 1)
      {
         CMFEM_ArrayInt_Set(&ess_bdr_u, 1, 1);
      }

      // 5. Define the current, reference, deformation, and pressure fields.
      x_gf = CMFEM_GridFunction_New(R_space);
      x_ref = CMFEM_GridFunction_New(R_space);
      x_def = CMFEM_GridFunction_New(R_space);
      p_gf = CMFEM_GridFunction_New(W_space);

      deform = CMFEM_VectorFunctionCoefficient_New(dim, initial_deformation, NULL);
      refconfig = CMFEM_VectorFunctionCoefficient_New(dim,
                                                      reference_configuration,
                                                      NULL);
      CMFEM_GridFunction_ProjectCoefficientVfc(x_gf, deform);
      CMFEM_GridFunction_ProjectCoefficientVfc(x_ref, refconfig);
      CMFEM_GridFunction_Assign(p_gf, 0.0);
      CMFEM_GridFunction_SetTrueVector(x_gf);
      CMFEM_GridFunction_SetTrueVector(p_gf);

      // 6. Initialize the incompressible neo-Hookean operator and solve the
      //    nonlinear system.
      oper = CMFEM_Ex19RubberOperator_New(R_space,
                                          W_space,
                                          &ess_bdr_u,
                                          &ess_bdr_p,
                                          newton_rel_tol,
                                          newton_abs_tol,
                                          newton_iter,
                                          mu);
      CMFEM_Ex19RubberOperator_Solve(oper, x_gf, p_gf);

      // 7. Compute the final deformation and save the output fields.
      CMFEM_GridFunction_CopyFromGf(x_def, x_gf);
      CMFEM_GridFunction_AddScaledGf(x_def, -1.0, x_ref);

      CMFEM_Mesh_SetNodalFESpace(mesh, R_space);
      CMFEM_Mesh_AddDisplacementToNodes(mesh, x_def);
      CMFEM_Mesh_Print(mesh, "deformed.mesh", 8);
      CMFEM_GridFunction_Save(p_gf, "pressure.sol", 8);
      CMFEM_GridFunction_Save(x_def, "deformation.sol", 8);

      // 8. Send the final fields to GLVis if requested.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, x_def, "localhost", 19916);
         CMFEM_SendSolutionToGLVis(mesh, p_gf, "localhost", 19916);
      }

      // 9. Free the used memory.
      CMFEM_Ex19RubberOperator_Delete(oper);
      CMFEM_VectorFunctionCoefficient_Delete(refconfig);
      CMFEM_VectorFunctionCoefficient_Delete(deform);
      CMFEM_GridFunction_Delete(p_gf);
      CMFEM_GridFunction_Delete(x_def);
      CMFEM_GridFunction_Delete(x_ref);
      CMFEM_GridFunction_Delete(x_gf);
      CMFEM_FiniteElementSpace_Delete(W_space);
      CMFEM_FiniteElementSpace_Delete(R_space);
      CMFEM_H1FeCollection_Delete(lin_coll);
      CMFEM_H1FeCollection_Delete(quad_coll);
      CMFEM_ArrayInt_Destroy(&ess_bdr_p);
      CMFEM_ArrayInt_Destroy(&ess_bdr_u);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
