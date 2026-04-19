//                              MFEM C Example 26
//
// Description: This example demonstrates a partial-assembly geometric
//              multigrid preconditioner built on a hierarchy of H1 spaces.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int geometric_refinements = 0;
   int order_refinements = 2;
   const char *device_config = "cpu";
   int visualization = 1;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-gr",
                                      "--geometric-refinements",
                                      &geometric_refinements);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-or",
                                      "--order-refinements",
                                      &order_refinements);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device",
                                         &device_config);
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

   // 2. Configure the device backend.
   CMFEM_Device *device = CMFEM_Device_New(device_config);
   CMFEM_Device_Print(device);

   // 3. Load the mesh and perform the initial uniform refinement sweep.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   const int dim = CMFEM_Mesh_Dimension(mesh);
   {
      const int ref_levels = (int)floor(log(5000.0 / (double)CMFEM_Mesh_GetNE(mesh)) /
                                        log(2.0) / (double)dim);
      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }
   }

   // 4. Build the H1 hierarchy and the diffusion multigrid preconditioner.
   _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
      CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
   CMFEM_ArrayInt_Assign(&ess_bdr, 1);
   CMFEM_DiffusionMultigrid *multigrid =
      CMFEM_DiffusionMultigrid_NewMeshGrOrAi(mesh,
                                             geometric_refinements,
                                             order_refinements,
                                             &ess_bdr);
   CMFEM_DiffusionMultigrid_SetCycleTypeV(multigrid, 1, 1);

   CMFEM_FiniteElementSpace *fespace =
      CMFEM_DiffusionMultigrid_GetFinestFESpace(multigrid);
   CMFEM_Mesh *fine_mesh = CMFEM_DiffusionMultigrid_GetFinestMesh(multigrid);
   printf("Number of finite element unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   // 5. Assemble the unit right-hand side on the finest space.
   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_LinearForm_AddDomainIntegratorDliCc(b, one);
   CMFEM_LinearForm_Assemble(b);

   // 6. Define the solution grid function and form the fine linear system.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
   CMFEM_GridFunction_Assign(x, 0.0);
   CMFEM_DiffusionMultigrid_FormFineLinearSystemOp(multigrid, x, b, &A, &X, &B);
   printf("Size of linear system: %d\n", CMFEM_OperatorPtr_Height(&A));

   // 7. Solve the fine linear system with PCG preconditioned by multigrid.
   CMFEM_DiffusionMultigrid_PCGOp(multigrid, &A, &B, &X, 1, 2000, 1e-12, 0.0);

   // 8. Recover the solution as a finite element grid function.
   CMFEM_DiffusionMultigrid_RecoverFineFEMSolution(multigrid, &X, b, x);

   // 9. Save the refined mesh and solution.
   CMFEM_Mesh_Print(fine_mesh, "refined.mesh", 8);
   CMFEM_GridFunction_Save(x, "sol.gf", 8);

   // 10. Send the solution to GLVis.
   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(fine_mesh, x, "localhost", 19916);
   }

   // 11. Free the allocated resources.
   CMFEM_Vector_Destroy(&X);
   CMFEM_Vector_Destroy(&B);
   CMFEM_OperatorPtr_Destroy(&A);
   CMFEM_GridFunction_Delete(x);
   CMFEM_ConstantCoefficient_Delete(one);
   CMFEM_LinearForm_Delete(b);
   CMFEM_DiffusionMultigrid_Delete(multigrid);
   CMFEM_ArrayInt_Destroy(&ess_bdr);
   CMFEM_Device_Delete(device);

   return 0;
}
