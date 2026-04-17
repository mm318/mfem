//                              MFEM C Example 1
//
// Description: This example code demonstrates the use of the CMFEM wrapper to
//              define a simple finite element discretization of the Poisson
//              problem -Delta u = 1 with homogeneous Dirichlet boundary
//              conditions.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int static_cond = 0;
   int pa = 0;
   int fa = 0;
   const char *device_config = "cpu";
   int visualization = 1;
   int dim;
   int ref_levels;
   int delete_fec = 0;
   int delete_fespace = 0;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh", &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-sc", "--static-condensation",
                                       "-no-sc", "--no-static-condensation",
                                       &static_cond);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-pa", "--partial-assembly",
                                       "-no-pa", "--no-partial-assembly",
                                       &pa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-fa", "--full-assembly",
                                       "-no-fa", "--no-full-assembly",
                                       &fa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device", &device_config);
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

   // 2. Enable hardware devices and programming models based on command line
   //    options.
   CMFEM_Device *device = CMFEM_Device_New(device_config);
   CMFEM_Device_Print(device);

   // 3. Read the mesh from the given mesh file.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   dim = CMFEM_Mesh_Dimension(mesh);
   // 4. Refine the mesh to increase the resolution.
   ref_levels = cmfem_uniform_refinement_levels(50000.0, CMFEM_Mesh_GetNE(mesh), dim);
   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }

   // 5. Define a finite element space on the mesh. If order < 1, use an
   //    isoparametric space when the mesh already has nodes.
   CMFEM_H1_FECollection *fec = NULL;
   CMFEM_FiniteElementSpace *fespace = NULL;
   if (order > 0)
   {
      fec = CMFEM_H1_FECollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, fec);
      delete_fec = 1;
      delete_fespace = 1;
   }
   else if (CMFEM_Mesh_HasNodes(mesh))
   {
      fespace = CMFEM_Mesh_GetNodesFESpace(mesh);
      printf("Using isoparametric FEs: %s\n", CMFEM_Mesh_GetNodesOwnFECName(mesh));
   }
   else
   {
      order = 1;
      fec = CMFEM_H1_FECollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, fec);
      delete_fec = 1;
      delete_fespace = 1;
   }

   printf("Number of finite element unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   // 6. Determine the list of true essential boundary dofs.
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
   if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
   {
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
         CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr, 0);
      CMFEM_Mesh_MarkExternalBoundaries(mesh, &ess_bdr);
      CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr, &ess_tdof_list);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side
   //    of the discrete system.
   CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_LinearForm_AddDomainIntegrator_DomainLFIntegrator_ConstantCoefficient(b, one);
   CMFEM_LinearForm_Assemble(b);

   // 8. Define the solution vector x as a finite element grid function and
   //    initialize it with zero, which satisfies the boundary conditions.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_Assign(x, 0.0);

   // 9. Set up the bilinear form a(.,.) corresponding to the Laplacian.
   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   if (pa)
   {
      CMFEM_BilinearForm_SetAssemblyLevelPartial(a);
   }
   if (fa)
   {
      CMFEM_BilinearForm_SetAssemblyLevelFull(a);
      CMFEM_BilinearForm_EnableSparseMatrixSorting(a, CMFEM_Device_IsEnabled());
   }
   CMFEM_BilinearForm_AddDomainIntegrator_DiffusionCoefficient(a, one);
   if (static_cond)
   {
      CMFEM_BilinearForm_EnableStaticCondensation(a);
   }
   CMFEM_BilinearForm_Assemble(a);

   // 10. Assemble the linear system A X = B.
   _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
   CMFEM_BilinearForm_FormLinearSystemOperator(a, &ess_tdof_list, x, b, &A, &X, &B);
   printf("Size of linear system: %d\n", CMFEM_OperatorPtr_Height(&A));

   // 11. Solve the linear system.
   if (!pa)
   {
      CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewOperator(&A);
      CMFEM_PCG_OperatorGSSmoother(&A, M, &B, &X, 1, 200, 1e-12, 0.0);
      CMFEM_GSSmoother_Delete(M);
   }
   else if (CMFEM_UsesTensorBasis(fespace))
   {
      CMFEM_OperatorJacobiSmoother *M =
         CMFEM_OperatorJacobiSmoother_NewBilinearForm(a, &ess_tdof_list);
      CMFEM_PCG_OperatorJacobiSmoother(&A, M, &B, &X, 1, 400, 1e-12, 0.0);
      CMFEM_OperatorJacobiSmoother_Delete(M);
   }
   else
   {
      CMFEM_CG_Operator(&A, &B, &X, 1, 400, 1e-12, 0.0);
   }

   // 12. Recover the solution as a finite element grid function.
   CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
   // 13. Save the refined mesh and the solution.
   CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
   CMFEM_GridFunction_Save(x, "sol.gf", 8);

   // 14. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
   }

   // 15. Free the used memory.
   CMFEM_OperatorPtr_Destroy(&A);
   CMFEM_Vector_Destroy(&B);
   CMFEM_Vector_Destroy(&X);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_GridFunction_Delete(x);
   CMFEM_LinearForm_Delete(b);
   CMFEM_ConstantCoefficient_Delete(one);
   CMFEM_ArrayInt_Destroy(&ess_tdof_list);
   if (delete_fespace)
   {
      CMFEM_FiniteElementSpace_Delete(fespace);
   }
   if (delete_fec)
   {
      CMFEM_H1_FECollection_Delete(fec);
   }
   CMFEM_Mesh_Delete(mesh);
   CMFEM_Device_Delete(device);
   return 0;
}
