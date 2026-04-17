//                              MFEM C Example 2
//
// Description: This example code solves a simple linear elasticity problem
//              describing a multi-material cantilever beam.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("beam-tri.mesh");
   int order = 1;
   int static_cond = 0;
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
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Read the mesh from the given mesh file.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   dim = CMFEM_Mesh_Dimension(mesh);

   if (CMFEM_Mesh_AttributesMax(mesh) < 2 || CMFEM_Mesh_BoundaryAttributesMax(mesh) < 2)
   {
      fprintf(stderr,
              "\nInput mesh should have at least two materials and two boundary attributes!\n\n");
      CMFEM_Mesh_Delete(mesh);
      return 3;
   }

   // 3. Select the order of the finite element discretization space. For NURBS
   //    meshes, increase the order by degree elevation.
   if (CMFEM_Mesh_HasNURBSext(mesh))
   {
      CMFEM_Mesh_DegreeElevate(mesh, order, order);
   }

   // 4. Refine the mesh to increase the resolution.
   ref_levels = cmfem_uniform_refinement_levels(5000.0, CMFEM_Mesh_GetNE(mesh), dim);
   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }

   // 5. Define a vector finite element space on the mesh.
   CMFEM_H1_FECollection *fec = NULL;
   CMFEM_FiniteElementSpace *fespace = NULL;
   if (CMFEM_Mesh_HasNURBSext(mesh))
   {
      fespace = CMFEM_Mesh_GetNodesFESpace(mesh);
   }
   else
   {
      fec = CMFEM_H1_FECollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1VDim(mesh, fec, dim);
      delete_fec = 1;
      delete_fespace = 1;
   }

   printf("Number of finite element unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));
   printf("Assembling: ");
   fflush(stdout);

   // 6. Determine the list of true essential boundary dofs.
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list =
      CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
   CMFEM_ArrayInt_Assign(&ess_tdof_list, 0);
   CMFEM_ArrayInt_Set(&ess_tdof_list, 0, 1);
   {
      _Alignas(max_align_t) CMFEM_ArrayInt true_dofs = CMFEM_ArrayInt_Construct();
      CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_tdof_list, &true_dofs);
      CMFEM_ArrayInt_Destroy(&ess_tdof_list);
      ess_tdof_list = true_dofs;
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side
   //    of the discrete system.
   CMFEM_VectorArrayCoefficient *force = CMFEM_VectorArrayCoefficient_New(dim);
   for (i = 0; i < dim - 1; i++)
   {
      CMFEM_VectorArrayCoefficient_SetConstantCoefficient(
         force, i, CMFEM_ConstantCoefficient_New(0.0));
   }
   {
      CMFEM_Vector *pull_force = CMFEM_Vector_NewSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
      CMFEM_Vector_Assign(pull_force, 0.0);
      CMFEM_Vector_Set(pull_force, 1, -1.0e-2);
      CMFEM_VectorArrayCoefficient_SetPWConstCoefficient(
         force, dim - 1, CMFEM_PWConstCoefficient_New(pull_force));
      CMFEM_Vector_Delete(pull_force);
   }

   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_LinearForm_AddBoundaryIntegrator_VectorBoundaryLFIntegrator(b, force);
   printf("r.h.s. ... ");
   fflush(stdout);
   CMFEM_LinearForm_Assemble(b);

   // 8. Define the solution vector x as a finite element grid function and
   //    initialize it with zero.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_Assign(x, 0.0);

   // 9. Set up the bilinear form a(.,.) corresponding to the linear elasticity
   //    operator with piece-wise constant lambda and mu coefficients.
   CMFEM_Vector *lambda_values = CMFEM_Vector_NewSize(CMFEM_Mesh_AttributesMax(mesh));
   CMFEM_Vector *mu_values = CMFEM_Vector_NewSize(CMFEM_Mesh_AttributesMax(mesh));
   CMFEM_Vector_Assign(lambda_values, 1.0);
   CMFEM_Vector_Assign(mu_values, 1.0);
   CMFEM_Vector_Set(lambda_values, 0, 50.0);
   CMFEM_Vector_Set(mu_values, 0, 50.0);
   CMFEM_PWConstCoefficient *lambda = CMFEM_PWConstCoefficient_New(lambda_values);
   CMFEM_PWConstCoefficient *mu = CMFEM_PWConstCoefficient_New(mu_values);
   CMFEM_Vector_Delete(lambda_values);
   CMFEM_Vector_Delete(mu_values);

   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   CMFEM_BilinearForm_AddDomainIntegrator_Elasticity(a, lambda, mu);
   printf("matrix ... ");
   fflush(stdout);
   if (static_cond)
   {
      CMFEM_BilinearForm_EnableStaticCondensation(a);
   }
   CMFEM_BilinearForm_Assemble(a);

   // 10. Assemble the linear system A X = B.
   _Alignas(max_align_t) CMFEM_SparseMatrix A = CMFEM_SparseMatrix_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
   CMFEM_BilinearForm_FormLinearSystemSparseMatrix(a, &ess_tdof_list, x, b, &A, &X, &B);
   printf("done.\n");
   printf("Size of linear system: %d\n", CMFEM_SparseMatrix_Height(&A));

   // 11. Define a simple symmetric Gauss-Seidel preconditioner and use it to
   //     solve the system Ax=b with PCG.
   CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewSparseMatrix(&A);
   CMFEM_PCG_SparseMatrixGSSmoother(&A, M, &B, &X, 1, 500, 1e-8, 0.0);

   // 12. Recover the solution as a finite element grid function.
   CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);

   // 13. For non-NURBS meshes, make the mesh curved based on the finite
   //     element space.
   if (!CMFEM_Mesh_HasNURBSext(mesh))
   {
      CMFEM_Mesh_SetNodalFESpace(mesh, fespace);
   }

   // 14. Save the displaced mesh and the inverted solution.
   CMFEM_Mesh_AddDisplacementToNodes(mesh, x);
   CMFEM_GridFunction_Scale(x, -1.0);
   CMFEM_Mesh_Print(mesh, "displaced.mesh", 8);
   CMFEM_GridFunction_Save(x, "sol.gf", 8);

   // 15. Send the displaced solution by socket to a GLVis server.
   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
   }

   // 16. Free the used memory.
   CMFEM_GSSmoother_Delete(M);
   CMFEM_SparseMatrix_Destroy(&A);
   CMFEM_Vector_Destroy(&B);
   CMFEM_Vector_Destroy(&X);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_PWConstCoefficient_Delete(lambda);
   CMFEM_PWConstCoefficient_Delete(mu);
   CMFEM_GridFunction_Delete(x);
   CMFEM_LinearForm_Delete(b);
   CMFEM_VectorArrayCoefficient_Delete(force);
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
   return 0;
}
