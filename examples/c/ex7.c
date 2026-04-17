//                                MFEM C Example 7
//
// Description: This example code demonstrates the use of MFEM to define a
//              triangulation of a unit sphere and a simple isoparametric
//              finite element discretization of the screened Poisson problem,
//              -Delta u + u = f.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

static double analytic_solution(const CMFEM_Vector *x, void *context)
{
   (void)context;
   const double x0 = CMFEM_Vector_Get(x, 0);
   const double x1 = CMFEM_Vector_Get(x, 1);
   const double x2 = CMFEM_Vector_Get(x, 2);
   const double l2 = x0 * x0 + x1 * x1 + x2 * x2;
   return x0 * x1 / l2;
}

static double analytic_rhs(const CMFEM_Vector *x, void *context)
{
   (void)context;
   const double x0 = CMFEM_Vector_Get(x, 0);
   const double x1 = CMFEM_Vector_Get(x, 1);
   const double x2 = CMFEM_Vector_Get(x, 2);
   const double l2 = x0 * x0 + x1 * x1 + x2 * x2;
   return 7.0 * x0 * x1 / l2;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int elem_type = 1;
   int ref_levels = 2;
   int amr = 0;
   int order = 2;
   int always_snap = 0;
   int visualization = 1;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_int_option(argc, argv, &i, "-e", "--elem",
                                          &elem_type);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-r", "--refine",
                                      &ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-amr",
                                      "--refine-locally", &amr);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-snap", "--always-snap",
                                       "-no-snap", "--snap-at-the-end",
                                       &always_snap);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Generate an initial high-order (surface) mesh on the unit sphere.
   {
      static const double tri_v[6][3] =
      {
         { 1.0,  0.0,  0.0}, { 0.0,  1.0,  0.0}, {-1.0,  0.0,  0.0},
         { 0.0, -1.0,  0.0}, { 0.0,  0.0,  1.0}, { 0.0,  0.0, -1.0}
      };
      static const int tri_e[8][3] =
      {
         {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4},
         {1, 0, 5}, {2, 1, 5}, {3, 2, 5}, {0, 3, 5}
      };
      static const double quad_v[8][3] =
      {
         {-1.0, -1.0, -1.0}, {+1.0, -1.0, -1.0}, {+1.0, +1.0, -1.0},
         {-1.0, +1.0, -1.0}, {-1.0, -1.0, +1.0}, {+1.0, -1.0, +1.0},
         {+1.0, +1.0, +1.0}, {-1.0, +1.0, +1.0}
      };
      static const int quad_e[6][4] =
      {
         {3, 2, 1, 0}, {0, 1, 5, 4}, {1, 2, 6, 5},
         {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}
      };

      const int num_vertices = (elem_type == 0) ? 6 : 8;
      const int num_elements = (elem_type == 0) ? 8 : 6;
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewDimensionVerticesElementsBoundarySpace(
                            2, num_vertices, num_elements, 0, 3);
      if (elem_type == 0)
      {
         for (i = 0; i < num_vertices; i++)
         {
            CMFEM_Mesh_AddVertex(mesh, tri_v[i]);
         }
         for (i = 0; i < num_elements; i++)
         {
            CMFEM_Mesh_AddTriangle(mesh, tri_e[i], i + 1);
         }
         CMFEM_Mesh_FinalizeTriMesh(mesh, 1, 1, 1);
      }
      else
      {
         for (i = 0; i < num_vertices; i++)
         {
            CMFEM_Mesh_AddVertex(mesh, quad_v[i]);
         }
         for (i = 0; i < num_elements; i++)
         {
            CMFEM_Mesh_AddQuad(mesh, quad_e[i], i + 1);
         }
         CMFEM_Mesh_FinalizeQuadMesh(mesh, 1, 1, 1);
      }

      // 3. Set the space for the high-order mesh nodes and refine while
      //    snapping nodes back to the unit sphere.
      {
         const int dim = CMFEM_Mesh_Dimension(mesh);
         const int space_dim = CMFEM_Mesh_SpaceDimension(mesh);
         CMFEM_H1_FECollection *fec =
            CMFEM_H1_FECollection_NewOrderDim(order, dim);
         CMFEM_FiniteElementSpace *nodal_fes =
            CMFEM_FiniteElementSpace_NewMeshH1VDim(mesh, fec, space_dim);
         CMFEM_Mesh_SetNodalFESpace(mesh, nodal_fes);

         for (i = 0; i <= ref_levels; i++)
         {
            if (i > 0)
            {
               CMFEM_Mesh_UniformRefinement(mesh);
            }
            if (always_snap || i == ref_levels)
            {
               CMFEM_Mesh_SnapNodesToUnitSphere(mesh);
            }
         }

         if (amr == 1)
         {
            int level;
            for (level = 0; level < 5; level++)
            {
               CMFEM_Mesh_RefineAtVertex3(mesh, 0.0, 0.0, 1.0);
            }
            CMFEM_Mesh_SnapNodesToUnitSphere(mesh);
         }
         else if (amr == 2)
         {
            int level;
            for (level = 0; level < 4; level++)
            {
               CMFEM_Mesh_RandomRefinement(mesh, 0.5);
            }
            CMFEM_Mesh_SnapNodesToUnitSphere(mesh);
         }

         // 4. Define a finite element space on the mesh using the same nodal
         //    finite elements.
         CMFEM_FiniteElementSpace *fespace =
            CMFEM_FiniteElementSpace_NewMeshH1(mesh, fec);
         printf("Number of unknowns: %d\n",
                CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

         // 5. Set up the linear form b(.) for the right-hand side.
         CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
         CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
         CMFEM_FunctionCoefficient *rhs_coef =
            CMFEM_FunctionCoefficient_New(analytic_rhs, NULL);
         CMFEM_FunctionCoefficient *sol_coef =
            CMFEM_FunctionCoefficient_New(analytic_solution, NULL);
         CMFEM_LinearForm_AddDomainIntegrator_DomainLFIntegrator_FunctionCoefficient(
            b, rhs_coef);
         CMFEM_LinearForm_Assemble(b);

         // 6. Define the solution vector x and initialize it with zero.
         CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
         CMFEM_GridFunction_Assign(x, 0.0);

         // 7. Set up the bilinear form a(.,.) with diffusion and mass terms.
         CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
         CMFEM_BilinearForm_AddDomainIntegrator_DiffusionCoefficient(a, one);
         CMFEM_BilinearForm_AddDomainIntegrator_MassCoefficient(a, one);

         // 8. Assemble the linear system, apply conforming constraints, etc.
         CMFEM_BilinearForm_Assemble(a);
         _Alignas(max_align_t) CMFEM_SparseMatrix A = CMFEM_SparseMatrix_Construct();
         _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
         _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
         _Alignas(max_align_t) CMFEM_ArrayInt empty_tdof_list =
            CMFEM_ArrayInt_Construct();
         CMFEM_BilinearForm_FormLinearSystemSparseMatrix(a, &empty_tdof_list, x,
                                                         b, &A, &X, &B);

         // 9. Solve the linear system with symmetric Gauss-Seidel
         //    preconditioning and PCG.
         {
            CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewSparseMatrix(&A);
            CMFEM_PCG_SparseMatrixGSSmoother(&A, M, &B, &X, 1, 200, 1e-12, 0.0);
            CMFEM_GSSmoother_Delete(M);
         }

         // 10. Recover the solution as a finite element grid function.
         CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);

         // 11. Compute and print the L2 norm of the error.
         printf("\nL2 norm of error: %.14g\n",
                CMFEM_GridFunction_ComputeL2ErrorFunctionCoefficient(x,
                                                                     sol_coef));

         // 12. Save the refined mesh and the solution.
         CMFEM_Mesh_Print(mesh, "sphere_refined.mesh", 8);
         CMFEM_GridFunction_Save(x, "sol.gf", 8);

         // 13. Send the solution by socket to a GLVis server.
         if (visualization)
         {
            CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
         }

         CMFEM_ArrayInt_Destroy(&empty_tdof_list);
         CMFEM_Vector_Destroy(&X);
         CMFEM_Vector_Destroy(&B);
         CMFEM_SparseMatrix_Destroy(&A);
         CMFEM_BilinearForm_Delete(a);
         CMFEM_GridFunction_Delete(x);
         CMFEM_FunctionCoefficient_Delete(sol_coef);
         CMFEM_FunctionCoefficient_Delete(rhs_coef);
         CMFEM_ConstantCoefficient_Delete(one);
         CMFEM_LinearForm_Delete(b);
         CMFEM_FiniteElementSpace_Delete(fespace);
         CMFEM_FiniteElementSpace_Delete(nodal_fes);
         CMFEM_H1_FECollection_Delete(fec);
      }

      // 14. Free the used memory.
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
