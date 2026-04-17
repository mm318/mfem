//                              MFEM C Example 0
//
// Description: This example code demonstrates the most basic usage of the
//              CMFEM wrapper to define a simple finite element discretization
//              of the Poisson problem -Delta u = 1 with zero Dirichlet
//              boundary conditions.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh", &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Read the mesh from the given mesh file, and refine once uniformly.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   CMFEM_Mesh_UniformRefinement(mesh);

   // 3. Define a finite element space on the mesh. Here we use H1 continuous
   //    high-order Lagrange finite elements of the given order.
   CMFEM_H1_FECollection *fec =
      CMFEM_H1_FECollection_NewOrderDim(order, CMFEM_Mesh_Dimension(mesh));
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, fec);
   printf("Number of unknowns: %d\n", CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   // 4. Extract the list of all the boundary DOFs. These will be marked as
   //    Dirichlet in order to enforce zero boundary conditions.
   _Alignas(max_align_t) CMFEM_ArrayInt boundary_dofs = CMFEM_ArrayInt_Construct();
   CMFEM_FiniteElementSpace_GetBoundaryTrueDofs(fespace, &boundary_dofs);

   // 5. Define the solution x as a finite element grid function in fespace. Set
   //    the initial guess to zero, which also sets the boundary conditions.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_Assign(x, 0.0);

   // 6. Set up the linear form b(.) corresponding to the right-hand side.
   CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_LinearForm_AddDomainIntegrator_DomainLFIntegrator_ConstantCoefficient(b, one);
   CMFEM_LinearForm_Assemble(b);

   // 7. Set up the bilinear form a(.,.) corresponding to the -Delta operator.
   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   CMFEM_BilinearForm_AddDomainIntegrator_Diffusion(a);
   CMFEM_BilinearForm_Assemble(a);

   // 8. Form the linear system A X = B. This includes eliminating boundary
   //    conditions, applying AMR constraints, and other transformations.
   _Alignas(max_align_t) CMFEM_SparseMatrix A = CMFEM_SparseMatrix_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
   CMFEM_BilinearForm_FormLinearSystemSparseMatrix(a, &boundary_dofs, x, b, &A, &X, &B);

   // 9. Solve the system using PCG with symmetric Gauss-Seidel preconditioner.
   CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewSparseMatrix(&A);
   CMFEM_PCG_SparseMatrixGSSmoother(&A, M, &B, &X, 1, 200, 1e-12, 0.0);

   // 10. Recover the solution x as a grid function and save to file.
   CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
   CMFEM_GridFunction_Save(x, "sol.gf", 16);
   CMFEM_Mesh_Save(mesh, "mesh.mesh", 16);

   CMFEM_GSSmoother_Delete(M);
   CMFEM_SparseMatrix_Destroy(&A);
   CMFEM_Vector_Destroy(&B);
   CMFEM_Vector_Destroy(&X);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_LinearForm_Delete(b);
   CMFEM_ConstantCoefficient_Delete(one);
   CMFEM_GridFunction_Delete(x);
   CMFEM_ArrayInt_Destroy(&boundary_dofs);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_H1_FECollection_Delete(fec);
   CMFEM_Mesh_Delete(mesh);
   return 0;
}
