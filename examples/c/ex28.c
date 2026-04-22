//                              MFEM C Example 28
//
// Description: This example demonstrates sliding normal constraints in a
//              small elasticity problem using MFEM's constrained solver path.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

static CMFEM_Mesh *build_trapezoid_mesh(double offset)
{
   CMFEM_Mesh *mesh;
   double vertex[2];
   int indices[4];
   int segment[2];

   if (offset >= 0.9)
   {
      fprintf(stderr, "offset is too large\n");
      return NULL;
   }

   mesh = CMFEM_Mesh_NewDimensionVerticesElementsBoundarySpace(2, 4, 1, 4, 2);

   vertex[0] = 0.0;
   vertex[1] = 0.0;
   CMFEM_Mesh_AddVertex(mesh, vertex);
   vertex[0] = 1.0;
   vertex[1] = 0.0;
   CMFEM_Mesh_AddVertex(mesh, vertex);
   vertex[0] = offset;
   vertex[1] = 1.0;
   CMFEM_Mesh_AddVertex(mesh, vertex);
   vertex[0] = 1.0;
   vertex[1] = 1.0;
   CMFEM_Mesh_AddVertex(mesh, vertex);

   indices[0] = 0;
   indices[1] = 1;
   indices[2] = 3;
   indices[3] = 2;
   CMFEM_Mesh_AddQuad(mesh, indices, 1);

   segment[0] = 0;
   segment[1] = 1;
   CMFEM_Mesh_AddBdrSegment(mesh, segment, 1);
   segment[0] = 1;
   segment[1] = 3;
   CMFEM_Mesh_AddBdrSegment(mesh, segment, 2);
   segment[0] = 2;
   segment[1] = 3;
   CMFEM_Mesh_AddBdrSegment(mesh, segment, 3);
   segment[0] = 0;
   segment[1] = 2;
   CMFEM_Mesh_AddBdrSegment(mesh, segment, 4);

   CMFEM_Mesh_FinalizeQuadMesh(mesh, 1, 0, 1);
   return mesh;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int order = 1;
   int visualization = 1;
   double offset = 0.3;
   int visit = 0;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "--offset", "--offset",
                                         &offset);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-visit", "--visit-datafiles",
                                       "-no-visit", "--no-visit-datafiles",
                                       &visit);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Build the trapezoidal mesh and uniformly refine it.
   {
      CMFEM_Mesh *mesh = build_trapezoid_mesh(offset);
      CMFEM_H1FeCollection *fec;
      CMFEM_FiniteElementSpace *fespace;
      CMFEM_GridFunction *x;
      CMFEM_VectorArrayCoefficient *force;
      CMFEM_LinearForm *b;
      CMFEM_Vector *lambda_values;
      CMFEM_Vector *mu_values;
      CMFEM_PWConstCoefficient *lambda_func;
      CMFEM_PWConstCoefficient *mu_func;
      CMFEM_BilinearForm *a;
      CMFEM_ArrayInt ess_bdr;
      CMFEM_ArrayInt ess_tdof_list;
      CMFEM_ArrayInt constraint_atts;
      CMFEM_ArrayInt lagrange_rowstarts;
      CMFEM_SparseMatrix *local_constraints;
      CMFEM_GSSmoother *M;
      CMFEM_SchurConstrainedSolver *solver;
      CMFEM_VisItDataCollection *visit_dc = NULL;
      int dim;

      _Alignas(max_align_t) CMFEM_SparseMatrix A = CMFEM_SparseMatrix_Construct();
      _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
      int ref_levels;

      if (mesh == NULL)
      {
         return 1;
      }
      dim = CMFEM_Mesh_Dimension(mesh);
      ref_levels = cmfem_uniform_refinement_levels(1000.0,
                                                   CMFEM_Mesh_GetNE(mesh),
                                                   dim);
      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      // 3. Define the vector H1 finite element space for displacement.
      fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1VDim(mesh, fec, dim);
      printf("Number of finite element unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(fespace));
      printf("Assembling matrix and r.h.s... ");
      fflush(stdout);

      // 4. Build the empty essential boundary list and the normal-direction
      //    constraints used by the Schur constrained solver.
      ess_bdr = CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr, 0);
      ess_tdof_list = CMFEM_ArrayInt_Construct();
      CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr,
                                                    &ess_tdof_list);

      constraint_atts = CMFEM_ArrayInt_ConstructSize(2);
      CMFEM_ArrayInt_Set(&constraint_atts, 0, 1);
      CMFEM_ArrayInt_Set(&constraint_atts, 1, 4);
      lagrange_rowstarts = CMFEM_ArrayInt_Construct();
      local_constraints = CMFEM_BuildNormalConstraints(
                             fespace, &constraint_atts, &lagrange_rowstarts);

      // 5. Set up the boundary push force on attribute 2.
      force = CMFEM_VectorArrayCoefficient_New(dim);
      for (i = 0; i < dim - 1; i++)
      {
         CMFEM_VectorArrayCoefficient_SetConstantCoefficient(
            force, i, CMFEM_ConstantCoefficient_New(0.0));
      }
      {
         CMFEM_Vector *push_force = CMFEM_Vector_NewSize(
                                       CMFEM_Mesh_BoundaryAttributesMax(mesh));
         CMFEM_Vector_Assign(push_force, 0.0);
         CMFEM_Vector_Set(push_force, 1, -5.0e-2);
         CMFEM_VectorArrayCoefficient_SetPWConstCoefficient(
            force, 0, CMFEM_PWConstCoefficient_New(push_force));
         CMFEM_Vector_Delete(push_force);
      }
      b = CMFEM_LinearForm_New(fespace);
      CMFEM_LinearForm_AddBoundaryIntegratorVbl(b, force);

      // 6. Assemble the elasticity bilinear form with piecewise-constant Lamé
      //    coefficients.
      lambda_values = CMFEM_Vector_NewSize(CMFEM_Mesh_AttributesMax(mesh));
      mu_values = CMFEM_Vector_NewSize(CMFEM_Mesh_AttributesMax(mesh));
      CMFEM_Vector_Assign(lambda_values, 1.0);
      CMFEM_Vector_Assign(mu_values, 1.0);
      lambda_func = CMFEM_PWConstCoefficient_New(lambda_values);
      mu_func = CMFEM_PWConstCoefficient_New(mu_values);
      CMFEM_Vector_Delete(lambda_values);
      CMFEM_Vector_Delete(mu_values);

      a = CMFEM_BilinearForm_New(fespace);
      CMFEM_BilinearForm_AddDomainIntegratorEi(a, lambda_func, mu_func);
      CMFEM_BilinearForm_Assemble(a);
      CMFEM_LinearForm_Assemble(b);

      // 7. Form the unconstrained linear system and then solve the additional
      //    normal constraints with the Schur-constrained solver.
      x = CMFEM_GridFunction_New(fespace);
      CMFEM_GridFunction_Assign(x, 0.0);
      CMFEM_BilinearForm_FormLinearSystemSm(a, &ess_tdof_list, x, b, &A, &X, &B);
      printf("done.\n");
      printf("Size of linear system: %d\n", CMFEM_SparseMatrix_Height(&A));

      M = CMFEM_GSSmoother_NewSm(&A);
      solver = CMFEM_SchurConstrainedSolver_NewSmSmGs(&A, local_constraints, M);
      CMFEM_SchurConstrainedSolver_SetRelTol(solver, 1e-5);
      CMFEM_SchurConstrainedSolver_SetMaxIter(solver, 2000);
      CMFEM_SchurConstrainedSolver_SetPrintLevel(solver, 1);
      CMFEM_SchurConstrainedSolver_Mult(solver, &B, &X);

      // 8. Recover the displacement field and deform the mesh for output.
      CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
      CMFEM_Mesh_SetNodalFESpace(mesh, fespace);
      CMFEM_Mesh_AddDisplacementToNodes(mesh, x);

      // 9. Save VisIt output for the displaced configuration if requested.
      if (visit)
      {
         visit_dc = CMFEM_VisItDataCollection_New("ex28", mesh);
         CMFEM_VisItDataCollection_SetLevelsOfDetail(visit_dc, 4);
         CMFEM_VisItDataCollection_RegisterFieldGf(visit_dc, "displacement", x);
         CMFEM_VisItDataCollection_Save(visit_dc);
      }

      // 10. Save the displaced mesh and the backward displacement field.
      CMFEM_GridFunction_Scale(x, -1.0);
      CMFEM_Mesh_Print(mesh, "displaced.mesh", 8);
      CMFEM_GridFunction_Save(x, "sol.gf", 8);

      // 11. Send the displaced solution to GLVis.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
      }

      if (visit_dc) { CMFEM_VisItDataCollection_Delete(visit_dc); }
      CMFEM_SchurConstrainedSolver_Delete(solver);
      CMFEM_GSSmoother_Delete(M);
      CMFEM_SparseMatrix_Delete(local_constraints);
      CMFEM_SparseMatrix_Destroy(&A);
      CMFEM_Vector_Destroy(&B);
      CMFEM_Vector_Destroy(&X);
      CMFEM_ArrayInt_Destroy(&lagrange_rowstarts);
      CMFEM_ArrayInt_Destroy(&constraint_atts);
      CMFEM_ArrayInt_Destroy(&ess_tdof_list);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_BilinearForm_Delete(a);
      CMFEM_PWConstCoefficient_Delete(mu_func);
      CMFEM_PWConstCoefficient_Delete(lambda_func);
      CMFEM_LinearForm_Delete(b);
      CMFEM_VectorArrayCoefficient_Delete(force);
      CMFEM_GridFunction_Delete(x);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_H1FeCollection_Delete(fec);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
