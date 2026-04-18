//                              MFEM C Example 29
//
// Description: This example solves a Poisson problem on a curved 2D surface
//              embedded in 3D, using an anisotropic 3x3 diffusion tensor.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

static double u_exact(const CMFEM_Vector *x, void *context)
{
   (void)context;
   return (0.25 * (2.0 + CMFEM_Vector_Get(x, 0)) - CMFEM_Vector_Get(x, 2)) *
          (CMFEM_Vector_Get(x, 2) +
           0.25 * (2.0 + CMFEM_Vector_Get(x, 0)));
}

static void du_exact(const CMFEM_Vector *x, CMFEM_Vector *du, void *context)
{
   (void)context;
   CMFEM_Vector_Set(du, 0,
                    0.125 * (2.0 + CMFEM_Vector_Get(x, 0)) *
                    CMFEM_Vector_Get(x, 1) * CMFEM_Vector_Get(x, 1));
   CMFEM_Vector_Set(du, 1,
                    -0.125 * (2.0 + CMFEM_Vector_Get(x, 0)) *
                    CMFEM_Vector_Get(x, 0) * CMFEM_Vector_Get(x, 1));
   CMFEM_Vector_Set(du, 2, -2.0 * CMFEM_Vector_Get(x, 2));
}

static void sigma_func(const CMFEM_Vector *x,
                       double *values,
                       int dim,
                       void *context)
{
   const double x0 = CMFEM_Vector_Get(x, 0);
   const double x1 = CMFEM_Vector_Get(x, 1);
   const double a = 17.0 - 2.0 * x0 * (1.0 + x0);
   int i;
   (void)context;

   for (i = 0; i < dim * dim; i++)
   {
      values[i] = 0.0;
   }

   values[0 * dim + 0] = 0.5 + x0 * x0 * (8.0 / a - 0.5);
   values[0 * dim + 1] = x0 * x1 * (8.0 / a - 0.5);
   values[1 * dim + 0] = values[0 * dim + 1];
   values[1 * dim + 1] = 0.5 * x0 * x0 + 8.0 * x1 * x1 / a;
   values[2 * dim + 2] = a / 32.0;
}

static void flux_exact(const CMFEM_Vector *x, CMFEM_Vector *f, void *context)
{
   double sigma[9];
   double du[3];
   int row;
   int col;
   (void)context;

   sigma_func(x, sigma, 3, NULL);
   du_exact(x, f, NULL);
   for (row = 0; row < 3; row++)
   {
      du[row] = CMFEM_Vector_Get(f, row);
   }

   for (row = 0; row < 3; row++)
   {
      double value = 0.0;
      for (col = 0; col < 3; col++)
      {
         value += sigma[row * 3 + col] * du[col];
      }
      CMFEM_Vector_Set(f, row, -value);
   }
}

static void transform_point(const CMFEM_Vector *x,
                            CMFEM_Vector *r,
                            void *context)
{
   const double tol = 1.0e-6;
   const double pi = 3.14159265358979323846;
   double theta = 0.0;
   (void)context;

   if (fabs(CMFEM_Vector_Get(x, 1) + 1.0) < tol)
   {
      theta = 0.25 * pi * (CMFEM_Vector_Get(x, 0) - 2.0);
   }
   else if (fabs(CMFEM_Vector_Get(x, 0) - 1.0) < tol)
   {
      theta = 0.25 * pi * CMFEM_Vector_Get(x, 1);
   }
   else if (fabs(CMFEM_Vector_Get(x, 1) - 1.0) < tol)
   {
      theta = 0.25 * pi * (2.0 - CMFEM_Vector_Get(x, 0));
   }
   else if (fabs(CMFEM_Vector_Get(x, 0) + 1.0) < tol)
   {
      theta = 0.25 * pi * (4.0 - CMFEM_Vector_Get(x, 1));
   }
   else
   {
      fprintf(stderr, "side not recognized %.17g %.17g %.17g\n",
              CMFEM_Vector_Get(x, 0),
              CMFEM_Vector_Get(x, 1),
              CMFEM_Vector_Get(x, 2));
   }

   CMFEM_Vector_Set(r, 0, cos(theta));
   CMFEM_Vector_Set(r, 1, sin(theta));
   CMFEM_Vector_Set(r, 2,
                    0.25 * (2.0 * CMFEM_Vector_Get(x, 2) - 1.0) *
                    (CMFEM_Vector_Get(r, 0) + 2.0));
}

static CMFEM_Mesh *get_mesh(int type)
{
   CMFEM_Mesh *mesh = NULL;

   if (type == 3)
   {
      const double vertices[12][3] =
      {
         {-1.0, -1.0, 0.0}, { 1.0, -1.0, 0.0}, { 1.0,  1.0, 0.0},
         {-1.0,  1.0, 0.0}, {-1.0, -1.0, 1.0}, { 1.0, -1.0, 1.0},
         { 1.0,  1.0, 1.0}, {-1.0,  1.0, 1.0}, { 0.0, -1.0, 0.5},
         { 1.0,  0.0, 0.5}, { 0.0,  1.0, 0.5}, {-1.0,  0.0, 0.5},
      };
      const int triangles[16][3] =
      {
         {0, 1, 8}, {1, 5, 8}, {5, 4, 8}, {4, 0, 8},
         {1, 2, 9}, {2, 6, 9}, {6, 5, 9}, {5, 1, 9},
         {2, 3, 10}, {3, 7, 10}, {7, 6, 10}, {6, 2, 10},
         {3, 0, 11}, {0, 4, 11}, {4, 7, 11}, {7, 3, 11},
      };
      const int bdr_segments[8][3] =
      {
         {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},
         {5, 4, 2}, {6, 5, 2}, {7, 6, 2}, {4, 7, 2},
      };
      int i;

      mesh = CMFEM_Mesh_NewDimensionVerticesElementsBoundarySpace(2, 12, 16, 8, 3);
      for (i = 0; i < 12; i++)
      {
         CMFEM_Mesh_AddVertex(mesh, vertices[i]);
      }
      for (i = 0; i < 16; i++)
      {
         CMFEM_Mesh_AddTriangle(mesh, triangles[i], 1);
      }
      for (i = 0; i < 8; i++)
      {
         CMFEM_Mesh_AddBdrSegment(mesh, bdr_segments[i], bdr_segments[i][2]);
      }
   }
   else if (type == 4)
   {
      const double vertices[8][3] =
      {
         {-1.0, -1.0, 0.0}, { 1.0, -1.0, 0.0}, { 1.0,  1.0, 0.0},
         {-1.0,  1.0, 0.0}, {-1.0, -1.0, 1.0}, { 1.0, -1.0, 1.0},
         { 1.0,  1.0, 1.0}, {-1.0,  1.0, 1.0},
      };
      const int quads[4][4] =
      {
         {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7},
      };
      const int bdr_segments[8][3] =
      {
         {0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 0, 1},
         {5, 4, 2}, {6, 5, 2}, {7, 6, 2}, {4, 7, 2},
      };
      int i;

      mesh = CMFEM_Mesh_NewDimensionVerticesElementsBoundarySpace(2, 8, 4, 8, 3);
      for (i = 0; i < 8; i++)
      {
         CMFEM_Mesh_AddVertex(mesh, vertices[i]);
      }
      for (i = 0; i < 4; i++)
      {
         CMFEM_Mesh_AddQuad(mesh, quads[i], 1);
      }
      for (i = 0; i < 8; i++)
      {
         CMFEM_Mesh_AddBdrSegment(mesh, bdr_segments[i], bdr_segments[i][2]);
      }
   }
   else
   {
      fprintf(stderr, "Unrecognized mesh type %d\n", type);
      return NULL;
   }

   CMFEM_Mesh_FinalizeTopology(mesh);
   return mesh;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   int order = 3;
   int mesh_type = 4;
   int mesh_order = 3;
   int ref_levels = 0;
   int static_cond = 0;
   int visualization = 1;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_int_option(argc, argv, &i, "-mt", "--mesh-type",
                                          &mesh_type);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-mo", "--mesh-order",
                                      &mesh_order);
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

   // 2. Construct a triangular or quadrilateral surface mesh and curve it into
   //    a skewed cylindrical surface embedded in 3D.
   CMFEM_Mesh *mesh = get_mesh(mesh_type);
   if (mesh == NULL)
   {
      return 2;
   }
   {
      const int dim = CMFEM_Mesh_Dimension(mesh);
      CMFEM_H1FeCollection *fec = NULL;
      CMFEM_FiniteElementSpace *fespace = NULL;
      CMFEM_ArrayInt *ess_tdof_list = NULL;
      CMFEM_LinearForm *b = NULL;
      CMFEM_BilinearForm *a = NULL;
      CMFEM_GridFunction *x = NULL;
      CMFEM_OperatorPtr A;
      CMFEM_Vector B;
      CMFEM_Vector X;
      CMFEM_ConstantCoefficient *one = NULL;
      CMFEM_FunctionCoefficient *u_coef = NULL;
      CMFEM_VectorFunctionCoefficient *flux_coef = NULL;
      CMFEM_MatrixFunctionCoefficient *sigma = NULL;
      CMFEM_DiffusionIntegrator *flux_integ = NULL;
      CMFEM_FiniteElementSpace *flux_fespace = NULL;
      CMFEM_GridFunction *flux = NULL;
      double error;
      double flux_error;

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      CMFEM_Mesh_SetCurvature(mesh, mesh_order);
      CMFEM_Mesh_Transform(mesh, transform_point, NULL);

      // 3. Define the H1 finite element space and essential boundary dofs.
      fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, fec);
      printf("Number of finite element unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

      ess_tdof_list = CMFEM_ArrayInt_New();
      if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
      {
         CMFEM_ArrayInt *ess_bdr = CMFEM_ArrayInt_NewSize(
                                      CMFEM_Mesh_BoundaryAttributesMax(mesh));
         CMFEM_ArrayInt_Assign(ess_bdr, 1);
         CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, ess_bdr,
                                                       ess_tdof_list);
         CMFEM_ArrayInt_Delete(ess_bdr);
      }

      // 4. Build the right-hand side and zero initial guess.
      one = CMFEM_ConstantCoefficient_New(1.0);
      b = CMFEM_LinearForm_New(fespace);
      CMFEM_LinearForm_AddDomainIntegratorDliCc(b, one);
      CMFEM_LinearForm_Assemble(b);

      x = CMFEM_GridFunction_New(fespace);
      CMFEM_GridFunction_Assign(x, 0.0);

      // 5. Assemble the anisotropic diffusion bilinear form and linear system.
      a = CMFEM_BilinearForm_New(fespace);
      sigma = CMFEM_MatrixFunctionCoefficient_New(3, sigma_func, NULL);
      flux_integ = CMFEM_DiffusionIntegrator_NewMfc(sigma);
      CMFEM_BilinearForm_AddDomainIntegratorDiMfc(a, sigma);
      if (static_cond)
      {
         CMFEM_BilinearForm_EnableStaticCondensation(a);
      }
      CMFEM_BilinearForm_Assemble(a);

      A = CMFEM_OperatorPtr_Construct();
      B = CMFEM_Vector_Construct();
      X = CMFEM_Vector_Construct();
      CMFEM_BilinearForm_FormLinearSystemOp(a, ess_tdof_list, x, b, &A, &X, &B);
      printf("Size of linear system: %d\n", CMFEM_OperatorPtr_Height(&A));

      // 6. Solve the linear system with PCG and Gauss-Seidel preconditioning.
      {
         CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewOp(&A);
         CMFEM_PCGOpGs(&A, M, &B, &X, 1, 200, 1e-12, 0.0);
         CMFEM_GSSmoother_Delete(M);
      }
      CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);

      // 7. Compute the solution and flux errors against the exact fields.
      u_coef = CMFEM_FunctionCoefficient_New(u_exact, NULL);
      error = CMFEM_GridFunction_ComputeL2ErrorFc(x, u_coef);
      printf("|u - u_h|_2 = %.14g\n", error);

      flux_fespace = CMFEM_FiniteElementSpace_NewMeshH1VDim(mesh, fec, 3);
      flux = CMFEM_GridFunction_New(flux_fespace);
      CMFEM_GridFunction_ComputeFluxDi(x, flux_integ, flux);
      CMFEM_GridFunction_Scale(flux, -1.0);

      flux_coef = CMFEM_VectorFunctionCoefficient_New(3, flux_exact, NULL);
      flux_error = CMFEM_GridFunction_ComputeL2ErrorVfc(flux, flux_coef);
      printf("|f - f_h|_2 = %.14g\n", flux_error);

      // 8. Save the refined mesh and scalar solution.
      CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
      CMFEM_GridFunction_Save(x, "sol.gf", 8);

      // 9. Send the solution and flux by socket to a GLVis server.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
         CMFEM_SendSolutionToGLVis(mesh, flux, "localhost", 19916);
      }

      CMFEM_GridFunction_Delete(flux);
      CMFEM_FiniteElementSpace_Delete(flux_fespace);
      CMFEM_VectorFunctionCoefficient_Delete(flux_coef);
      CMFEM_FunctionCoefficient_Delete(u_coef);
      CMFEM_DiffusionIntegrator_Delete(flux_integ);
      CMFEM_MatrixFunctionCoefficient_Delete(sigma);
      CMFEM_OperatorPtr_Destroy(&A);
      CMFEM_Vector_Destroy(&B);
      CMFEM_Vector_Destroy(&X);
      CMFEM_GridFunction_Delete(x);
      CMFEM_BilinearForm_Delete(a);
      CMFEM_LinearForm_Delete(b);
      CMFEM_ArrayInt_Delete(ess_tdof_list);
      CMFEM_ConstantCoefficient_Delete(one);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_H1FeCollection_Delete(fec);
   }

   CMFEM_Mesh_Delete(mesh);
   return 0;
}
