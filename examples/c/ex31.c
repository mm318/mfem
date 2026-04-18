//                              MFEM C Example 31
//
// Description: This example solves a restricted-space electromagnetic
//              diffusion problem in 1D, 2D, or 3D using Nedelec elements.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

static double freq_ = 1.0;
static double kappa_ = 0.0;
static int dim_ = 0;
static const double cmfem_pi_ = 3.14159265358979323846;
static const double cmfem_sqrt1_2_ = 0.70710678118654752440;

static void e_exact(const CMFEM_Vector *x, CMFEM_Vector *e, void *context)
{
   (void)context;
   if (dim_ == 1)
   {
      CMFEM_Vector_Set(e, 0, 1.1 * sin(kappa_ * CMFEM_Vector_Get(x, 0)));
      CMFEM_Vector_Set(e, 1, 1.2 * sin(kappa_ * CMFEM_Vector_Get(x, 0) +
                                       0.4 * cmfem_pi_));
      CMFEM_Vector_Set(e, 2, 1.3 * sin(kappa_ * CMFEM_Vector_Get(x, 0) +
                                       0.9 * cmfem_pi_));
   }
   else if (dim_ == 2)
   {
      const double phase = kappa_ * cmfem_sqrt1_2_ *
                           (CMFEM_Vector_Get(x, 0) + CMFEM_Vector_Get(x, 1));
      CMFEM_Vector_Set(e, 0, 1.1 * sin(phase));
      CMFEM_Vector_Set(e, 1, 1.2 * sin(phase + 0.4 * cmfem_pi_));
      CMFEM_Vector_Set(e, 2, 1.3 * sin(phase + 0.9 * cmfem_pi_));
   }
   else
   {
      const double phase = kappa_ * cmfem_sqrt1_2_ *
                           (CMFEM_Vector_Get(x, 0) + CMFEM_Vector_Get(x, 1));
      const double factor = cos(kappa_ * CMFEM_Vector_Get(x, 2));
      CMFEM_Vector_Set(e, 0, 1.1 * sin(phase) * factor);
      CMFEM_Vector_Set(e, 1, 1.2 * sin(phase + 0.4 * cmfem_pi_) * factor);
      CMFEM_Vector_Set(e, 2, 1.3 * sin(phase + 0.9 * cmfem_pi_) * factor);
   }
}

static void curl_e_exact(const CMFEM_Vector *x, CMFEM_Vector *de, void *context)
{
   (void)context;
   if (dim_ == 1)
   {
      const double c4 = cos(kappa_ * CMFEM_Vector_Get(x, 0) + 0.4 * cmfem_pi_);
      const double c9 = cos(kappa_ * CMFEM_Vector_Get(x, 0) + 0.9 * cmfem_pi_);

      CMFEM_Vector_Set(de, 0, 0.0);
      CMFEM_Vector_Set(de, 1, -1.3 * c9 * kappa_);
      CMFEM_Vector_Set(de, 2, 1.2 * c4 * kappa_);
   }
   else if (dim_ == 2)
   {
      const double phase = kappa_ * cmfem_sqrt1_2_ *
                           (CMFEM_Vector_Get(x, 0) + CMFEM_Vector_Get(x, 1));
      const double c0 = cos(phase);
      const double c4 = cos(phase + 0.4 * cmfem_pi_);
      const double c9 = cos(phase + 0.9 * cmfem_pi_);
      const double scale = kappa_ * cmfem_sqrt1_2_;

      CMFEM_Vector_Set(de, 0, 1.3 * c9 * scale);
      CMFEM_Vector_Set(de, 1, -1.3 * c9 * scale);
      CMFEM_Vector_Set(de, 2, (1.2 * c4 - 1.1 * c0) * scale);
   }
   else
   {
      const double phase = kappa_ * cmfem_sqrt1_2_ *
                           (CMFEM_Vector_Get(x, 0) + CMFEM_Vector_Get(x, 1));
      const double s0 = sin(phase);
      const double c0 = cos(phase);
      const double s4 = sin(phase + 0.4 * cmfem_pi_);
      const double c4 = cos(phase + 0.4 * cmfem_pi_);
      const double c9 = cos(phase + 0.9 * cmfem_pi_);
      const double sk = sin(kappa_ * CMFEM_Vector_Get(x, 2));
      const double ck = cos(kappa_ * CMFEM_Vector_Get(x, 2));

      CMFEM_Vector_Set(de, 0,
                       kappa_ * (1.2 * s4 * sk +
                                 1.3 * cmfem_sqrt1_2_ * c9 * ck));
      CMFEM_Vector_Set(de, 1,
                       kappa_ * (-1.1 * s0 * sk -
                                 1.3 * cmfem_sqrt1_2_ * c9 * ck));
      CMFEM_Vector_Set(de, 2,
                       -kappa_ * cmfem_sqrt1_2_ * (1.1 * c0 - 1.2 * c4) * ck);
   }
}

static void f_exact(const CMFEM_Vector *x, CMFEM_Vector *f, void *context)
{
   (void)context;
   if (dim_ == 1)
   {
      const double s0 = sin(kappa_ * CMFEM_Vector_Get(x, 0));
      const double s4 = sin(kappa_ * CMFEM_Vector_Get(x, 0) + 0.4 * cmfem_pi_);
      const double s9 = sin(kappa_ * CMFEM_Vector_Get(x, 0) + 0.9 * cmfem_pi_);

      CMFEM_Vector_Set(f, 0, 2.2 * s0 + 1.2 * cmfem_sqrt1_2_ * s4);
      CMFEM_Vector_Set(f, 1,
                       1.2 * (2.0 + kappa_ * kappa_) * s4 +
                       cmfem_sqrt1_2_ * (1.1 * s0 + 1.3 * s9));
      CMFEM_Vector_Set(f, 2,
                       1.3 * (2.0 + kappa_ * kappa_) * s9 +
                       1.2 * cmfem_sqrt1_2_ * s4);
   }
   else if (dim_ == 2)
   {
      const double phase = kappa_ * cmfem_sqrt1_2_ *
                           (CMFEM_Vector_Get(x, 0) + CMFEM_Vector_Get(x, 1));
      const double s0 = sin(phase);
      const double s4 = sin(phase + 0.4 * cmfem_pi_);
      const double s9 = sin(phase + 0.9 * cmfem_pi_);

      CMFEM_Vector_Set(f, 0,
                       0.55 * (4.0 + kappa_ * kappa_) * s0 +
                       0.6 * (1.4142135623730950488 - kappa_ * kappa_) * s4);
      CMFEM_Vector_Set(f, 1,
                       0.55 * (1.4142135623730950488 - kappa_ * kappa_) * s0 +
                       0.6 * (4.0 + kappa_ * kappa_) * s4 +
                       0.65 * 1.4142135623730950488 * s9);
      CMFEM_Vector_Set(f, 2,
                       0.6 * 1.4142135623730950488 * s4 +
                       1.3 * (2.0 + kappa_ * kappa_) * s9);
   }
   else
   {
      const double phase = kappa_ * cmfem_sqrt1_2_ *
                           (CMFEM_Vector_Get(x, 0) + CMFEM_Vector_Get(x, 1));
      const double s0 = sin(phase);
      const double c0 = cos(phase);
      const double s4 = sin(phase + 0.4 * cmfem_pi_);
      const double c4 = cos(phase + 0.4 * cmfem_pi_);
      const double s9 = sin(phase + 0.9 * cmfem_pi_);
      const double c9 = cos(phase + 0.9 * cmfem_pi_);
      const double sk = sin(kappa_ * CMFEM_Vector_Get(x, 2));
      const double ck = cos(kappa_ * CMFEM_Vector_Get(x, 2));

      CMFEM_Vector_Set(f, 0,
                       0.55 * (4.0 + 3.0 * kappa_ * kappa_) * s0 * ck +
                       0.6 * (1.4142135623730950488 - kappa_ * kappa_) * s4 * ck -
                       0.65 * 1.4142135623730950488 * kappa_ * kappa_ * c9 * sk);
      CMFEM_Vector_Set(f, 1,
                       0.55 * (1.4142135623730950488 - kappa_ * kappa_) * s0 * ck +
                       0.6 * (4.0 + 3.0 * kappa_ * kappa_) * s4 * ck +
                       0.65 * 1.4142135623730950488 * s9 * ck -
                       0.65 * 1.4142135623730950488 * kappa_ * kappa_ * c9 * sk);
      CMFEM_Vector_Set(f, 2,
                       0.6 * 1.4142135623730950488 * s4 * ck -
                       1.4142135623730950488 * kappa_ * kappa_ *
                       (0.55 * c0 + 0.6 * c4) * sk +
                       1.3 * (2.0 + kappa_ * kappa_) * s9 * ck);
   }
}

static void sigma_func(const CMFEM_Vector *x,
                       double *values,
                       int dim,
                       void *context)
{
   int i;
   (void)x;
   (void)context;
   for (i = 0; i < dim * dim; i++)
   {
      values[i] = 0.0;
   }
   values[0 * dim + 0] = 2.0;
   values[1 * dim + 1] = 2.0;
   values[2 * dim + 2] = 2.0;
   values[0 * dim + 1] = cmfem_sqrt1_2_;
   values[1 * dim + 0] = cmfem_sqrt1_2_;
   values[1 * dim + 2] = cmfem_sqrt1_2_;
   values[2 * dim + 1] = cmfem_sqrt1_2_;
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("inline-quad.mesh");
   int ref_levels = 2;
   int order = 1;
   int visualization = 1;
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

      parsed = cmfem_parse_double_option(argc, argv, &i, "-f", "--frequency",
                                         &freq_);
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
   kappa_ = freq_ * cmfem_pi_;

   // 2. Read and refine the mesh.
   {
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      dim_ = CMFEM_Mesh_Dimension(mesh);
      CMFEM_H1FeCollection *h1_fec = NULL;
      CMFEM_NdFeCollection *nd_fec = NULL;
      CMFEM_NdR1dFeCollection *ndr1d_fec = NULL;
      CMFEM_NdR2dFeCollection *ndr2d_fec = NULL;
      CMFEM_FiniteElementSpace *fespace = NULL;
      CMFEM_GridFunction *sol = NULL;
      CMFEM_LinearForm *b = NULL;
      CMFEM_BilinearForm *a = NULL;
      CMFEM_ConstantCoefficient *muinv = NULL;
      CMFEM_MatrixFunctionCoefficient *sigma = NULL;
      CMFEM_VectorFunctionCoefficient *f = NULL;
      CMFEM_VectorFunctionCoefficient *e = NULL;
      CMFEM_VectorFunctionCoefficient *curl_e = NULL;
      CMFEM_ArrayInt *ess_tdof_list = NULL;
      CMFEM_OperatorPtr A;
      CMFEM_Vector B;
      CMFEM_Vector X;
      double error;

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      // 3. Define the restricted Nedelec space matching the mesh dimension.
      if (dim_ == 1)
      {
         ndr1d_fec = CMFEM_NdR1dFeCollection_NewOrderDim(order, dim_);
         fespace = CMFEM_FiniteElementSpace_NewMeshNdR1d(mesh, ndr1d_fec);
      }
      else if (dim_ == 2)
      {
         ndr2d_fec = CMFEM_NdR2dFeCollection_NewOrderDim(order, dim_);
         fespace = CMFEM_FiniteElementSpace_NewMeshNdR2d(mesh, ndr2d_fec);
      }
      else
      {
         nd_fec = CMFEM_NdFeCollection_NewOrderDim(order, dim_);
         fespace = CMFEM_FiniteElementSpace_NewMeshNd(mesh, nd_fec);
      }
      printf("Number of H(Curl) unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

      // 4. Mark all boundary attributes as essential.
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

      // 5. Build the right-hand side, initial boundary projection, and the
      //    electromagnetic diffusion bilinear form.
      f = CMFEM_VectorFunctionCoefficient_New(3, f_exact, NULL);
      b = CMFEM_LinearForm_New(fespace);
      CMFEM_LinearForm_AddDomainIntegratorVfd(b, f);
      CMFEM_LinearForm_Assemble(b);

      sol = CMFEM_GridFunction_New(fespace);
      e = CMFEM_VectorFunctionCoefficient_New(3, e_exact, NULL);
      curl_e = CMFEM_VectorFunctionCoefficient_New(3, curl_e_exact, NULL);
      CMFEM_GridFunction_ProjectCoefficientVfc(sol, e);

      muinv = CMFEM_ConstantCoefficient_New(1.0);
      sigma = CMFEM_MatrixFunctionCoefficient_New(3, sigma_func, NULL);
      a = CMFEM_BilinearForm_New(fespace);
      CMFEM_BilinearForm_AddDomainIntegratorCci(a, muinv);
      CMFEM_BilinearForm_AddDomainIntegratorVmiMfc(a, sigma);
      CMFEM_BilinearForm_Assemble(a);

      // 6. Assemble and solve the linear system.
      A = CMFEM_OperatorPtr_Construct();
      B = CMFEM_Vector_Construct();
      X = CMFEM_Vector_Construct();
      CMFEM_BilinearForm_FormLinearSystemOp(a, ess_tdof_list, sol, b, &A, &X, &B);
      {
         CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewOp(&A);
         CMFEM_PCGOpGs(&A, M, &B, &X, 1, 500, 1e-12, 0.0);
         CMFEM_GSSmoother_Delete(M);
      }
      CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, sol);

      // 7. Compute and print the H(Curl) error.
      error = CMFEM_GridFunction_ComputeHCurlErrorVfcVfc(sol, e, curl_e);
      printf("\n|| E_h - E ||_{H(Curl)} = %.14g\n\n", error);

      // 8. Save the refined mesh and the solution.
      CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
      CMFEM_GridFunction_Save(sol, "sol.gf", 8);

      // 9. Send the solution to GLVis. The full C++ visualization split by
      //    components is not mirrored here yet.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, sol, "localhost", 19916);
      }

      CMFEM_OperatorPtr_Destroy(&A);
      CMFEM_Vector_Destroy(&B);
      CMFEM_Vector_Destroy(&X);
      CMFEM_ArrayInt_Delete(ess_tdof_list);
      CMFEM_VectorFunctionCoefficient_Delete(curl_e);
      CMFEM_VectorFunctionCoefficient_Delete(e);
      CMFEM_VectorFunctionCoefficient_Delete(f);
      CMFEM_MatrixFunctionCoefficient_Delete(sigma);
      CMFEM_ConstantCoefficient_Delete(muinv);
      CMFEM_BilinearForm_Delete(a);
      CMFEM_GridFunction_Delete(sol);
      CMFEM_LinearForm_Delete(b);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_NdFeCollection_Delete(nd_fec);
      CMFEM_NdR1dFeCollection_Delete(ndr1d_fec);
      CMFEM_NdR2dFeCollection_Delete(ndr2d_fec);
      CMFEM_H1FeCollection_Delete(h1_fec);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
