//                               MFEM C Example 15
//
// Description: Building on Example 6, this example demonstrates dynamic AMR.
//              The mesh is adapted to a time-dependent solution by refinement
//              as well as by derefinement.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

static int g_problem = 0;
static int g_nfeatures = 1;
static const double alpha = 0.02;
static const double pi_value = 3.14159265358979323846;

static double front(double x, double y, double z, double t, int dim)
{
   const double r = sqrt(x * x + y * y + z * z);
   (void)dim;
   return exp(-0.5 * pow((r - t) / alpha, 2.0));
}

static double front_laplace(double x, double y, double z, double t, int dim)
{
   const double x2 = x * x;
   const double y2 = y * y;
   const double z2 = z * z;
   const double t2 = t * t;
   const double r = sqrt(x2 + y2 + z2);
   const double a2 = alpha * alpha;
   const double a4 = a2 * a2;
   return -exp(-0.5 * pow((r - t) / alpha, 2.0)) / a4 *
          (-2.0 * t * (x2 + y2 + z2 - (dim - 1) * a2 / 2.0) / r +
           x2 + y2 + z2 + t2 - dim * a2);
}

static double ball(double x, double y, double z, double t, int dim)
{
   const double r = sqrt(x * x + y * y + z * z);
   (void)dim;
   return -atan(2.0 * (r - t) / alpha);
}

static double ball_laplace(double x, double y, double z, double t, int dim)
{
   const double x2 = x * x;
   const double y2 = y * y;
   const double z2 = z * z;
   const double t2 = 4.0 * t * t;
   const double r = sqrt(x2 + y2 + z2);
   const double a2 = alpha * alpha;
   const double den = pow(-a2 - 4.0 * (x2 + y2 + z2 - 2.0 * r * t) - t2, 2.0);
   if (dim == 2)
   {
      return 2.0 * alpha * (a2 + t2 - 4.0 * x2 - 4.0 * y2) / r / den;
   }
   return 4.0 * alpha * (a2 + t2 - 4.0 * r * t) / r / den;
}

static double composite_func(const CMFEM_Vector *pt,
                             double time,
                             double (*f0)(double, double, double, double, int),
                             double (*f1)(double, double, double, double, int))
{
   const int dim = CMFEM_Vector_Size(pt);
   const double x = CMFEM_Vector_Get(pt, 0);
   const double y = CMFEM_Vector_Get(pt, 1);
   const double z = dim == 3 ? CMFEM_Vector_Get(pt, 2) : 0.0;
   int i;

   if (g_problem == 0)
   {
      if (g_nfeatures <= 1)
      {
         return f0(x, y, z, time, dim);
      }

      {
         double sum = 0.0;
         for (i = 0; i < g_nfeatures; i++)
         {
            const double x0 = 0.5 * cos(2.0 * pi_value * i / g_nfeatures);
            const double y0 = 0.5 * sin(2.0 * pi_value * i / g_nfeatures);
            sum += f0(x - x0, y - y0, z, time, dim);
         }
         return sum;
      }
   }

   {
      double sum = 0.0;
      for (i = 0; i < g_nfeatures; i++)
      {
         const double x0 =
            0.5 * cos(2.0 * pi_value * i / g_nfeatures + pi_value * time);
         const double y0 =
            0.5 * sin(2.0 * pi_value * i / g_nfeatures + pi_value * time);
         sum += f1(x - x0, y - y0, z, 0.25, dim);
      }
      return sum;
   }
}

static double bdr_func(const CMFEM_Vector *pt, double time, void *context)
{
   (void)context;
   return composite_func(pt, time, front, ball);
}

static double rhs_func(const CMFEM_Vector *pt, double time, void *context)
{
   (void)context;
   return composite_func(pt, time, front_laplace, ball_laplace);
}

static void update_problem(CMFEM_FiniteElementSpace *fespace,
                           CMFEM_GridFunction *x,
                           CMFEM_BilinearForm *a,
                           CMFEM_LinearForm *b)
{
   CMFEM_FiniteElementSpace_Update(fespace);
   CMFEM_GridFunction_Update(x);
   CMFEM_FiniteElementSpace_UpdatesFinished(fespace);
   CMFEM_BilinearForm_Update(a);
   CMFEM_LinearForm_Update(b);
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   g_problem = 0;
   g_nfeatures = 1;
   const char *mesh_file = CMFEM_ExamplesDataPath("star-hilbert.mesh");
   int order = 2;
   double t_final = 1.0;
   double max_elem_error = 5.0e-3;
   double hysteresis = 0.15;
   int ref_levels = 0;
   int nc_limit = 3;
   int visualization = 1;
   int visit = 0;
   int which_estimator = 0;
   int dim;
   int sdim;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-p", "--problem",
                                      &g_problem);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-n", "--nfeatures",
                                      &g_nfeatures);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-e", "--max-err",
                                         &max_elem_error);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-y", "--hysteresis",
                                         &hysteresis);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-r", "--ref-levels",
                                      &ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-l", "--nc-limit",
                                      &nc_limit);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-tf", "--t-final",
                                         &t_final);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-est", "--estimator",
                                      &which_estimator);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-visit", "--visit-datafiles",
                                       "-no-visit", "--no-visit-datafiles",
                                       &visit);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Read and prepare the mesh.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   dim = CMFEM_Mesh_Dimension(mesh);
   sdim = CMFEM_Mesh_SpaceDimension(mesh);

   if (CMFEM_Mesh_HasNURBSext(mesh))
   {
      CMFEM_Mesh_UniformRefinement(mesh);
      if (ref_levels > 0) { ref_levels--; }
      CMFEM_Mesh_SetCurvature(mesh, 2);
   }
   CMFEM_Mesh_EnsureNCMesh(mesh, 1);
   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }
   CMFEM_Mesh_Finalize(mesh, 1, 0);

   // 3. Mark all boundary attributes as essential.
   if (CMFEM_Mesh_BoundaryAttributesSize(mesh) <= 0)
   {
      fprintf(stderr, "Boundary attributes required in the mesh.\n");
      CMFEM_Mesh_Delete(mesh);
      return 1;
   }
   _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
      CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
   CMFEM_ArrayInt_Assign(&ess_bdr, 1);

   // 4. Define the finite element space and the Poisson forms.
   CMFEM_H1_FECollection *fec = CMFEM_H1_FECollection_NewOrderDim(order, dim);
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh,
                                                                          fec);
   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);

   CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_FunctionCoefficient *bdr =
      CMFEM_FunctionCoefficient_NewTimeDependent(bdr_func, NULL);
   CMFEM_FunctionCoefficient *rhs =
      CMFEM_FunctionCoefficient_NewTimeDependent(rhs_func, NULL);
   CMFEM_BilinearForm_AddDomainIntegrator_DiffusionCoefficient(a, one);
   CMFEM_LinearForm_AddDomainIntegrator_DomainLFIntegrator_FunctionCoefficient(
      b, rhs);

   // 5. Maintain the solution over AMR iterations.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_Assign(x, 0.0);

   // 6. Optional VisIt output.
   CMFEM_VisItDataCollection *visit_dc = NULL;
   int vis_cycle = 0;
   if (visit)
   {
      visit_dc = CMFEM_VisItDataCollection_New("Example15", mesh);
      CMFEM_VisItDataCollection_RegisterFieldGridFunction(visit_dc, "solution",
                                                          x);
   }

   // 7. Set up the error estimator and refinement operators.
   CMFEM_DiffusionIntegrator *diff_integ =
      CMFEM_DiffusionIntegrator_NewConstantCoefficient(one);
   CMFEM_ZienkiewiczZhuEstimator *zz_estimator = NULL;
   CMFEM_KellyErrorEstimator *kelly_estimator = NULL;
   CMFEM_ThresholdRefiner *refiner = NULL;
   CMFEM_ThresholdDerefiner *derefiner = NULL;
   CMFEM_L2_FECollection *l2_flux_fec = NULL;

   if (which_estimator == 1)
   {
      l2_flux_fec = CMFEM_L2_FECollection_NewOrderDim(order, dim);
      CMFEM_FiniteElementSpace *flux_fes =
         CMFEM_FiniteElementSpace_NewMeshL2VDim(mesh, l2_flux_fec, sdim);
      kelly_estimator =
         CMFEM_KellyErrorEstimator_NewDiffusionIntegratorGridFunctionFESpace(
            diff_integ, x, flux_fes);
      refiner = CMFEM_ThresholdRefiner_NewKellyErrorEstimator(kelly_estimator);
      derefiner = CMFEM_ThresholdDerefiner_NewKellyErrorEstimator(
         kelly_estimator);
   }
   else
   {
      CMFEM_FiniteElementSpace *flux_fes =
         CMFEM_FiniteElementSpace_NewMeshH1VDim(mesh, fec, sdim);
      zz_estimator =
         CMFEM_ZienkiewiczZhuEstimator_NewDiffusionIntegratorGridFunctionFESpace(
            diff_integ, x, flux_fes);
      refiner = CMFEM_ThresholdRefiner_NewZienkiewiczZhuEstimator(zz_estimator);
      derefiner = CMFEM_ThresholdDerefiner_NewZienkiewiczZhuEstimator(
         zz_estimator);
   }

   CMFEM_ThresholdRefiner_SetTotalErrorFraction(refiner, 0.0);
   CMFEM_ThresholdRefiner_SetLocalErrorGoal(refiner, max_elem_error);
   CMFEM_ThresholdRefiner_PreferConformingRefinement(refiner);
   CMFEM_ThresholdRefiner_SetNCLimit(refiner, nc_limit);
   CMFEM_ThresholdDerefiner_SetThreshold(derefiner, hysteresis * max_elem_error);
   CMFEM_ThresholdDerefiner_SetNCLimit(derefiner, nc_limit);

   // 8. The outer time loop.
   {
      double time;
      CMFEM_GridFunction_Assign(x, 0.0);
      for (time = 0.0; time < t_final + 1.0e-10; time += 0.01)
      {
         printf("\nTime %.14g\n\nRefinement:\n", time);

         CMFEM_FunctionCoefficient_SetTime(bdr, time);
         CMFEM_FunctionCoefficient_SetTime(rhs, time);

         CMFEM_ThresholdRefiner_Reset(refiner);
         CMFEM_ThresholdDerefiner_Reset(derefiner);

         // 9. Inner refinement loop.
         {
            int ref_it;
            for (ref_it = 1; ; ref_it++)
            {
               printf("Iteration: %d, number of unknowns: %d\n",
                      ref_it,
                      CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

               CMFEM_BilinearForm_Assemble(a);
               CMFEM_LinearForm_Assemble(b);
               CMFEM_GridFunction_ProjectBdrCoefficient_FunctionCoefficient(
                  x, bdr, &ess_bdr);

               _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list =
                  CMFEM_ArrayInt_Construct();
               CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr,
                                                             &ess_tdof_list);

               _Alignas(max_align_t) CMFEM_SparseMatrix A =
                  CMFEM_SparseMatrix_Construct();
               _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
               _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
               CMFEM_BilinearForm_FormLinearSystemSparseMatrix(a,
                                                               &ess_tdof_list,
                                                               x,
                                                               b,
                                                               &A,
                                                               &X,
                                                               &B);

               {
                  CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewSparseMatrix(&A);
                  CMFEM_PCG_SparseMatrixGSSmoother(&A, M, &B, &X, 0, 500,
                                                   1.0e-12, 0.0);
                  CMFEM_GSSmoother_Delete(M);
               }

               CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);

               if (visualization)
               {
                  CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
               }
               if (visit)
               {
                  CMFEM_VisItDataCollection_SetCycle(visit_dc, vis_cycle++);
                  CMFEM_VisItDataCollection_SetTime(visit_dc, time);
                  CMFEM_VisItDataCollection_Save(visit_dc);
               }

               CMFEM_SparseMatrix_Destroy(&A);
               CMFEM_Vector_Destroy(&B);
               CMFEM_Vector_Destroy(&X);
               CMFEM_ArrayInt_Destroy(&ess_tdof_list);

               if (!CMFEM_ThresholdRefiner_Apply(refiner, mesh) ||
                   CMFEM_ThresholdRefiner_Stop(refiner))
               {
                  break;
               }

               update_problem(fespace, x, a, b);
            }
         }

         // 10. Derefine if indicated.
         if (CMFEM_ThresholdDerefiner_Apply(derefiner, mesh))
         {
            printf("\nDerefined elements.\n");
            update_problem(fespace, x, a, b);
         }

         CMFEM_BilinearForm_Update(a);
         CMFEM_LinearForm_Update(b);
      }
   }

   if (visit_dc)
   {
      CMFEM_VisItDataCollection_Delete(visit_dc);
   }
   CMFEM_ThresholdDerefiner_Delete(derefiner);
   CMFEM_ThresholdRefiner_Delete(refiner);
   if (zz_estimator) { CMFEM_ZienkiewiczZhuEstimator_Delete(zz_estimator); }
   if (kelly_estimator) { CMFEM_KellyErrorEstimator_Delete(kelly_estimator); }
   if (l2_flux_fec) { CMFEM_L2_FECollection_Delete(l2_flux_fec); }
   CMFEM_DiffusionIntegrator_Delete(diff_integ);
   CMFEM_GridFunction_Delete(x);
   CMFEM_FunctionCoefficient_Delete(rhs);
   CMFEM_FunctionCoefficient_Delete(bdr);
   CMFEM_ConstantCoefficient_Delete(one);
   CMFEM_LinearForm_Delete(b);
   CMFEM_BilinearForm_Delete(a);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_H1_FECollection_Delete(fec);
   CMFEM_ArrayInt_Destroy(&ess_bdr);
   CMFEM_Mesh_Delete(mesh);
   return 0;
}
