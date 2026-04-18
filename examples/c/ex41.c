//                               MFEM C Example 41
//
// Description: This example solves the time-dependent advection-diffusion
//              equation with an IMEX time integrator. The PDE-specific split
//              operator remains in C, while the IMEX solver is accessed through
//              the generic CMFEM wrapper layer.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct IMEXEvolutionContext
{
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list;
   _Alignas(max_align_t) CMFEM_Vector z;
   _Alignas(max_align_t) CMFEM_SparseMatrix mmat;
   _Alignas(max_align_t) CMFEM_SparseMatrix smat;
   CMFEM_BilinearForm *m_form;
   CMFEM_BilinearForm *k_form;
   CMFEM_BilinearForm *s_form;
   CMFEM_DSmoother *m_prec;
   CMFEM_SparseMatrix *implicit_operator;
   CMFEM_DSmoother *implicit_prec;
   double current_dt;
} IMEXEvolutionContext;

static int problem_ = 0;
static double bb_min_[3] = { 0.0, 0.0, 0.0 };
static double bb_max_[3] = { 1.0, 1.0, 1.0 };
static const double cmfem_pi_ = 3.141592653589793238462643383279502884;

static void map_to_reference_box(const CMFEM_Vector *x, double *X, int dim)
{
   int i;
   for (i = 0; i < dim; i++)
   {
      const double center = 0.5 * (bb_min_[i] + bb_max_[i]);
      X[i] = 2.0 * (CMFEM_Vector_Get(x, i) - center) / (bb_max_[i] - bb_min_[i]);
   }
}

static void velocity_function(const CMFEM_Vector *x,
                              CMFEM_Vector *v,
                              void *context)
{
   int d;
   const int dim = CMFEM_Vector_Size(x);
   double X[3] = { 0.0, 0.0, 0.0 };
   (void)context;

   map_to_reference_box(x, X, dim);
   for (d = 0; d < dim; d++)
   {
      CMFEM_Vector_Set(v, d, 0.0);
   }

   switch (problem_)
   {
      case 0:
      {
         switch (dim)
         {
            case 1:
               CMFEM_Vector_Set(v, 0, 1.0);
               break;
            case 2:
               CMFEM_Vector_Set(v, 0, sqrt(2.0 / 3.0));
               CMFEM_Vector_Set(v, 1, sqrt(1.0 / 3.0));
               break;
            case 3:
               CMFEM_Vector_Set(v, 0, sqrt(3.0 / 6.0));
               CMFEM_Vector_Set(v, 1, sqrt(2.0 / 6.0));
               CMFEM_Vector_Set(v, 2, sqrt(1.0 / 6.0));
               break;
         }
         break;
      }
      case 1:
      case 2:
      {
         const double w = cmfem_pi_ / 2.0;
         switch (dim)
         {
            case 1:
               CMFEM_Vector_Set(v, 0, 1.0);
               break;
            case 2:
               CMFEM_Vector_Set(v, 0, w * X[1]);
               CMFEM_Vector_Set(v, 1, -w * X[0]);
               break;
            case 3:
               CMFEM_Vector_Set(v, 0, w * X[1]);
               CMFEM_Vector_Set(v, 1, -w * X[0]);
               CMFEM_Vector_Set(v, 2, 0.0);
               break;
         }
         break;
      }
      case 3:
      {
         const double w = cmfem_pi_ / 2.0;
         double d_factor = fmax((X[0] + 1.0) * (1.0 - X[0]), 0.0) *
                           fmax((X[1] + 1.0) * (1.0 - X[1]), 0.0);
         d_factor *= d_factor;
         switch (dim)
         {
            case 1:
               CMFEM_Vector_Set(v, 0, 1.0);
               break;
            case 2:
               CMFEM_Vector_Set(v, 0, d_factor * w * X[1]);
               CMFEM_Vector_Set(v, 1, -d_factor * w * X[0]);
               break;
            case 3:
               CMFEM_Vector_Set(v, 0, d_factor * w * X[1]);
               CMFEM_Vector_Set(v, 1, -d_factor * w * X[0]);
               CMFEM_Vector_Set(v, 2, 0.0);
               break;
         }
         break;
      }
   }
}

static double initial_condition(const CMFEM_Vector *x, void *context)
{
   const int dim = CMFEM_Vector_Size(x);
   double X[3] = { 0.0, 0.0, 0.0 };
   (void)context;

   map_to_reference_box(x, X, dim);

   switch (problem_)
   {
      case 0:
      case 1:
      {
         switch (dim)
         {
            case 1:
               return exp(-40.0 * pow(X[0] - 0.5, 2.0));
            case 2:
            case 3:
            {
               double rx = 0.45;
               double ry = 0.25;
               const double cx = 0.0;
               const double cy = -0.2;
               const double w = 10.0;
               if (dim == 3)
               {
                  const double scale = 1.0 + 0.25 * cos(2.0 * cmfem_pi_ * X[2]);
                  rx *= scale;
                  ry *= scale;
               }
               return (erfc(w * (X[0] - cx - rx)) *
                       erfc(-w * (X[0] - cx + rx)) *
                       erfc(w * (X[1] - cy - ry)) *
                       erfc(-w * (X[1] - cy + ry))) / 16.0;
            }
         }
         break;
      }
      case 2:
      {
         const double rho = hypot(X[0], X[1]);
         const double phi = atan2(X[1], X[0]);
         return pow(sin(cmfem_pi_ * rho), 2.0) * sin(3.0 * phi);
      }
      case 3:
         return sin(cmfem_pi_ * X[0]) * sin(cmfem_pi_ * X[1]);
   }

   return 0.0;
}

static void imex_evolution_reset_implicit(IMEXEvolutionContext *oper)
{
   if (oper->implicit_prec != NULL)
   {
      CMFEM_DSmoother_Delete(oper->implicit_prec);
      oper->implicit_prec = NULL;
   }
   if (oper->implicit_operator != NULL)
   {
      CMFEM_SparseMatrix_Delete(oper->implicit_operator);
      oper->implicit_operator = NULL;
   }
   oper->current_dt = -1.0;
}

static void imex_evolution_init(IMEXEvolutionContext *oper,
                                CMFEM_BilinearForm *m_form,
                                CMFEM_BilinearForm *k_form,
                                CMFEM_BilinearForm *s_form,
                                int size)
{
   oper->ess_tdof_list = CMFEM_ArrayInt_Construct();
   oper->z = CMFEM_Vector_ConstructSize(size);
   oper->mmat = CMFEM_SparseMatrix_Construct();
   oper->smat = CMFEM_SparseMatrix_Construct();
   oper->m_form = m_form;
   oper->k_form = k_form;
   oper->s_form = s_form;
   oper->m_prec = NULL;
   oper->implicit_operator = NULL;
   oper->implicit_prec = NULL;
   oper->current_dt = -1.0;

   CMFEM_BilinearForm_FormSystemMatrixSm(oper->m_form,
                                         &oper->ess_tdof_list,
                                         &oper->mmat);
   CMFEM_BilinearForm_FormSystemMatrixSm(oper->s_form,
                                         &oper->ess_tdof_list,
                                         &oper->smat);
   oper->m_prec = CMFEM_DSmoother_NewSm(&oper->mmat);
}

static void imex_evolution_destroy(IMEXEvolutionContext *oper)
{
   imex_evolution_reset_implicit(oper);
   if (oper->m_prec != NULL)
   {
      CMFEM_DSmoother_Delete(oper->m_prec);
   }
   CMFEM_SparseMatrix_Destroy(&oper->smat);
   CMFEM_SparseMatrix_Destroy(&oper->mmat);
   CMFEM_Vector_Destroy(&oper->z);
   CMFEM_ArrayInt_Destroy(&oper->ess_tdof_list);
}

static void imex_mult_term_1(const CMFEM_Vector *x,
                             CMFEM_Vector *y,
                             double time,
                             void *context)
{
   IMEXEvolutionContext *oper = (IMEXEvolutionContext *)context;
   (void)time;

   // Evaluate y = M^{-1}(K x).
   CMFEM_BilinearForm_Mult(oper->k_form, x, &oper->z);
   CMFEM_PCGSmDs(&oper->mmat,
                 oper->m_prec,
                 &oper->z,
                 y,
                 0,
                 100,
                 1.0e-9,
                 0.0);
}

static void imex_implicit_solve_term_2(double dt,
                                       const CMFEM_Vector *x,
                                       CMFEM_Vector *k,
                                       double time,
                                       void *context)
{
   IMEXEvolutionContext *oper = (IMEXEvolutionContext *)context;
   const double tolerance = 1.0e-14 * (1.0 + fabs(dt));
   (void)time;

   // Solve (M + dt S) k = -S x.
   CMFEM_BilinearForm_Mult(oper->s_form, x, &oper->z);
   CMFEM_Vector_Neg(&oper->z);

   if (oper->implicit_operator == NULL ||
       fabs(dt - oper->current_dt) > tolerance)
   {
      imex_evolution_reset_implicit(oper);
      oper->implicit_operator = CMFEM_AddSmSm(1.0, &oper->mmat, dt, &oper->smat);
      oper->implicit_prec = CMFEM_DSmoother_NewSm(oper->implicit_operator);
      oper->current_dt = dt;
   }

   CMFEM_PCGSmDs(oper->implicit_operator,
                 oper->implicit_prec,
                 &oper->z,
                 k,
                 0,
                 100,
                 1.0e-9,
                 0.0);
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("periodic-square.mesh");
   int ref_levels = 2;
   int order = 3;
   int ode_solver_type = 64;
   double t_final = 10.0;
   double dt = 0.01;
   int paraview = 0;
   int cg = 0;
   int vis_steps = 50;
   double diffusion_term = 0.01;
   double kappa = -1.0;
   double sigma = -1.0;
   int visualization = 1;
   int visit = 0;
   int binary = 0;
   const int precision = 8;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-p", "--problem",
                                      &problem_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-r", "--refine",
                                      &ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-s", "--ode-solver",
                                      &ode_solver_type);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-tf", "--t-final",
                                         &t_final);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-dt", "--time-step",
                                         &dt);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-dc",
                                         "--diffusion-coeff", &diffusion_term);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-paraview", "--paraview-datafiles",
                                       "-no-paraview", "--no-paraview-datafiles",
                                       &paraview);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-vs",
                                      "--visualization-steps", &vis_steps);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-binary", "--binary-datafiles",
                                       "-ascii", "--ascii-datafiles",
                                       &binary);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-visit", "--visit-datafiles",
                                       "-no-visit", "--no-visit-datafiles",
                                       &visit);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-cg", "--continuous-galerkin",
                                       "-dg", "--discontinuous-galerkin",
                                       &cg);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   if (kappa < 0.0)
   {
      kappa = (double)((order + 1) * (order + 1));
   }
   if (paraview)
   {
      fprintf(stderr,
              "ParaView output is not yet supported in examples/c/ex41.c.\n");
      return 1;
   }
   if (visit && binary)
   {
      fprintf(stderr,
              "Binary VisIt output is not yet supported in examples/c/ex41.c.\n");
      return 1;
   }

   // 2. Read and uniformly refine the mesh.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   const int dim = CMFEM_Mesh_Dimension(mesh);
   CMFEM_H1FeCollection *h1_fec = NULL;
   CMFEM_DgFeCollection *dg_fec = NULL;
   CMFEM_FiniteElementSpace *fespace = NULL;
   CMFEM_VectorFunctionCoefficient *velocity = NULL;
   CMFEM_ConstantCoefficient *diff_coeff = NULL;
   CMFEM_ConstantCoefficient *one = NULL;
   CMFEM_BilinearForm *m_form = NULL;
   CMFEM_BilinearForm *k_form = NULL;
   CMFEM_BilinearForm *s_form = NULL;
   CMFEM_FunctionCoefficient *u0 = NULL;
   CMFEM_GridFunction *u_gf = NULL;
   _Alignas(max_align_t) CMFEM_Vector u = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector bb_min = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector bb_max = CMFEM_Vector_Construct();
   IMEXEvolutionContext oper;
   CMFEM_TimeDependentOperator *oper_adapter = NULL;
   CMFEM_ODESolver *ode_solver = NULL;
   CMFEM_VisItDataCollection *visit_dc = NULL;

   ode_solver = CMFEM_ODESolver_NewImexType(ode_solver_type);
   if (ode_solver == NULL)
   {
      fprintf(stderr, "Unsupported IMEX ODE solver type: %d\n", ode_solver_type);
      CMFEM_Vector_Destroy(&bb_max);
      CMFEM_Vector_Destroy(&bb_min);
      CMFEM_Vector_Destroy(&u);
      CMFEM_Mesh_Delete(mesh);
      return 1;
   }

   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }
   if (CMFEM_Mesh_HasNURBSext(mesh))
   {
      CMFEM_Mesh_SetCurvature(mesh, order > 1 ? order : 1);
   }
   CMFEM_Mesh_GetBoundingBox(mesh, &bb_min, &bb_max, order > 1 ? order : 1);
   for (i = 0; i < dim && i < 3; i++)
   {
      bb_min_[i] = CMFEM_Vector_Get(&bb_min, i);
      bb_max_[i] = CMFEM_Vector_Get(&bb_max, i);
   }

   // 3. Define the continuous or discontinuous finite element space.
   if (cg)
   {
      h1_fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, h1_fec);
   }
   else
   {
      dg_fec = CMFEM_DgFeCollection_NewOrderDimBasis(order,
                                                     dim,
                                                     CMFEM_BASIS_GAUSS_LOBATTO);
      fespace = CMFEM_FiniteElementSpace_NewMeshDg(mesh, dg_fec);
   }

   printf("Number of unknowns: %d\n", CMFEM_FiniteElementSpace_GetVSize(fespace));

   // 4. Assemble the mass, advection, and diffusion forms.
   velocity = CMFEM_VectorFunctionCoefficient_New(dim, velocity_function, NULL);
   diff_coeff = CMFEM_ConstantCoefficient_New(diffusion_term);
   one = CMFEM_ConstantCoefficient_New(1.0);
   m_form = CMFEM_BilinearForm_New(fespace);
   k_form = CMFEM_BilinearForm_New(fespace);
   s_form = CMFEM_BilinearForm_New(fespace);

   CMFEM_BilinearForm_AddDomainIntegratorMiCc(m_form, one);
   CMFEM_BilinearForm_AddDomainIntegratorCviVfc(k_form, velocity, -1.0);
   CMFEM_BilinearForm_AddDomainIntegratorDiCc(s_form, diff_coeff);
   if (!cg)
   {
      CMFEM_BilinearForm_AddInteriorFaceIntegratorNdtVfc(k_form, velocity, -1.0);
      CMFEM_BilinearForm_AddBdrFaceIntegratorNdtVfc(k_form, velocity, -1.0);
      CMFEM_BilinearForm_AddInteriorFaceIntegratorDgd(s_form,
                                                      diff_coeff,
                                                      sigma,
                                                      kappa);
      CMFEM_BilinearForm_AddBdrFaceIntegratorDgd(s_form,
                                                 diff_coeff,
                                                 sigma,
                                                 kappa);
   }

   CMFEM_BilinearForm_AssembleSkipZeros(m_form, 0);
   CMFEM_BilinearForm_AssembleSkipZeros(k_form, 0);
   CMFEM_BilinearForm_AssembleSkipZeros(s_form, 0);
   CMFEM_BilinearForm_FinalizeSkipZeros(m_form, 0);
   CMFEM_BilinearForm_FinalizeSkipZeros(k_form, 0);
   CMFEM_BilinearForm_FinalizeSkipZeros(s_form, 0);

   // 5. Define and project the initial condition.
   u0 = CMFEM_FunctionCoefficient_New(initial_condition, NULL);
   u_gf = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_ProjectCoefficientFc(u_gf, u0);
   CMFEM_GridFunction_GetTrueDofs(u_gf, &u);

   if (visit)
   {
      visit_dc = CMFEM_VisItDataCollection_New("Example41", mesh);
      CMFEM_VisItDataCollection_RegisterFieldGf(visit_dc, "solution", u_gf);
      CMFEM_VisItDataCollection_SetPrecision(visit_dc, precision);
      CMFEM_VisItDataCollection_SetCycle(visit_dc, 0);
      CMFEM_VisItDataCollection_SetTime(visit_dc, 0.0);
      CMFEM_VisItDataCollection_Save(visit_dc);
   }

   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, u_gf, "localhost", 19916);
   }

   // 6. Define the IMEX evolution operator and initialize the IMEX ODE solver.
   imex_evolution_init(&oper,
                       m_form,
                       k_form,
                       s_form,
                       CMFEM_FiniteElementSpace_GetTrueVSize(fespace));
   oper_adapter = CMFEM_TimeDependentOperator_NewImex(
                     CMFEM_FiniteElementSpace_GetTrueVSize(fespace),
                     imex_mult_term_1,
                     &oper,
                     imex_implicit_solve_term_2,
                     &oper);
   CMFEM_ODESolver_Init(ode_solver, oper_adapter);

   // 7. Perform time integration.
   {
      double t = 0.0;
      int ti = 0;
      int done = 0;

      while (!done)
      {
         double dt_real = dt;
         if (t + dt_real > t_final)
         {
            dt_real = t_final - t;
         }
         CMFEM_ODESolver_Step(ode_solver, &u, &t, &dt_real);
         ti++;

         done = (t >= t_final - 1.0e-8 * dt);
         if (done || (vis_steps > 0 && (ti % vis_steps) == 0))
         {
            printf("time step: %d, time: %.17g\n", ti, t);
            CMFEM_GridFunction_SetFromTrueDofs(u_gf, &u);

            if (visualization)
            {
               CMFEM_SendSolutionToGLVis(mesh, u_gf, "localhost", 19916);
            }
            if (visit)
            {
               CMFEM_VisItDataCollection_SetCycle(visit_dc, ti);
               CMFEM_VisItDataCollection_SetTime(visit_dc, t);
               CMFEM_VisItDataCollection_Save(visit_dc);
            }
         }
      }
   }

   // 8. Free the used memory.
   if (visit_dc != NULL)
   {
      CMFEM_VisItDataCollection_Delete(visit_dc);
   }
   CMFEM_TimeDependentOperator_Delete(oper_adapter);
   imex_evolution_destroy(&oper);
   CMFEM_GridFunction_Delete(u_gf);
   CMFEM_FunctionCoefficient_Delete(u0);
   CMFEM_BilinearForm_Delete(s_form);
   CMFEM_BilinearForm_Delete(k_form);
   CMFEM_BilinearForm_Delete(m_form);
   CMFEM_ConstantCoefficient_Delete(one);
   CMFEM_ConstantCoefficient_Delete(diff_coeff);
   CMFEM_VectorFunctionCoefficient_Delete(velocity);
   CMFEM_FiniteElementSpace_Delete(fespace);
   if (dg_fec != NULL)
   {
      CMFEM_DgFeCollection_Delete(dg_fec);
   }
   if (h1_fec != NULL)
   {
      CMFEM_H1FeCollection_Delete(h1_fec);
   }
   CMFEM_ODESolver_Delete(ode_solver);
   CMFEM_Vector_Destroy(&bb_max);
   CMFEM_Vector_Destroy(&bb_min);
   CMFEM_Vector_Destroy(&u);
   CMFEM_Mesh_Delete(mesh);
   return 0;
}
