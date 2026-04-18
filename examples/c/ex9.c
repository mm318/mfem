//                               MFEM C Example 9
//
// Description: This example solves the time-dependent advection equation with
//              a discontinuous Galerkin discretization. The PDE-specific DG
//              operator remains in C, while time stepping goes through the
//              generic CMFEM ODE wrapper layer.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum DGAssemblyMode
{
   DG_ASSEMBLY_LEGACY = 0,
   DG_ASSEMBLY_PARTIAL = 1,
   DG_ASSEMBLY_ELEMENT = 2,
   DG_ASSEMBLY_FULL = 3,
};

typedef struct DGEvolutionContext
{
   int assembly_mode;
   int implicit_variable_type;
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list;
   _Alignas(max_align_t) CMFEM_Vector z;
   _Alignas(max_align_t) CMFEM_Vector b;
   _Alignas(max_align_t) CMFEM_SparseMatrix mmat;
   _Alignas(max_align_t) CMFEM_SparseMatrix kmat;
   CMFEM_BilinearForm *m_form;
   CMFEM_BilinearForm *k_form;
   CMFEM_DSmoother *m_prec_ds;
   CMFEM_OperatorJacobiSmoother *m_prec_ojs;
   CMFEM_SparseMatrix *dg_operator;
   CMFEM_GSSmoother *dg_prec;
   double current_dt;
} DGEvolutionContext;

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

static double inflow_condition(const CMFEM_Vector *x, void *context)
{
   (void)x;
   (void)context;
   return 0.0;
}

static void dg_evolution_reset_implicit(DGEvolutionContext *oper)
{
   if (oper->dg_prec != NULL)
   {
      CMFEM_GSSmoother_Delete(oper->dg_prec);
      oper->dg_prec = NULL;
   }
   if (oper->dg_operator != NULL)
   {
      CMFEM_SparseMatrix_Delete(oper->dg_operator);
      oper->dg_operator = NULL;
   }
   oper->current_dt = -1.0;
}

static void dg_evolution_init(DGEvolutionContext *oper,
                              CMFEM_BilinearForm *m_form,
                              CMFEM_BilinearForm *k_form,
                              const CMFEM_LinearForm *boundary_flow,
                              int size,
                              int assembly_mode)
{
   oper->assembly_mode = assembly_mode;
   oper->implicit_variable_type = CMFEM_IMPLICIT_VARIABLE_SLOPE;
   oper->ess_tdof_list = CMFEM_ArrayInt_Construct();
   oper->z = CMFEM_Vector_ConstructSize(size);
   oper->b = CMFEM_Vector_Construct();
   oper->mmat = CMFEM_SparseMatrix_Construct();
   oper->kmat = CMFEM_SparseMatrix_Construct();
   oper->m_form = m_form;
   oper->k_form = k_form;
   oper->m_prec_ds = NULL;
   oper->m_prec_ojs = NULL;
   oper->dg_operator = NULL;
   oper->dg_prec = NULL;
   oper->current_dt = -1.0;

   CMFEM_LinearForm_CopyToVector(boundary_flow, &oper->b);

   if (assembly_mode == DG_ASSEMBLY_LEGACY)
   {
      CMFEM_BilinearForm_FormSystemMatrixSm(oper->m_form,
                                            &oper->ess_tdof_list,
                                            &oper->mmat);
      CMFEM_BilinearForm_FormSystemMatrixSm(oper->k_form,
                                            &oper->ess_tdof_list,
                                            &oper->kmat);
      oper->m_prec_ds = CMFEM_DSmoother_NewSm(&oper->mmat);
   }
   else
   {
      oper->m_prec_ojs = CMFEM_OperatorJacobiSmoother_NewBf(oper->m_form,
                                                            &oper->ess_tdof_list);
   }
}

static void dg_evolution_destroy(DGEvolutionContext *oper)
{
   dg_evolution_reset_implicit(oper);
   if (oper->m_prec_ojs != NULL)
   {
      CMFEM_OperatorJacobiSmoother_Delete(oper->m_prec_ojs);
   }
   if (oper->m_prec_ds != NULL)
   {
      CMFEM_DSmoother_Delete(oper->m_prec_ds);
   }
   CMFEM_SparseMatrix_Destroy(&oper->kmat);
   CMFEM_SparseMatrix_Destroy(&oper->mmat);
   CMFEM_Vector_Destroy(&oper->b);
   CMFEM_Vector_Destroy(&oper->z);
   CMFEM_ArrayInt_Destroy(&oper->ess_tdof_list);
}

static void dg_evolution_mult(const CMFEM_Vector *x,
                              CMFEM_Vector *y,
                              double time,
                              void *context)
{
   DGEvolutionContext *oper = (DGEvolutionContext *)context;
   (void)time;

   // Compute y = M^{-1}(K x + b).
   CMFEM_BilinearForm_Mult(oper->k_form, x, &oper->z);
   CMFEM_Vector_Add(&oper->z, &oper->b);

   if (oper->assembly_mode == DG_ASSEMBLY_LEGACY)
   {
      CMFEM_PCGSmDs(&oper->mmat,
                    oper->m_prec_ds,
                    &oper->z,
                    y,
                    0,
                    100,
                    1.0e-9,
                    0.0);
   }
   else
   {
      CMFEM_PCGBfOjs(oper->m_form,
                     oper->m_prec_ojs,
                     &oper->z,
                     y,
                     0,
                     100,
                     1.0e-9,
                     0.0);
   }
}

static void dg_evolution_implicit_solve(double dt,
                                        const CMFEM_Vector *x,
                                        CMFEM_Vector *k,
                                        double time,
                                        void *context)
{
   DGEvolutionContext *oper = (DGEvolutionContext *)context;
   (void)time;

   if (oper->assembly_mode != DG_ASSEMBLY_LEGACY)
   {
      fprintf(stderr,
              "Implicit time integration is not supported with partial, element, or full assembly in examples/c/ex9.c.\n");
      abort();
   }

   // Construct the current right-hand side for a stage-state or stage-slope
   // solve, then solve (M - dt K) * k = rhs.
   if (oper->implicit_variable_type == CMFEM_IMPLICIT_VARIABLE_STATE)
   {
      CMFEM_BilinearForm_Mult(oper->m_form, x, &oper->z);
   }
   else
   {
      CMFEM_BilinearForm_Mult(oper->k_form, x, &oper->z);
   }
   CMFEM_Vector_Add(&oper->z, &oper->b);

   if (oper->dg_operator == NULL ||
       fabs(dt - oper->current_dt) > 1.0e-14 * (1.0 + fabs(dt)))
   {
      dg_evolution_reset_implicit(oper);
      oper->dg_operator = CMFEM_AddSmSm(1.0, &oper->mmat, -dt, &oper->kmat);
      oper->dg_prec = CMFEM_GSSmoother_NewSm(oper->dg_operator);
      oper->current_dt = dt;
   }

   CMFEM_GMRESSmGs(oper->dg_operator,
                   oper->dg_prec,
                   &oper->z,
                   k,
                   0,
                   100,
                   50,
                   1.0e-9,
                   0.0);
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("periodic-hexagon.mesh");
   int ref_levels = 2;
   int order = 3;
   int pa = 0;
   int ea = 0;
   int fa = 0;
   const char *device_config = "cpu";
   int ode_solver_type = 4;
   double t_final = 10.0;
   double dt = 0.01;
   int visualization = 1;
   int visit = 0;
   int paraview = 0;
   int binary = 0;
   int vis_steps = 5;
   int solve_implicit_state = 0;
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

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-pa", "--partial-assembly",
                                       "-no-pa", "--no-partial-assembly",
                                       &pa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-ea", "--element-assembly",
                                       "-no-ea", "--no-element-assembly",
                                       &ea);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-fa", "--full-assembly",
                                       "-no-fa", "--no-full-assembly",
                                       &fa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device",
                                         &device_config);
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

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-imp-state", "--implicit-state",
                                       "-imp-slope", "--implicit-slope",
                                       &solve_implicit_state);
      if (parsed == 1) { continue; }

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

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-paraview", "--paraview-datafiles",
                                       "-no-paraview", "--no-paraview-datafiles",
                                       &paraview);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-binary", "--binary-datafiles",
                                       "-ascii", "--ascii-datafiles",
                                       &binary);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-vs",
                                      "--visualization-steps", &vis_steps);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   if (paraview)
   {
      fprintf(stderr,
              "ParaView output is not yet supported in examples/c/ex9.c.\n");
      return 1;
   }
   if (visit && binary)
   {
      fprintf(stderr,
              "Binary VisIt output is not yet supported in examples/c/ex9.c.\n");
      return 1;
   }

   // 2. Configure the device and read the mesh. This example can handle
   //    geometrically periodic meshes.
   CMFEM_Device *device = CMFEM_Device_New(device_config);
   CMFEM_Device_Print(device);

   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   const int dim = CMFEM_Mesh_Dimension(mesh);
   CMFEM_ODESolver *ode_solver = NULL;
   CMFEM_DgFeCollection *fec = NULL;
   CMFEM_FiniteElementSpace *fespace = NULL;
   CMFEM_VectorFunctionCoefficient *velocity = NULL;
   CMFEM_FunctionCoefficient *inflow = NULL;
   CMFEM_FunctionCoefficient *u0 = NULL;
   CMFEM_BilinearForm *m_form = NULL;
   CMFEM_BilinearForm *k_form = NULL;
   CMFEM_LinearForm *boundary_flow = NULL;
   CMFEM_GridFunction *u_gf = NULL;
   _Alignas(max_align_t) CMFEM_Vector u = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector bb_min = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector bb_max = CMFEM_Vector_Construct();
   DGEvolutionContext oper;
   CMFEM_TimeDependentOperator *oper_adapter = NULL;
   CMFEM_VisItDataCollection *visit_dc = NULL;
   const int assembly_mode = pa ? DG_ASSEMBLY_PARTIAL
                             : (ea ? DG_ASSEMBLY_ELEMENT
                                : (fa ? DG_ASSEMBLY_FULL
                                   : DG_ASSEMBLY_LEGACY));

   // 3. Define the ODE solver and uniformly refine the mesh.
   ode_solver = CMFEM_ODESolver_NewType(ode_solver_type);
   if (ode_solver == NULL)
   {
      fprintf(stderr, "Unsupported ODE solver type: %d\n", ode_solver_type);
      CMFEM_Vector_Destroy(&bb_max);
      CMFEM_Vector_Destroy(&bb_min);
      CMFEM_Vector_Destroy(&u);
      CMFEM_Mesh_Delete(mesh);
      CMFEM_Device_Delete(device);
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

   // 4. Define the discontinuous finite element space of the requested order.
   fec = CMFEM_DgFeCollection_NewOrderDimBasis(order,
                                               dim,
                                               CMFEM_BASIS_GAUSS_LOBATTO);
   fespace = CMFEM_FiniteElementSpace_NewMeshDg(mesh, fec);

   printf("Number of unknowns: %d\n", CMFEM_FiniteElementSpace_GetVSize(fespace));

   // 5. Assemble the DG mass, advection, and boundary-flow forms.
   velocity = CMFEM_VectorFunctionCoefficient_New(dim, velocity_function, NULL);
   inflow = CMFEM_FunctionCoefficient_New(inflow_condition, NULL);
   u0 = CMFEM_FunctionCoefficient_New(initial_condition, NULL);
   m_form = CMFEM_BilinearForm_New(fespace);
   k_form = CMFEM_BilinearForm_New(fespace);
   boundary_flow = CMFEM_LinearForm_New(fespace);

   if (assembly_mode == DG_ASSEMBLY_PARTIAL)
   {
      CMFEM_BilinearForm_SetAssemblyLevelPartial(m_form);
      CMFEM_BilinearForm_SetAssemblyLevelPartial(k_form);
   }
   else if (assembly_mode == DG_ASSEMBLY_ELEMENT)
   {
      CMFEM_BilinearForm_SetAssemblyLevelElement(m_form);
      CMFEM_BilinearForm_SetAssemblyLevelElement(k_form);
   }
   else if (assembly_mode == DG_ASSEMBLY_FULL)
   {
      CMFEM_BilinearForm_SetAssemblyLevelFull(m_form);
      CMFEM_BilinearForm_SetAssemblyLevelFull(k_form);
   }

   {
      CMFEM_ConstantCoefficient *one = CMFEM_ConstantCoefficient_New(1.0);
      CMFEM_BilinearForm_AddDomainIntegratorMiCc(m_form, one);
      CMFEM_ConstantCoefficient_Delete(one);
   }
   CMFEM_BilinearForm_AddDomainIntegratorCviVfc(k_form, velocity, -1.0);
   CMFEM_BilinearForm_AddInteriorFaceIntegratorNdtVfc(k_form, velocity, -1.0);
   CMFEM_BilinearForm_AddBdrFaceIntegratorNdtVfc(k_form, velocity, -1.0);
   CMFEM_LinearForm_AddBdrFaceIntegratorBfiFcVfc(boundary_flow,
                                                 inflow,
                                                 velocity,
                                                 -1.0);

   CMFEM_BilinearForm_Assemble(m_form);
   CMFEM_BilinearForm_AssembleSkipZeros(k_form, 0);
   CMFEM_LinearForm_Assemble(boundary_flow);
   CMFEM_BilinearForm_Finalize(m_form);
   CMFEM_BilinearForm_FinalizeSkipZeros(k_form, 0);

   // 6. Project the initial condition, save it, and initialize output.
   u_gf = CMFEM_GridFunction_New(fespace);
   CMFEM_GridFunction_ProjectCoefficientFc(u_gf, u0);
   CMFEM_GridFunction_GetTrueDofs(u_gf, &u);
   CMFEM_Mesh_Print(mesh, "ex9.mesh", precision);
   CMFEM_GridFunction_Save(u_gf, "ex9-init.gf", precision);

   if (visit)
   {
      visit_dc = CMFEM_VisItDataCollection_New("Example9", mesh);
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

   // 7. Build the time-dependent DG evolution operator and choose whether the
   //    implicit callbacks solve for stage states or stage slopes.
   dg_evolution_init(&oper,
                     m_form,
                     k_form,
                     boundary_flow,
                     CMFEM_FiniteElementSpace_GetTrueVSize(fespace),
                     assembly_mode);
   oper_adapter = CMFEM_TimeDependentOperator_New(
                     CMFEM_FiniteElementSpace_GetTrueVSize(fespace),
                     CMFEM_TIME_DEPENDENT_OPERATOR_IMPLICIT,
                     dg_evolution_mult,
                     &oper,
                     dg_evolution_implicit_solve,
                     &oper);
   oper.implicit_variable_type = solve_implicit_state
                                 ? CMFEM_IMPLICIT_VARIABLE_STATE
                                 : CMFEM_IMPLICIT_VARIABLE_SLOPE;
   CMFEM_TimeDependentOperator_SetImplicitVariableType(
      oper_adapter,
      oper.implicit_variable_type);

   // 8. Perform time integration.
   {
      double t = 0.0;
      int ti = 0;
      int done = 0;

      CMFEM_ODESolver_Init(ode_solver, oper_adapter);

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

   // 9. Save the final solution for later visualization.
   CMFEM_GridFunction_SetFromTrueDofs(u_gf, &u);
   CMFEM_GridFunction_Save(u_gf, "ex9-final.gf", precision);

   // 10. Free the used memory.
   if (visit_dc != NULL)
   {
      CMFEM_VisItDataCollection_Delete(visit_dc);
   }
   CMFEM_TimeDependentOperator_Delete(oper_adapter);
   dg_evolution_destroy(&oper);
   CMFEM_GridFunction_Delete(u_gf);
   CMFEM_LinearForm_Delete(boundary_flow);
   CMFEM_BilinearForm_Delete(k_form);
   CMFEM_BilinearForm_Delete(m_form);
   CMFEM_FunctionCoefficient_Delete(u0);
   CMFEM_FunctionCoefficient_Delete(inflow);
   CMFEM_VectorFunctionCoefficient_Delete(velocity);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_DgFeCollection_Delete(fec);
   CMFEM_ODESolver_Delete(ode_solver);
   CMFEM_Vector_Destroy(&bb_max);
   CMFEM_Vector_Destroy(&bb_min);
   CMFEM_Vector_Destroy(&u);
   CMFEM_Mesh_Delete(mesh);
   CMFEM_Device_Delete(device);

   return 0;
}
