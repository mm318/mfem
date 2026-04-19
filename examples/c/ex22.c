//                              MFEM C Example 22
//
// Description: This example solves simple complex-valued damped harmonic
//              oscillator problems in H1, H(curl), or H(div) using MFEM's
//              complex sesquilinear forms.

#include "common.h"
#include "cmfem.h"

#include <complex.h>
#include <stdio.h>
#include <string.h>

static double mu_ = 1.0;
static double epsilon_ = 1.0;
static double sigma_ = 20.0;
static double omega_ = 10.0;

typedef struct DefaultVectorContext
{
   int dim;
   int component;
   double value;
} DefaultVectorContext;

static int check_for_inline_mesh(const char *mesh_file)
{
   const char *name = strrchr(mesh_file, '/');
   if (name)
   {
      name++;
   }
   else
   {
      name = mesh_file;
   }
   return strncmp(name, "inline-", 7) == 0;
}

static double complex u0_exact(const CMFEM_Vector *x)
{
   const int dim = CMFEM_Vector_Size(x);
   const double complex alpha = epsilon_ * omega_ - I * sigma_;
   const double complex kappa = csqrt(mu_ * omega_ * alpha);
   return cexp(-I * kappa * CMFEM_Vector_Get(x, dim - 1));
}

static double u0_real_exact(const CMFEM_Vector *x, void *context)
{
   (void)context;
   return creal(u0_exact(x));
}

static double u0_imag_exact(const CMFEM_Vector *x, void *context)
{
   (void)context;
   return cimag(u0_exact(x));
}

static double scalar_one(const CMFEM_Vector *x, void *context)
{
   (void)x;
   (void)context;
   return 1.0;
}

static double scalar_zero(const CMFEM_Vector *x, void *context)
{
   (void)x;
   (void)context;
   return 0.0;
}

static void default_vector_value(const CMFEM_Vector *x, CMFEM_Vector *value,
                                 void *context)
{
   const DefaultVectorContext *ctx = (const DefaultVectorContext *)context;
   int i;
   (void)x;
   for (i = 0; i < ctx->dim; i++)
   {
      CMFEM_Vector_Set(value, i, 0.0);
   }
   CMFEM_Vector_Set(value, ctx->component, ctx->value);
}

static void u1_real_exact(const CMFEM_Vector *x, CMFEM_Vector *value,
                          void *context)
{
   const int dim = CMFEM_Vector_Size(x);
   int i;
   (void)context;
   for (i = 0; i < dim; i++)
   {
      CMFEM_Vector_Set(value, i, 0.0);
   }
   CMFEM_Vector_Set(value, 0, u0_real_exact(x, NULL));
}

static void u1_imag_exact(const CMFEM_Vector *x, CMFEM_Vector *value,
                          void *context)
{
   const int dim = CMFEM_Vector_Size(x);
   int i;
   (void)context;
   for (i = 0; i < dim; i++)
   {
      CMFEM_Vector_Set(value, i, 0.0);
   }
   CMFEM_Vector_Set(value, 0, u0_imag_exact(x, NULL));
}

static void u2_real_exact(const CMFEM_Vector *x, CMFEM_Vector *value,
                          void *context)
{
   const int dim = CMFEM_Vector_Size(x);
   int i;
   (void)context;
   for (i = 0; i < dim; i++)
   {
      CMFEM_Vector_Set(value, i, 0.0);
   }
   CMFEM_Vector_Set(value, dim - 1, u0_real_exact(x, NULL));
}

static void u2_imag_exact(const CMFEM_Vector *x, CMFEM_Vector *value,
                          void *context)
{
   const int dim = CMFEM_Vector_Size(x);
   int i;
   (void)context;
   for (i = 0; i < dim; i++)
   {
      CMFEM_Vector_Set(value, i, 0.0);
   }
   CMFEM_Vector_Set(value, dim - 1, u0_imag_exact(x, NULL));
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("inline-quad.mesh");
   int ref_levels = 0;
   int order = 1;
   int prob = 0;
   double freq = -1.0;
   double a_coef = 0.0;
   int visualization = 1;
   int herm_conv = 1;
   int exact_sol = 1;
   int pa = 0;
   const char *device_config = "cpu";
   enum CMFEM_ComplexConvention convention = CMFEM_ComplexConventionHermitian;
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

      parsed = cmfem_parse_int_option(argc, argv, &i, "-p", "--problem-type",
                                      &prob);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-a",
                                         "--stiffness-coef", &a_coef);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-b", "--mass-coef",
                                         &epsilon_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-c",
                                         "--damping-coef", &sigma_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-mu",
                                         "--permeability", &mu_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-eps",
                                         "--permittivity", &epsilon_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-sigma",
                                         "--conductivity", &sigma_);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-f", "--frequency",
                                         &freq);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-herm", "--hermitian",
                                       "-no-herm", "--no-hermitian",
                                       &herm_conv);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-pa", "--partial-assembly",
                                       "-no-pa", "--no-partial-assembly",
                                       &pa);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device",
                                         &device_config);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   if (prob < 0 || prob > 2)
   {
      fprintf(stderr, "Unrecognized problem type: %d\n", prob);
      return 1;
   }

   if (a_coef != 0.0)
   {
      mu_ = 1.0 / a_coef;
   }
   if (freq > 0.0)
   {
      omega_ = 2.0 * 3.14159265358979323846264338327950288 * freq;
   }

   exact_sol = check_for_inline_mesh(mesh_file);
   if (exact_sol)
   {
      printf("Identified a mesh with known exact solution\n");
   }
   convention = herm_conv ? CMFEM_ComplexConventionHermitian :
                CMFEM_ComplexConventionBlockSymmetric;

   // 2. Enable hardware devices and read/refine the mesh.
   {
      CMFEM_Device *device = CMFEM_Device_New(device_config);
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      const int dim = CMFEM_Mesh_Dimension(mesh);
      void *fec = NULL;
      CMFEM_FiniteElementSpace *fespace = NULL;
      _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr = CMFEM_ArrayInt_Construct();
      CMFEM_ComplexLinearForm *b = NULL;
      CMFEM_ComplexGridFunction *u = NULL;
      CMFEM_ComplexGridFunction *u_exact = NULL;
      CMFEM_FunctionCoefficient *u0_r = NULL;
      CMFEM_FunctionCoefficient *u0_i = NULL;
      CMFEM_VectorFunctionCoefficient *u1_r = NULL;
      CMFEM_VectorFunctionCoefficient *u1_i = NULL;
      CMFEM_VectorFunctionCoefficient *u2_r = NULL;
      CMFEM_VectorFunctionCoefficient *u2_i = NULL;
      CMFEM_ConstantCoefficient *stiffness_coef = NULL;
      CMFEM_ConstantCoefficient *mass_coef = NULL;
      CMFEM_ConstantCoefficient *loss_coef = NULL;
      CMFEM_ConstantCoefficient *neg_mass_coef = NULL;
      CMFEM_FunctionCoefficient *one_coef = NULL;
      CMFEM_FunctionCoefficient *zero_coef = NULL;
      DefaultVectorContext one_vec_ctx;
      DefaultVectorContext zero_vec_ctx;
      CMFEM_VectorFunctionCoefficient *one_vec_coef = NULL;
      CMFEM_VectorFunctionCoefficient *zero_vec_coef = NULL;
      CMFEM_SesquilinearForm *a = NULL;
      CMFEM_BilinearForm *pc_op = NULL;
      _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
      _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector U = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt block_offsets = CMFEM_ArrayInt_Construct();
      CMFEM_BlockDiagonalPreconditioner *bdp = NULL;
      CMFEM_GSSmoother *pc_gs = NULL;
      CMFEM_DSmoother *pc_ds = NULL;
      CMFEM_OperatorJacobiSmoother *pc_ojs = NULL;
      CMFEM_ScaledOperator *pc_i = NULL;
      _Alignas(max_align_t) CMFEM_SparseMatrix pc_matrix = {0};
      int have_pc_matrix = 0;
      CMFEM_GridFunction *u_real = NULL;
      CMFEM_GridFunction *u_imag = NULL;

      CMFEM_Device_Print(device);

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      if (dim == 1 && prob != 0)
      {
         printf("Switching to problem type 0, H1 basis functions, for 1 dimensional mesh.\n");
         prob = 0;
      }

      switch (prob)
      {
         case 0:
            fec = CMFEM_H1FeCollection_NewOrderDim(order, dim);
            fespace = CMFEM_FiniteElementSpace_NewMeshH1(mesh, fec);
            break;
         case 1:
            fec = CMFEM_NdFeCollection_NewOrderDim(order, dim);
            fespace = CMFEM_FiniteElementSpace_NewMeshNd(mesh, fec);
            break;
         case 2:
            fec = CMFEM_RtFeCollection_NewOrderDim(order - 1, dim);
            fespace = CMFEM_FiniteElementSpace_NewMeshRt(mesh, fec);
            break;
      }

      printf("Number of finite element unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

      if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
      {
         CMFEM_ArrayInt_Destroy(&ess_bdr);
         ess_bdr = CMFEM_ArrayInt_ConstructSize(
                      CMFEM_Mesh_BoundaryAttributesMax(mesh));
         CMFEM_ArrayInt_Assign(&ess_bdr, 1);
         CMFEM_FiniteElementSpace_GetEssentialTrueDofs(
            fespace,
            &ess_bdr,
            &ess_tdof_list);
      }

      b = CMFEM_ComplexLinearForm_NewFesConv(fespace, convention);
      CMFEM_ComplexLinearForm_Assign(b, 0.0, 0.0);

      u = CMFEM_ComplexGridFunction_New(fespace);
      CMFEM_ComplexGridFunction_Assign(u, 0.0, 0.0);
      if (exact_sol)
      {
         u_exact = CMFEM_ComplexGridFunction_New(fespace);
         CMFEM_ComplexGridFunction_Assign(u_exact, 0.0, 0.0);
      }

      u0_r = CMFEM_FunctionCoefficient_New(u0_real_exact, NULL);
      u0_i = CMFEM_FunctionCoefficient_New(u0_imag_exact, NULL);
      u1_r = CMFEM_VectorFunctionCoefficient_New(dim, u1_real_exact, NULL);
      u1_i = CMFEM_VectorFunctionCoefficient_New(dim, u1_imag_exact, NULL);
      u2_r = CMFEM_VectorFunctionCoefficient_New(dim, u2_real_exact, NULL);
      u2_i = CMFEM_VectorFunctionCoefficient_New(dim, u2_imag_exact, NULL);

      zero_coef = CMFEM_FunctionCoefficient_New(scalar_zero, NULL);
      one_coef = CMFEM_FunctionCoefficient_New(scalar_one, NULL);
      one_vec_ctx.dim = dim;
      one_vec_ctx.component = prob == 2 ? dim - 1 : 0;
      one_vec_ctx.value = 1.0;
      zero_vec_ctx.dim = dim;
      zero_vec_ctx.component = 0;
      zero_vec_ctx.value = 0.0;
      one_vec_coef = CMFEM_VectorFunctionCoefficient_New(
                        dim, default_vector_value, &one_vec_ctx);
      zero_vec_coef = CMFEM_VectorFunctionCoefficient_New(
                         dim, default_vector_value, &zero_vec_ctx);

      switch (prob)
      {
         case 0:
            if (exact_sol)
            {
               CMFEM_ComplexGridFunction_ProjectBdrCoefficientFcFcAi(
                  u, u0_r, u0_i, &ess_bdr);
               CMFEM_ComplexGridFunction_ProjectCoefficientFcFc(u_exact,
                                                                u0_r, u0_i);
            }
            else
            {
               CMFEM_ComplexGridFunction_ProjectBdrCoefficientFcFcAi(
                  u, one_coef, zero_coef, &ess_bdr);
            }
            break;
         case 1:
            if (exact_sol)
            {
               CMFEM_ComplexGridFunction_ProjectBdrCoefficientTangentVfcVfcAi(
                  u, u1_r, u1_i, &ess_bdr);
               CMFEM_ComplexGridFunction_ProjectCoefficientVfcVfc(u_exact,
                                                                  u1_r, u1_i);
            }
            else
            {
               CMFEM_ComplexGridFunction_ProjectBdrCoefficientTangentVfcVfcAi(
                  u, one_vec_coef, zero_vec_coef, &ess_bdr);
            }
            break;
         case 2:
            if (exact_sol)
            {
               CMFEM_ComplexGridFunction_ProjectBdrCoefficientNormalVfcVfcAi(
                  u, u2_r, u2_i, &ess_bdr);
               CMFEM_ComplexGridFunction_ProjectCoefficientVfcVfc(u_exact,
                                                                  u2_r, u2_i);
            }
            else
            {
               CMFEM_ComplexGridFunction_ProjectBdrCoefficientNormalVfcVfcAi(
                  u, one_vec_coef, zero_vec_coef, &ess_bdr);
            }
            break;
      }

      stiffness_coef = CMFEM_ConstantCoefficient_New(1.0 / mu_);
      mass_coef = CMFEM_ConstantCoefficient_New(-omega_ * omega_ * epsilon_);
      loss_coef = CMFEM_ConstantCoefficient_New(omega_ * sigma_);
      neg_mass_coef = CMFEM_ConstantCoefficient_New(omega_ * omega_ * epsilon_);

      a = CMFEM_SesquilinearForm_NewFesConv(fespace, convention);
      pc_op = CMFEM_BilinearForm_New(fespace);
      if (pa)
      {
         CMFEM_SesquilinearForm_SetAssemblyLevelPartial(a);
         CMFEM_BilinearForm_SetAssemblyLevelPartial(pc_op);
      }

      switch (prob)
      {
         case 0:
            CMFEM_SesquilinearForm_AddDomainIntegratorDiCcNull(
               a, stiffness_coef);
            CMFEM_SesquilinearForm_AddDomainIntegratorMiCcMiCc(
               a, mass_coef, loss_coef);
            CMFEM_BilinearForm_AddDomainIntegratorDiCc(pc_op, stiffness_coef);
            CMFEM_BilinearForm_AddDomainIntegratorMiCc(pc_op, mass_coef);
            CMFEM_BilinearForm_AddDomainIntegratorMiCc(pc_op, loss_coef);
            break;
         case 1:
            CMFEM_SesquilinearForm_AddDomainIntegratorCciCcNull(
               a, stiffness_coef);
            CMFEM_SesquilinearForm_AddDomainIntegratorVmiCcVmiCc(
               a, mass_coef, loss_coef);
            CMFEM_BilinearForm_AddDomainIntegratorCci(pc_op, stiffness_coef);
            CMFEM_BilinearForm_AddDomainIntegratorVmi(pc_op, neg_mass_coef);
            CMFEM_BilinearForm_AddDomainIntegratorVmi(pc_op, loss_coef);
            break;
         case 2:
            CMFEM_SesquilinearForm_AddDomainIntegratorDdiCcNull(
               a, stiffness_coef);
            CMFEM_SesquilinearForm_AddDomainIntegratorVmiCcVmiCc(
               a, mass_coef, loss_coef);
            CMFEM_BilinearForm_AddDomainIntegratorDdi(pc_op, stiffness_coef);
            CMFEM_BilinearForm_AddDomainIntegratorVmi(pc_op, mass_coef);
            CMFEM_BilinearForm_AddDomainIntegratorVmi(pc_op, loss_coef);
            break;
      }

      CMFEM_SesquilinearForm_Assemble(a);
      CMFEM_BilinearForm_Assemble(pc_op);
      CMFEM_SesquilinearForm_FormLinearSystemOp(
         a, &ess_tdof_list, u, b, &A, &U, &B);

      printf("Size of linear system: %d\n\n", CMFEM_OperatorPtr_Height(&A));

      CMFEM_ArrayInt_Destroy(&block_offsets);
      block_offsets = CMFEM_ArrayInt_ConstructSize(3);
      CMFEM_ArrayInt_Set(&block_offsets, 0, 0);
      CMFEM_ArrayInt_Set(&block_offsets, 1, CMFEM_OperatorPtr_Height(&A) / 2);
      CMFEM_ArrayInt_Set(&block_offsets, 2, CMFEM_OperatorPtr_Height(&A));
      bdp = CMFEM_BlockDiagonalPreconditioner_NewAi(&block_offsets);

      if (pa)
      {
         pc_ojs = CMFEM_OperatorJacobiSmoother_NewBf(pc_op, &ess_tdof_list);
         pc_i = CMFEM_ScaledOperator_NewOjs(
                   pc_ojs,
                   (convention == CMFEM_ComplexConventionHermitian) ?
                   (prob != 1 ? 1.0 : -1.0) :
                   (prob != 1 ? -1.0 : 1.0));
         CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockOjs(bdp, 0, pc_ojs);
         CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockSop(bdp, 1, pc_i);
      }
      else
      {
         CMFEM_BilinearForm_SetDiagonalPolicyOne(pc_op);
         pc_matrix = CMFEM_SparseMatrix_Construct();
         have_pc_matrix = 1;
         CMFEM_BilinearForm_FormSystemMatrixSm(pc_op, &ess_tdof_list, &pc_matrix);
         if (prob == 1)
         {
            pc_gs = CMFEM_GSSmoother_NewSm(&pc_matrix);
            pc_i = CMFEM_ScaledOperator_NewGs(
                      pc_gs,
                      convention == CMFEM_ComplexConventionHermitian ? -1.0 : 1.0);
            CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockGs(bdp, 0, pc_gs);
         }
         else
         {
            pc_ds = CMFEM_DSmoother_NewSm(&pc_matrix);
            pc_i = CMFEM_ScaledOperator_NewDs(
                      pc_ds,
                      convention == CMFEM_ComplexConventionHermitian ? 1.0 : -1.0);
            CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockDs(bdp, 0, pc_ds);
         }
         CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockSop(bdp, 1, pc_i);
      }

      CMFEM_GMRESOpBdp(&A, bdp, &B, &U, 1, 1000, 200, 1.0e-12, 0.0);
      CMFEM_SesquilinearForm_RecoverFEMSolution(a, &U, b, u);

      if (exact_sol)
      {
         double err_r = -1.0;
         double err_i = -1.0;
         switch (prob)
         {
            case 0:
               err_r = CMFEM_ComplexGridFunction_ComputeL2ErrorRealFc(u, u0_r);
               err_i = CMFEM_ComplexGridFunction_ComputeL2ErrorImagFc(u, u0_i);
               break;
            case 1:
               err_r = CMFEM_ComplexGridFunction_ComputeL2ErrorRealVfc(u, u1_r);
               err_i = CMFEM_ComplexGridFunction_ComputeL2ErrorImagVfc(u, u1_i);
               break;
            case 2:
               err_r = CMFEM_ComplexGridFunction_ComputeL2ErrorRealVfc(u, u2_r);
               err_i = CMFEM_ComplexGridFunction_ComputeL2ErrorImagVfc(u, u2_i);
               break;
         }
         printf("\n|| Re (u_h - u) ||_{L^2} = %.8g\n", err_r);
         printf("|| Im (u_h - u) ||_{L^2} = %.8g\n\n", err_i);
      }

      CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
      CMFEM_ComplexGridFunction_SaveReal(u, "sol_r.gf", 8);
      CMFEM_ComplexGridFunction_SaveImag(u, "sol_i.gf", 8);
      CMFEM_ComplexGridFunction_Save(u, "sol_z.gf", 8);

      if (visualization)
      {
         u_real = CMFEM_GridFunction_New(fespace);
         u_imag = CMFEM_GridFunction_New(fespace);
         CMFEM_ComplexGridFunction_CopyRealToGf(u, u_real);
         CMFEM_ComplexGridFunction_CopyImagToGf(u, u_imag);
         CMFEM_SendSolutionToGLVis(mesh, u_real, "localhost", 19916);
         CMFEM_SendSolutionToGLVis(mesh, u_imag, "localhost", 19916);
      }

      CMFEM_GridFunction_Delete(u_imag);
      CMFEM_GridFunction_Delete(u_real);
      if (have_pc_matrix)
      {
         CMFEM_SparseMatrix_Destroy(&pc_matrix);
      }
      CMFEM_ScaledOperator_Delete(pc_i);
      CMFEM_GSSmoother_Delete(pc_gs);
      CMFEM_DSmoother_Delete(pc_ds);
      CMFEM_OperatorJacobiSmoother_Delete(pc_ojs);
      CMFEM_BlockDiagonalPreconditioner_Delete(bdp);
      CMFEM_OperatorPtr_Destroy(&A);
      CMFEM_Vector_Destroy(&B);
      CMFEM_Vector_Destroy(&U);
      CMFEM_ArrayInt_Destroy(&block_offsets);
      CMFEM_SesquilinearForm_Delete(a);
      CMFEM_BilinearForm_Delete(pc_op);
      CMFEM_FunctionCoefficient_Delete(one_coef);
      CMFEM_FunctionCoefficient_Delete(zero_coef);
      CMFEM_ConstantCoefficient_Delete(stiffness_coef);
      CMFEM_ConstantCoefficient_Delete(mass_coef);
      CMFEM_ConstantCoefficient_Delete(loss_coef);
      CMFEM_ConstantCoefficient_Delete(neg_mass_coef);
      CMFEM_VectorFunctionCoefficient_Delete(one_vec_coef);
      CMFEM_VectorFunctionCoefficient_Delete(zero_vec_coef);
      CMFEM_FunctionCoefficient_Delete(u0_r);
      CMFEM_FunctionCoefficient_Delete(u0_i);
      CMFEM_VectorFunctionCoefficient_Delete(u1_r);
      CMFEM_VectorFunctionCoefficient_Delete(u1_i);
      CMFEM_VectorFunctionCoefficient_Delete(u2_r);
      CMFEM_VectorFunctionCoefficient_Delete(u2_i);
      CMFEM_ComplexGridFunction_Delete(u_exact);
      CMFEM_ComplexGridFunction_Delete(u);
      CMFEM_ComplexLinearForm_Delete(b);
      CMFEM_ArrayInt_Destroy(&ess_tdof_list);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_FiniteElementSpace_Delete(fespace);
      switch (prob)
      {
         case 0:
            CMFEM_H1FeCollection_Delete(fec);
            break;
         case 1:
            CMFEM_NdFeCollection_Delete(fec);
            break;
         case 2:
            CMFEM_RtFeCollection_Delete(fec);
            break;
      }
      CMFEM_Mesh_Delete(mesh);
      CMFEM_Device_Delete(device);
   }

   return 0;
}
