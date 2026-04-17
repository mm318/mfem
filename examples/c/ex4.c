//                              MFEM C Example 4
//
// Description: This example code solves a simple 2D/3D H(div) diffusion
//              problem corresponding to the second order definite equation
//              -grad(alpha div F) + beta F = f with boundary condition F dot n
//              equal to a given normal field.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

static double cmfem_ex4_freq = 1.0;
static double cmfem_ex4_kappa = 0.0;
static const double cmfem_ex4_pi = 3.14159265358979323846;

// Exact solution, F, and right-hand side, f.
static void cmfem_ex4_F_exact(const CMFEM_Vector *p, CMFEM_Vector *F,
                              void *context)
{
   double x = CMFEM_Vector_Get(p, 0);
   double y = CMFEM_Vector_Get(p, 1);
   int dim = CMFEM_Vector_Size(p);
   (void)context;

   CMFEM_Vector_Set(F, 0, cos(cmfem_ex4_kappa * x) * sin(cmfem_ex4_kappa * y));
   CMFEM_Vector_Set(F, 1, cos(cmfem_ex4_kappa * y) * sin(cmfem_ex4_kappa * x));
   if (dim == 3)
   {
      CMFEM_Vector_Set(F, 2, 0.0);
   }
}

static void cmfem_ex4_f_exact(const CMFEM_Vector *p, CMFEM_Vector *f,
                              void *context)
{
   double x = CMFEM_Vector_Get(p, 0);
   double y = CMFEM_Vector_Get(p, 1);
   double temp = 1.0 + 2.0 * cmfem_ex4_kappa * cmfem_ex4_kappa;
   int dim = CMFEM_Vector_Size(p);
   (void)context;

   CMFEM_Vector_Set(f, 0, temp * cos(cmfem_ex4_kappa * x) * sin(
                       cmfem_ex4_kappa * y));
   CMFEM_Vector_Set(f, 1, temp * cos(cmfem_ex4_kappa * y) * sin(
                       cmfem_ex4_kappa * x));
   if (dim == 3)
   {
      CMFEM_Vector_Set(f, 2, 0.0);
   }
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int set_bc = 1;
   int static_cond = 0;
   int hybridization = 0;
   int pa = 0;
   int ea = 0;
   const char *device_config = "cpu";
   int visualization = 1;
   int dim;
   int sdim;
   int ref_levels;
   int i;

   for (i = 1; i < argc; i++)
   {
      int parsed = cmfem_parse_string_option(argc, argv, &i, "-m", "--mesh",
                                             &mesh_file);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-o", "--order", &order);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-bc", "--impose-bc",
                                       "-no-bc", "--dont-impose-bc",
                                       &set_bc);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-f", "--frequency",
                                         &cmfem_ex4_freq);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-sc", "--static-condensation",
                                       "-no-sc", "--no-static-condensation",
                                       &static_cond);
      if (parsed == 1) { continue; }

      parsed = cmfem_parse_bool_option(&i, argv,
                                       "-hb", "--hybridization",
                                       "-no-hb", "--no-hybridization",
                                       &hybridization);
      if (parsed == 1) { continue; }

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

      parsed = cmfem_parse_string_option(argc, argv, &i, "-d", "--device",
                                         &device_config);
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

   cmfem_ex4_kappa = cmfem_ex4_freq * cmfem_ex4_pi;

   // 2. Enable hardware devices and programming models based on command line
   //    options.
   CMFEM_Device *device = CMFEM_Device_New(device_config);
   CMFEM_Device_Print(device);

   // 3. Read the mesh from the given mesh file.
   CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
   dim = CMFEM_Mesh_Dimension(mesh);
   sdim = CMFEM_Mesh_SpaceDimension(mesh);
   // 4. Refine the mesh to increase the resolution.
   ref_levels = cmfem_uniform_refinement_levels(25000.0, CMFEM_Mesh_GetNE(mesh),
                                                dim);
   for (i = 0; i < ref_levels; i++)
   {
      CMFEM_Mesh_UniformRefinement(mesh);
   }

   // 5. Define a Raviart-Thomas finite element space on the mesh.
   CMFEM_RtFeCollection *fec = CMFEM_RtFeCollection_NewOrderDim(order - 1, dim);
   CMFEM_FiniteElementSpace *fespace = CMFEM_FiniteElementSpace_NewMeshRt(mesh,
                                                                          fec);
   printf("Number of finite element unknowns: %d\n",
          CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

   // 6. Determine the list of true essential boundary dofs.
   _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
   if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
   {
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
         CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr, set_bc ? 1 : 0);
      CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr,
                                                    &ess_tdof_list);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side
   //    of the discrete system.
   CMFEM_VectorFunctionCoefficient *f =
      CMFEM_VectorFunctionCoefficient_New(sdim, cmfem_ex4_f_exact, NULL);
   CMFEM_LinearForm *b = CMFEM_LinearForm_New(fespace);
   CMFEM_LinearForm_AddDomainIntegratorVfd(b, f);
   CMFEM_LinearForm_Assemble(b);

   // 8. Define the solution vector x as a finite element grid function and
   //    initialize it by projecting the exact solution.
   CMFEM_GridFunction *x = CMFEM_GridFunction_New(fespace);
   CMFEM_VectorFunctionCoefficient *F =
      CMFEM_VectorFunctionCoefficient_New(sdim, cmfem_ex4_F_exact, NULL);
   CMFEM_GridFunction_ProjectCoefficientVfc(x, F);

   // 9. Set up the bilinear form corresponding to the H(div) diffusion
   //    operator.
   CMFEM_ConstantCoefficient *alpha = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_ConstantCoefficient *beta = CMFEM_ConstantCoefficient_New(1.0);
   CMFEM_BilinearForm *a = CMFEM_BilinearForm_New(fespace);
   if (pa)
   {
      CMFEM_BilinearForm_SetAssemblyLevelPartial(a);
   }
   if (ea)
   {
      CMFEM_BilinearForm_SetAssemblyLevelElement(a);
   }
   CMFEM_BilinearForm_AddDomainIntegratorDdi(a, alpha);
   CMFEM_BilinearForm_AddDomainIntegratorVmi(a, beta);

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying static condensation, hybridization, or other required
   //     transformations.
   CMFEM_DgInterfaceFeCollection *hfec = NULL;
   CMFEM_FiniteElementSpace *hfes = NULL;
   if (static_cond)
   {
      CMFEM_BilinearForm_EnableStaticCondensation(a);
   }
   else if (hybridization)
   {
      hfec = CMFEM_DgInterfaceFeCollection_NewOrderDim(order - 1, dim);
      hfes = CMFEM_FiniteElementSpace_NewMeshDgInterface(mesh, hfec);
      CMFEM_BilinearForm_EnableHybridization(a, hfes, &ess_tdof_list);
   }
   CMFEM_BilinearForm_Assemble(a);

   // 11. Form and solve the linear system A X = B.
   _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
   _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
   _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
   CMFEM_BilinearForm_FormLinearSystemOp(a, &ess_tdof_list, x, b, &A, &X,
                                         &B);

   printf("Size of linear system: %d\n", CMFEM_OperatorPtr_Height(&A));

   if (!pa && (!ea || hybridization))
   {
      CMFEM_GSSmoother *M = CMFEM_GSSmoother_NewOp(&A);
      CMFEM_PCGOpGs(&A, M, &B, &X, 1, 10000, 1e-20, 0.0);
      CMFEM_GSSmoother_Delete(M);
   }
   else if (CMFEM_UsesTensorBasis(fespace))
   {
      CMFEM_OperatorJacobiSmoother *M =
         CMFEM_OperatorJacobiSmoother_NewBf(a, &ess_tdof_list);
      CMFEM_PCGOpOjs(&A, M, &B, &X, 1, 10000, 1e-20, 0.0);
      CMFEM_OperatorJacobiSmoother_Delete(M);
   }
   else
   {
      CMFEM_CGOp(&A, &B, &X, 1, 10000, 1e-20, 0.0);
   }

   // 12. Recover the solution as a finite element grid function.
   CMFEM_BilinearForm_RecoverFEMSolution(a, &X, b, x);
   // 13. Compute and print the L^2 norm of the error.
   printf("\n|| F_h - F ||_{L^2} = %.16g\n\n",
          CMFEM_GridFunction_ComputeL2ErrorVfc(x, F));

   // 14. Save the refined mesh and the solution.
   CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
   CMFEM_GridFunction_Save(x, "sol.gf", 8);

   // 15. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      CMFEM_SendSolutionToGLVis(mesh, x, "localhost", 19916);
   }

   // 16. Free the used memory.
   CMFEM_OperatorPtr_Destroy(&A);
   CMFEM_Vector_Destroy(&B);
   CMFEM_Vector_Destroy(&X);
   if (hfes)
   {
      CMFEM_FiniteElementSpace_Delete(hfes);
   }
   if (hfec)
   {
      CMFEM_DgInterfaceFeCollection_Delete(hfec);
   }
   CMFEM_BilinearForm_Delete(a);
   CMFEM_ConstantCoefficient_Delete(alpha);
   CMFEM_ConstantCoefficient_Delete(beta);
   CMFEM_GridFunction_Delete(x);
   CMFEM_VectorFunctionCoefficient_Delete(F);
   CMFEM_LinearForm_Delete(b);
   CMFEM_VectorFunctionCoefficient_Delete(f);
   CMFEM_ArrayInt_Destroy(&ess_tdof_list);
   CMFEM_FiniteElementSpace_Delete(fespace);
   CMFEM_RtFeCollection_Delete(fec);
   CMFEM_Mesh_Delete(mesh);
   CMFEM_Device_Delete(device);
   return 0;
}
