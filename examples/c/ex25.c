//                              MFEM C Example 25
//
// Description: This example solves a complex-valued Maxwell problem with a
//              perfectly matched layer using Nedelec finite elements.

#include "common.h"
#include "cmfem.h"

#include <stdio.h>

enum
{
   ex25_beam = 0,
   ex25_disc = 1,
   ex25_lshape = 2,
   ex25_fichera = 3,
   ex25_load_src = 4
};

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = NULL;
   int order = 1;
   int ref_levels = 3;
   int prob = ex25_load_src;
   double freq = 5.0;
   double mu = 1.0;
   double epsilon = 1.0;
   double omega;
   int herm_conv = 1;
   int visualization = 1;
   int pa = 0;
   int exact_known = 0;
   const char *device_config = "cpu";
   enum CMFEM_ComplexConvention convention;
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

      parsed = cmfem_parse_int_option(argc, argv, &i, "-prob", "--problem",
                                      &prob);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_int_option(argc, argv, &i, "-ref",
                                      "--refinements", &ref_levels);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-mu",
                                         "--permeability", &mu);
      if (parsed == 1) { continue; }
      if (parsed == 0) { return 1; }

      parsed = cmfem_parse_double_option(argc, argv, &i, "-eps",
                                         "--permittivity", &epsilon);
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

   if (prob > ex25_load_src) { prob = ex25_load_src; }

   if (!mesh_file)
   {
      exact_known = 1;
      switch (prob)
      {
         case ex25_beam:
            mesh_file = CMFEM_ExamplesDataPath("beam-hex.mesh");
            break;
         case ex25_disc:
            mesh_file = CMFEM_ExamplesDataPath("square-disc.mesh");
            break;
         case ex25_lshape:
            mesh_file = CMFEM_ExamplesDataPath("l-shape.mesh");
            break;
         case ex25_fichera:
            mesh_file = CMFEM_ExamplesDataPath("fichera.mesh");
            break;
         case ex25_load_src:
         default:
            exact_known = 0;
            mesh_file = CMFEM_ExamplesDataPath("inline-quad.mesh");
            break;
      }
   }

   omega = 2.0 * 3.14159265358979323846264338327950288 * freq;
   convention = herm_conv ? CMFEM_ComplexConventionHermitian :
                CMFEM_ComplexConventionBlockSymmetric;

   // 2. Configure the device, mesh, and PML helper.
   {
      CMFEM_Device *device = CMFEM_Device_New(device_config);
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      void *pml = NULL;
      const int dim = CMFEM_Mesh_Dimension(mesh);
      CMFEM_NdFeCollection *fec = NULL;
      CMFEM_FiniteElementSpace *fespace = NULL;
      _Alignas(max_align_t) CMFEM_ArrayInt ess_tdof_list = CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr = CMFEM_ArrayInt_Construct();
      CMFEM_ComplexLinearForm *b = NULL;
      CMFEM_ComplexGridFunction *x = NULL;
      CMFEM_SesquilinearForm *a = NULL;
      CMFEM_BilinearForm *prec = NULL;
      _Alignas(max_align_t) CMFEM_OperatorPtr A = CMFEM_OperatorPtr_Construct();
      _Alignas(max_align_t) CMFEM_Vector B = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector X = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_ArrayInt offsets = CMFEM_ArrayInt_Construct();
      CMFEM_BlockDiagonalPreconditioner *bdp = NULL;
      CMFEM_OperatorJacobiSmoother *pc_ojs = NULL;
      CMFEM_GSSmoother *pc_gs = NULL;
      CMFEM_ScaledOperator *pc_i = NULL;
      _Alignas(max_align_t) CMFEM_SparseMatrix pc_matrix = {0};
      int have_pc_matrix = 0;
      double rel_re = 0.0;
      double rel_im = 0.0;
      double total = 0.0;
      CMFEM_GridFunction *x_real = NULL;
      CMFEM_GridFunction *x_imag = NULL;

      CMFEM_Device_Print(device);
      pml = CMFEM_Ex25_NewPml(mesh, prob);

      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      CMFEM_Ex25_SetAttributesPml(pml, mesh);

      fec = CMFEM_NdFeCollection_NewOrderDim(order, dim);
      fespace = CMFEM_FiniteElementSpace_NewMeshNd(mesh, fec);
      printf("Number of finite element unknowns: %d\n",
             CMFEM_FiniteElementSpace_GetTrueVSize(fespace));

      if (CMFEM_Mesh_BoundaryAttributesSize(mesh) > 0)
      {
         CMFEM_ArrayInt_Destroy(&ess_bdr);
         ess_bdr = CMFEM_ArrayInt_ConstructSize(
                      CMFEM_Mesh_BoundaryAttributesMax(mesh));
         CMFEM_Ex25_BuildEssentialBoundaryPml(pml, mesh, &ess_bdr);
      }
      CMFEM_FiniteElementSpace_GetEssentialTrueDofs(fespace, &ess_bdr,
                                                    &ess_tdof_list);

      b = CMFEM_ComplexLinearForm_NewFesConv(fespace, convention);
      CMFEM_Ex25_AssembleRhsPml(pml, mu, epsilon, omega, b);

      x = CMFEM_ComplexGridFunction_New(fespace);
      CMFEM_ComplexGridFunction_Assign(x, 0.0, 0.0);
      CMFEM_Ex25_ProjectBoundaryPml(pml, mu, epsilon, omega, x, &ess_bdr);

      a = CMFEM_SesquilinearForm_NewFesConv(fespace, convention);
      if (pa) { CMFEM_SesquilinearForm_SetAssemblyLevelPartial(a); }
      CMFEM_Ex25_AddOperatorIntegratorsPml(pml, mu, epsilon, omega, a);
      CMFEM_SesquilinearForm_AssembleSkipZeros(a, 0);

      CMFEM_SesquilinearForm_FormLinearSystemOp(a, &ess_tdof_list, x, b,
                                                &A, &X, &B);

      prec = CMFEM_BilinearForm_New(fespace);
      if (pa) { CMFEM_BilinearForm_SetAssemblyLevelPartial(prec); }
      CMFEM_Ex25_AddPreconditionerIntegratorsPml(pml, mu, epsilon, omega, prec);
      CMFEM_BilinearForm_Assemble(prec);

      CMFEM_ArrayInt_Destroy(&offsets);
      offsets = CMFEM_ArrayInt_ConstructSize(3);
      CMFEM_ArrayInt_Set(&offsets, 0, 0);
      CMFEM_ArrayInt_Set(&offsets, 1, CMFEM_FiniteElementSpace_GetTrueVSize(fespace));
      CMFEM_ArrayInt_Set(&offsets, 2,
                         2 * CMFEM_FiniteElementSpace_GetTrueVSize(fespace));
      bdp = CMFEM_BlockDiagonalPreconditioner_NewAi(&offsets);

      if (pa)
      {
         pc_ojs = CMFEM_OperatorJacobiSmoother_NewBf(prec, &ess_tdof_list);
         pc_i = CMFEM_ScaledOperator_NewOjs(
                   pc_ojs,
                   convention == CMFEM_ComplexConventionHermitian ? -1.0 : 1.0);
         CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockOjs(bdp, 0, pc_ojs);
         CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockSop(bdp, 1, pc_i);
      }
      else
      {
         CMFEM_BilinearForm_SetDiagonalPolicyOne(prec);
         pc_matrix = CMFEM_SparseMatrix_Construct();
         have_pc_matrix = 1;
         CMFEM_BilinearForm_FormSystemMatrixSm(prec, &ess_tdof_list, &pc_matrix);
         pc_gs = CMFEM_GSSmoother_NewSm(&pc_matrix);
         pc_i = CMFEM_ScaledOperator_NewGs(
                   pc_gs,
                   convention == CMFEM_ComplexConventionHermitian ? -1.0 : 1.0);
         CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockGs(bdp, 0, pc_gs);
         CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockSop(bdp, 1, pc_i);
      }

      CMFEM_GMRESOpBdp(&A, bdp, &B, &X, 1, pa ? 5000 : 2000, 200, 1.0e-5, 0.0);
      CMFEM_SesquilinearForm_RecoverFEMSolution(a, &X, b, x);

      if (exact_known)
      {
         CMFEM_Ex25_ComputeErrorsPml(pml, mu, epsilon, omega, x, order,
                                     &rel_re, &rel_im, &total);
         printf("\n Relative Error (Re part): || E_h - E || / ||E|| = %.8g\n",
                rel_re);
         printf(" Relative Error (Im part): || E_h - E || / ||E|| = %.8g\n",
                rel_im);
         printf(" Total Error: %.8g\n\n", total);
      }

      CMFEM_Mesh_Print(mesh, "ex25.mesh", 8);
      CMFEM_ComplexGridFunction_SaveReal(x, "ex25-sol_r.gf", 8);
      CMFEM_ComplexGridFunction_SaveImag(x, "ex25-sol_i.gf", 8);

      if (visualization)
      {
         x_real = CMFEM_GridFunction_New(fespace);
         x_imag = CMFEM_GridFunction_New(fespace);
         CMFEM_ComplexGridFunction_CopyRealToGf(x, x_real);
         CMFEM_ComplexGridFunction_CopyImagToGf(x, x_imag);
         CMFEM_SendSolutionToGLVis(mesh, x_real, "localhost", 19916);
         CMFEM_SendSolutionToGLVis(mesh, x_imag, "localhost", 19916);
      }

      CMFEM_GridFunction_Delete(x_imag);
      CMFEM_GridFunction_Delete(x_real);
      if (have_pc_matrix) { CMFEM_SparseMatrix_Destroy(&pc_matrix); }
      CMFEM_ScaledOperator_Delete(pc_i);
      CMFEM_GSSmoother_Delete(pc_gs);
      CMFEM_OperatorJacobiSmoother_Delete(pc_ojs);
      CMFEM_BlockDiagonalPreconditioner_Delete(bdp);
      CMFEM_ArrayInt_Destroy(&offsets);
      CMFEM_OperatorPtr_Destroy(&A);
      CMFEM_Vector_Destroy(&B);
      CMFEM_Vector_Destroy(&X);
      CMFEM_BilinearForm_Delete(prec);
      CMFEM_SesquilinearForm_Delete(a);
      CMFEM_ComplexGridFunction_Delete(x);
      CMFEM_ComplexLinearForm_Delete(b);
      CMFEM_ArrayInt_Destroy(&ess_tdof_list);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_FiniteElementSpace_Delete(fespace);
      CMFEM_NdFeCollection_Delete(fec);
      CMFEM_Ex25_DeletePml(pml);
      CMFEM_Mesh_Delete(mesh);
      CMFEM_Device_Delete(device);
   }

   return 0;
}
