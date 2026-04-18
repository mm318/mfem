//                              MFEM C Example 8
//
// Description: This example demonstrates a primal DPG discretization of the
//              Poisson problem using a 2x2 block normal equation system.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int visualization = 1;
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
                                       "-vis", "--visualization",
                                       "-no-vis", "--no-visualization",
                                       &visualization);
      if (parsed == 1) { continue; }

      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      return 1;
   }

   // 2. Read the mesh from the given mesh file and uniformly refine it.
   {
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      CMFEM_H1FeCollection *x0_fec;
      CMFEM_RtTraceFeCollection *xhat_fec;
      CMFEM_L2FeCollection *test_fec;
      CMFEM_FiniteElementSpace *x0_space;
      CMFEM_FiniteElementSpace *xhat_space;
      CMFEM_FiniteElementSpace *test_space;
      CMFEM_ConstantCoefficient *one;
      CMFEM_LinearForm *F;
      CMFEM_MixedBilinearForm *B0;
      CMFEM_MixedBilinearForm *Bhat;
      CMFEM_BilinearForm *Sinv;
      CMFEM_BilinearForm *S0;
      CMFEM_SparseMatrix *matB0;
      CMFEM_SparseMatrix *matBhat;
      CMFEM_SparseMatrix *matSinv;
      CMFEM_SparseMatrix *matS0;
      CMFEM_SparseMatrix *Shat;
      CMFEM_BlockOperator *Bop;
      CMFEM_RAPOperator *A;
      CMFEM_CGSolver *S0inv;
      CMFEM_CGSolver *Shatinv;
      CMFEM_BlockDiagonalPreconditioner *P;
      CMFEM_GridFunction *x0;
      int dim;
      int mesh_generator;
      int ref_levels;
      int trial_order;
      int trace_order;
      int test_order;
      int s0;
      int s1;
      int s_test;

      _Alignas(max_align_t) CMFEM_ArrayInt offsets =
         CMFEM_ArrayInt_ConstructSize(3);
      _Alignas(max_align_t) CMFEM_ArrayInt offsets_test =
         CMFEM_ArrayInt_ConstructSize(2);
      _Alignas(max_align_t) CMFEM_ArrayInt ess_bdr =
         CMFEM_ArrayInt_Construct();
      _Alignas(max_align_t) CMFEM_Vector x = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector b = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector zero_x0 = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector SinvF = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector Fvec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector LSres = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector negF = CMFEM_Vector_Construct();

      dim = CMFEM_Mesh_Dimension(mesh);
      mesh_generator = CMFEM_Mesh_MeshGenerator(mesh);
      ref_levels = (int)floor(log(10000.0 / CMFEM_Mesh_GetNE(mesh)) /
                              log(2.0) / dim);
      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      // 3. Define the trial, trace, and test DPG spaces.
      trial_order = order;
      trace_order = order - 1;
      test_order = order;
      if (dim == 2 && (order % 2 == 0 || ((mesh_generator & 2) != 0 && order > 1)))
      {
         test_order++;
      }
      if (test_order < trial_order)
      {
         fprintf(stderr,
                 "Warning, test space not enriched enough to handle primal"
                 " trial space\n");
      }
      x0_fec = CMFEM_H1FeCollection_NewOrderDim(trial_order, dim);
      xhat_fec = CMFEM_RtTraceFeCollection_NewOrderDim(trace_order, dim);
      test_fec = CMFEM_L2FeCollection_NewOrderDim(test_order, dim);

      x0_space = CMFEM_FiniteElementSpace_NewMeshH1(mesh, x0_fec);
      xhat_space = CMFEM_FiniteElementSpace_NewMeshRtTrace(mesh, xhat_fec);
      test_space = CMFEM_FiniteElementSpace_NewMeshL2Vdim(mesh, test_fec, 1);

      s0 = CMFEM_FiniteElementSpace_GetVSize(x0_space);
      s1 = CMFEM_FiniteElementSpace_GetVSize(xhat_space);
      s_test = CMFEM_FiniteElementSpace_GetVSize(test_space);

      CMFEM_ArrayInt_Set(&offsets, 0, 0);
      CMFEM_ArrayInt_Set(&offsets, 1, s0);
      CMFEM_ArrayInt_Set(&offsets, 2, s0 + s1);
      CMFEM_ArrayInt_Set(&offsets_test, 0, 0);
      CMFEM_ArrayInt_Set(&offsets_test, 1, s_test);

      printf("\nNumber of Unknowns:\n");
      printf(" Trial space,     X0   : %d (order %d)\n", s0, trial_order);
      printf(" Interface space, Xhat : %d (order %d)\n", s1, trace_order);
      printf(" Test space,      Y    : %d (order %d)\n\n", s_test, test_order);

      CMFEM_Vector_Destroy(&x);
      CMFEM_Vector_Destroy(&b);
      CMFEM_Vector_Destroy(&zero_x0);
      CMFEM_Vector_Destroy(&SinvF);
      CMFEM_Vector_Destroy(&Fvec);
      CMFEM_Vector_Destroy(&LSres);
      CMFEM_Vector_Destroy(&negF);
      x = CMFEM_Vector_ConstructSize(s0 + s1);
      b = CMFEM_Vector_ConstructSize(s0 + s1);
      zero_x0 = CMFEM_Vector_ConstructSize(s0);
      SinvF = CMFEM_Vector_ConstructSize(s_test);
      Fvec = CMFEM_Vector_ConstructSize(s_test);
      LSres = CMFEM_Vector_ConstructSize(s_test);
      negF = CMFEM_Vector_ConstructSize(s_test);
      CMFEM_Vector_Assign(&x, 0.0);
      CMFEM_Vector_Assign(&b, 0.0);
      CMFEM_Vector_Assign(&zero_x0, 0.0);

      // 4. Build the right-hand side linear form on the test space.
      one = CMFEM_ConstantCoefficient_New(1.0);
      F = CMFEM_LinearForm_New(test_space);
      CMFEM_LinearForm_AddDomainIntegratorDliCc(F, one);
      CMFEM_LinearForm_Assemble(F);

      // 5. Assemble the mixed blocks B0 and Bhat, and the test/trial norms.
      ess_bdr = CMFEM_ArrayInt_ConstructSize(CMFEM_Mesh_BoundaryAttributesMax(mesh));
      CMFEM_ArrayInt_Assign(&ess_bdr, 1);

      B0 = CMFEM_MixedBilinearForm_New(x0_space, test_space);
      CMFEM_MixedBilinearForm_AddDomainIntegratorDiCc(B0, one);
      CMFEM_MixedBilinearForm_Assemble(B0);
      CMFEM_MixedBilinearForm_EliminateTrialEssentialBCAiVecLf(
         B0, &ess_bdr, &zero_x0, F);
      CMFEM_MixedBilinearForm_Finalize(B0);

      Bhat = CMFEM_MixedBilinearForm_New(xhat_space, test_space);
      CMFEM_MixedBilinearForm_AddTraceFaceIntegratorTji(Bhat);
      CMFEM_MixedBilinearForm_Assemble(Bhat);
      CMFEM_MixedBilinearForm_Finalize(Bhat);

      Sinv = CMFEM_BilinearForm_New(test_space);
      CMFEM_BilinearForm_AddDomainIntegratorIiDiMiCc(Sinv, one);
      CMFEM_BilinearForm_Assemble(Sinv);
      CMFEM_BilinearForm_Finalize(Sinv);

      S0 = CMFEM_BilinearForm_New(x0_space);
      CMFEM_BilinearForm_AddDomainIntegratorDiCc(S0, one);
      CMFEM_BilinearForm_Assemble(S0);
      CMFEM_BilinearForm_EliminateEssentialBCAi(S0, &ess_bdr);
      CMFEM_BilinearForm_Finalize(S0);

      matB0 = CMFEM_MixedBilinearForm_SpMat(B0);
      matBhat = CMFEM_MixedBilinearForm_SpMat(Bhat);
      matSinv = CMFEM_BilinearForm_SpMat(Sinv);
      matS0 = CMFEM_BilinearForm_SpMat(S0);

      // 6. Form the 1x2 block operator B and the normal-equation operator A.
      Bop = CMFEM_BlockOperator_NewAiAi(&offsets_test, &offsets);
      CMFEM_BlockOperator_SetBlockSm(Bop, 0, 0, matB0);
      CMFEM_BlockOperator_SetBlockSm(Bop, 0, 1, matBhat);
      A = CMFEM_RAPOperator_NewBopSmBop(Bop, matSinv, Bop);

      CMFEM_LinearForm_CopyToVector(F, &Fvec);
      CMFEM_SparseMatrix_Mult(matSinv, &Fvec, &SinvF);
      CMFEM_BlockOperator_MultTranspose(Bop, &SinvF, &b);

      // 7. Build the block-diagonal preconditioner for the normal equation.
      Shat = CMFEM_RAPSmSmSm(matBhat, matSinv, matBhat);

      S0inv = CMFEM_CGSolver_New();
      CMFEM_CGSolver_SetOperatorSm(S0inv, matS0);
      CMFEM_CGSolver_SetPrintLevel(S0inv, -1);
      CMFEM_CGSolver_SetRelTol(S0inv, 1e-3);
      CMFEM_CGSolver_SetMaxIter(S0inv, 200);
      CMFEM_CGSolver_SetIterativeMode(S0inv, 0);

      Shatinv = CMFEM_CGSolver_New();
      CMFEM_CGSolver_SetOperatorSm(Shatinv, Shat);
      CMFEM_CGSolver_SetPrintLevel(Shatinv, -1);
      CMFEM_CGSolver_SetRelTol(Shatinv, 1e-3);
      CMFEM_CGSolver_SetMaxIter(Shatinv, 200);
      CMFEM_CGSolver_SetIterativeMode(Shatinv, 0);

      P = CMFEM_BlockDiagonalPreconditioner_NewAi(&offsets);
      CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockCgs(P, 0, S0inv);
      CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockCgs(P, 1, Shatinv);

      // 8. Solve the normal equation system with PCG.
      CMFEM_PCGRAPOperatorBdp(A, P, &b, &x, 1, 200, 1e-12, 0.0);

      // 9. Check the weighted residual norm for the DPG least-squares system.
      CMFEM_BlockOperator_Mult(Bop, &x, &LSres);
      CMFEM_Vector_Destroy(&negF);
      negF = CMFEM_Vector_Copy(&Fvec);
      CMFEM_Vector_Neg(&negF);
      CMFEM_Vector_Add(&LSres, &negF);
      printf("\n|| B0*x0 + Bhat*xhat - F ||_{S^-1} = %.8g\n",
             sqrt(CMFEM_SparseMatrix_InnerProduct(matSinv, &LSres, &LSres)));

      // 10. Extract the primal solution block into a GridFunction for output.
      x0 = CMFEM_GridFunction_New(x0_space);
      for (i = 0; i < s0; i++)
      {
         CMFEM_GridFunction_Set(x0, i, CMFEM_Vector_Get(&x, i));
      }

      // 11. Save the refined mesh and the primal solution.
      CMFEM_Mesh_Print(mesh, "refined.mesh", 8);
      CMFEM_GridFunction_Save(x0, "sol.gf", 8);

      // 12. Send the primal solution to GLVis.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, x0, "localhost", 19916);
      }

      CMFEM_GridFunction_Delete(x0);
      CMFEM_BlockDiagonalPreconditioner_Delete(P);
      CMFEM_CGSolver_Delete(Shatinv);
      CMFEM_CGSolver_Delete(S0inv);
      CMFEM_SparseMatrix_Delete(Shat);
      CMFEM_RAPOperator_Delete(A);
      CMFEM_BlockOperator_Delete(Bop);
      CMFEM_BilinearForm_Delete(S0);
      CMFEM_BilinearForm_Delete(Sinv);
      CMFEM_MixedBilinearForm_Delete(Bhat);
      CMFEM_MixedBilinearForm_Delete(B0);
      CMFEM_LinearForm_Delete(F);
      CMFEM_ConstantCoefficient_Delete(one);
      CMFEM_FiniteElementSpace_Delete(test_space);
      CMFEM_FiniteElementSpace_Delete(xhat_space);
      CMFEM_FiniteElementSpace_Delete(x0_space);
      CMFEM_L2FeCollection_Delete(test_fec);
      CMFEM_RtTraceFeCollection_Delete(xhat_fec);
      CMFEM_H1FeCollection_Delete(x0_fec);
      CMFEM_ArrayInt_Destroy(&ess_bdr);
      CMFEM_ArrayInt_Destroy(&offsets_test);
      CMFEM_ArrayInt_Destroy(&offsets);
      CMFEM_Vector_Destroy(&Fvec);
      CMFEM_Vector_Destroy(&SinvF);
      CMFEM_Vector_Destroy(&zero_x0);
      CMFEM_Vector_Destroy(&b);
      CMFEM_Vector_Destroy(&x);
      CMFEM_Vector_Destroy(&negF);
      CMFEM_Vector_Destroy(&LSres);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
