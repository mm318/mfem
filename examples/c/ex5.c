//                              MFEM C Example 5
//
// Description: This example solves a mixed Darcy problem using Raviart-Thomas
//              elements for the velocity and discontinuous L2 elements for the
//              pressure.

#include "common.h"
#include "cmfem.h"

#include <math.h>
#include <stdio.h>

static void uFun_ex(const CMFEM_Vector *x, CMFEM_Vector *u, void *context)
{
   double xi = CMFEM_Vector_Get(x, 0);
   double yi = CMFEM_Vector_Get(x, 1);
   double zi = 0.0;
   (void)context;

   if (CMFEM_Vector_Size(x) == 3)
   {
      zi = CMFEM_Vector_Get(x, 2);
   }

   CMFEM_Vector_Set(u, 0, -exp(xi) * sin(yi) * cos(zi));
   CMFEM_Vector_Set(u, 1, -exp(xi) * cos(yi) * cos(zi));

   if (CMFEM_Vector_Size(x) == 3)
   {
      CMFEM_Vector_Set(u, 2, exp(xi) * sin(yi) * sin(zi));
   }
}

static double pFun_ex(const CMFEM_Vector *x, void *context)
{
   double xi = CMFEM_Vector_Get(x, 0);
   double yi = CMFEM_Vector_Get(x, 1);
   double zi = 0.0;
   (void)context;

   if (CMFEM_Vector_Size(x) == 3)
   {
      zi = CMFEM_Vector_Get(x, 2);
   }

   return exp(xi) * sin(yi) * cos(zi);
}

static void fFun(const CMFEM_Vector *x, CMFEM_Vector *f, void *context)
{
   int i;
   (void)x;
   (void)context;
   for (i = 0; i < CMFEM_Vector_Size(f); i++)
   {
      CMFEM_Vector_Set(f, i, 0.0);
   }
}

static double gFun(const CMFEM_Vector *x, void *context)
{
   (void)context;
   if (CMFEM_Vector_Size(x) == 3)
   {
      return -pFun_ex(x, NULL);
   }
   return 0.0;
}

static double f_natural(const CMFEM_Vector *x, void *context)
{
   (void)context;
   return -pFun_ex(x, NULL);
}

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = CMFEM_ExamplesDataPath("star.mesh");
   int order = 1;
   int pa = 0;
   const char *device_config = "cpu";
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
                                       "-pa", "--partial-assembly",
                                       "-no-pa", "--no-partial-assembly",
                                       &pa);
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

   if (pa)
   {
      fprintf(stderr,
              "The C port of ex5 does not support partial assembly yet.\n");
      return 2;
   }

   // 2. Enable the requested device configuration.
   {
      CMFEM_Device *device = CMFEM_Device_New(device_config);
      CMFEM_Device_Print(device);
      CMFEM_Device_Delete(device);
   }

   // 3. Read the mesh from the given mesh file and uniformly refine it.
   {
      CMFEM_Mesh *mesh = CMFEM_Mesh_NewFile(mesh_file, 1, 1);
      CMFEM_RtFeCollection *hdiv_coll;
      CMFEM_L2FeCollection *l2_coll;
      CMFEM_FiniteElementSpace *r_space;
      CMFEM_FiniteElementSpace *w_space;
      CMFEM_ConstantCoefficient *k;
      CMFEM_VectorFunctionCoefficient *fcoeff;
      CMFEM_FunctionCoefficient *fnatcoeff;
      CMFEM_FunctionCoefficient *gcoeff;
      CMFEM_VectorFunctionCoefficient *ucoeff;
      CMFEM_FunctionCoefficient *pcoeff;
      CMFEM_LinearForm *fform;
      CMFEM_LinearForm *gform;
      CMFEM_BilinearForm *m_varf;
      CMFEM_MixedBilinearForm *b_varf;
      CMFEM_SparseMatrix *m_matrix;
      CMFEM_SparseMatrix *b_matrix;
      CMFEM_SparseMatrix *bt_matrix;
      CMFEM_SparseMatrix *minv_bt;
      CMFEM_SparseMatrix *s_matrix;
      CMFEM_BlockOperator *darcy_op;
      CMFEM_BlockDiagonalPreconditioner *darcy_prec;
      CMFEM_DSmoother *inv_m;
      CMFEM_GSSmoother *inv_s;
      CMFEM_GridFunction *u;
      CMFEM_GridFunction *p;
      CMFEM_VisItDataCollection *visit_dc;
      CMFEM_ParaViewDataCollection *paraview_dc;
      int dim;
      int ref_levels;
      int r_size;
      int w_size;
      int quadrature_order;
      double err_u;
      double norm_u;
      double err_p;
      double norm_p;

      _Alignas(max_align_t) CMFEM_ArrayInt block_offsets =
         CMFEM_ArrayInt_ConstructSize(3);
      _Alignas(max_align_t) CMFEM_Vector rhs = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector x = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector fvec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector gvec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector md = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector u_vec = CMFEM_Vector_Construct();
      _Alignas(max_align_t) CMFEM_Vector p_vec = CMFEM_Vector_Construct();

      dim = CMFEM_Mesh_Dimension(mesh);
      ref_levels =
         cmfem_uniform_refinement_levels(10000.0, CMFEM_Mesh_GetNE(mesh), dim);
      for (i = 0; i < ref_levels; i++)
      {
         CMFEM_Mesh_UniformRefinement(mesh);
      }

      // 4. Define finite element spaces for the velocity and pressure.
      hdiv_coll = CMFEM_RtFeCollection_NewOrderDim(order, dim);
      l2_coll = CMFEM_L2FeCollection_NewOrderDim(order, dim);
      r_space = CMFEM_FiniteElementSpace_NewMeshRt(mesh, hdiv_coll);
      w_space = CMFEM_FiniteElementSpace_NewMeshL2Vdim(mesh, l2_coll, 1);

      r_size = CMFEM_FiniteElementSpace_GetVSize(r_space);
      w_size = CMFEM_FiniteElementSpace_GetVSize(w_space);
      CMFEM_ArrayInt_Set(&block_offsets, 0, 0);
      CMFEM_ArrayInt_Set(&block_offsets, 1, r_size);
      CMFEM_ArrayInt_Set(&block_offsets, 2, r_size + w_size);

      printf("***********************************************************\n");
      printf("dim(R) = %d\n", r_size);
      printf("dim(W) = %d\n", w_size);
      printf("dim(R+W) = %d\n", r_size + w_size);
      printf("***********************************************************\n");

      // 5. Define the exact solution and forcing data.
      k = CMFEM_ConstantCoefficient_New(1.0);
      fcoeff = CMFEM_VectorFunctionCoefficient_New(dim, fFun, NULL);
      fnatcoeff = CMFEM_FunctionCoefficient_New(f_natural, NULL);
      gcoeff = CMFEM_FunctionCoefficient_New(gFun, NULL);
      ucoeff = CMFEM_VectorFunctionCoefficient_New(dim, uFun_ex, NULL);
      pcoeff = CMFEM_FunctionCoefficient_New(pFun_ex, NULL);

      // 6. Assemble the right-hand side linear forms.
      fform = CMFEM_LinearForm_New(r_space);
      CMFEM_LinearForm_AddDomainIntegratorVfd(fform, fcoeff);
      CMFEM_LinearForm_AddBoundaryIntegratorVfbfFc(fform, fnatcoeff);
      CMFEM_LinearForm_Assemble(fform);

      gform = CMFEM_LinearForm_New(w_space);
      CMFEM_LinearForm_AddDomainIntegratorDliFc(gform, gcoeff);
      CMFEM_LinearForm_Assemble(gform);

      fvec = CMFEM_Vector_ConstructSize(r_size);
      gvec = CMFEM_Vector_ConstructSize(w_size);
      rhs = CMFEM_Vector_ConstructSize(r_size + w_size);
      x = CMFEM_Vector_ConstructSize(r_size + w_size);
      md = CMFEM_Vector_ConstructSize(r_size);
      u_vec = CMFEM_Vector_ConstructSize(r_size);
      p_vec = CMFEM_Vector_ConstructSize(w_size);
      CMFEM_Vector_Assign(&rhs, 0.0);
      CMFEM_Vector_Assign(&x, 0.0);
      CMFEM_LinearForm_CopyToVector(fform, &fvec);
      CMFEM_LinearForm_CopyToVector(gform, &gvec);
      for (i = 0; i < r_size; i++)
      {
         CMFEM_Vector_Set(&rhs, i, CMFEM_Vector_Get(&fvec, i));
      }
      for (i = 0; i < w_size; i++)
      {
         CMFEM_Vector_Set(&rhs, r_size + i, CMFEM_Vector_Get(&gvec, i));
      }

      // 7. Assemble the mixed Darcy operator.
      m_varf = CMFEM_BilinearForm_New(r_space);
      CMFEM_BilinearForm_AddDomainIntegratorVmi(m_varf, k);
      CMFEM_BilinearForm_Assemble(m_varf);
      CMFEM_BilinearForm_Finalize(m_varf);

      b_varf = CMFEM_MixedBilinearForm_New(r_space, w_space);
      CMFEM_MixedBilinearForm_AddDomainIntegratorVfdiCc(b_varf, k);
      CMFEM_MixedBilinearForm_Assemble(b_varf);
      CMFEM_MixedBilinearForm_Finalize(b_varf);

      m_matrix = CMFEM_BilinearForm_SpMat(m_varf);
      b_matrix = CMFEM_MixedBilinearForm_SpMat(b_varf);
      bt_matrix = CMFEM_TransposeSm(b_matrix);
      CMFEM_SparseMatrix_Scale(b_matrix, -1.0);
      CMFEM_SparseMatrix_Scale(bt_matrix, -1.0);

      darcy_op = CMFEM_BlockOperator_NewAiAi(&block_offsets, &block_offsets);
      CMFEM_BlockOperator_SetBlockSm(darcy_op, 0, 0, m_matrix);
      CMFEM_BlockOperator_SetBlockSm(darcy_op, 0, 1, bt_matrix);
      CMFEM_BlockOperator_SetBlockSm(darcy_op, 1, 0, b_matrix);

      // 8. Build the block-diagonal preconditioner.
      CMFEM_SparseMatrix_GetDiag(m_matrix, &md);
      minv_bt = CMFEM_TransposeSm(b_matrix);
      for (i = 0; i < r_size; i++)
      {
         double diag_i = CMFEM_Vector_Get(&md, i);
         if (diag_i != 0.0)
         {
            CMFEM_SparseMatrix_ScaleRow(minv_bt, i, 1.0 / diag_i);
         }
      }
      s_matrix = CMFEM_MultSmSm(b_matrix, minv_bt);

      inv_m = CMFEM_DSmoother_NewSm(m_matrix);
      inv_s = CMFEM_GSSmoother_NewSm(s_matrix);
      CMFEM_DSmoother_SetIterativeMode(inv_m, 0);
      CMFEM_GSSmoother_SetIterativeMode(inv_s, 0);
      darcy_prec = CMFEM_BlockDiagonalPreconditioner_NewAi(&block_offsets);
      CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockDs(darcy_prec, 0, inv_m);
      CMFEM_BlockDiagonalPreconditioner_SetDiagonalBlockGs(darcy_prec, 1, inv_s);

      // 9. Solve the linear system with MINRES.
      CMFEM_MINRESBopBdp(darcy_op, darcy_prec, &rhs, &x, 1, 1000, 1e-6, 1e-10);

      // 10. Recover the velocity and pressure grid functions and compute errors.
      for (i = 0; i < r_size; i++)
      {
         CMFEM_Vector_Set(&u_vec, i, CMFEM_Vector_Get(&x, i));
      }
      for (i = 0; i < w_size; i++)
      {
         CMFEM_Vector_Set(&p_vec, i, CMFEM_Vector_Get(&x, r_size + i));
      }

      u = CMFEM_GridFunction_New(r_space);
      p = CMFEM_GridFunction_New(w_space);
      CMFEM_GridFunction_SetFromVec(u, &u_vec);
      CMFEM_GridFunction_SetFromVec(p, &p_vec);

      quadrature_order = order * 2 + 1;
      if (quadrature_order < 2)
      {
         quadrature_order = 2;
      }

      err_u = CMFEM_GridFunction_ComputeL2ErrorVfcOrder(u, ucoeff,
                                                        quadrature_order);
      norm_u = CMFEM_ComputeLpNormVfcMeshOrder(2.0, ucoeff, mesh,
                                               quadrature_order);
      err_p = CMFEM_GridFunction_ComputeL2ErrorFcOrder(p, pcoeff,
                                                       quadrature_order);
      norm_p = CMFEM_ComputeLpNormFcMeshOrder(2.0, pcoeff, mesh,
                                              quadrature_order);

      printf("|| u_h - u_ex || / || u_ex || = %.8g\n", err_u / norm_u);
      printf("|| p_h - p_ex || / || p_ex || = %.8g\n", err_p / norm_p);

      // 11. Save the mesh and solution fields.
      CMFEM_Mesh_Print(mesh, "ex5.mesh", 8);
      CMFEM_GridFunction_Save(u, "sol_u.gf", 8);
      CMFEM_GridFunction_Save(p, "sol_p.gf", 8);

      // 12. Save data collections for VisIt and ParaView.
      visit_dc = CMFEM_VisItDataCollection_New("Example5", mesh);
      CMFEM_VisItDataCollection_RegisterFieldGf(visit_dc, "velocity", u);
      CMFEM_VisItDataCollection_RegisterFieldGf(visit_dc, "pressure", p);
      CMFEM_VisItDataCollection_Save(visit_dc);

      paraview_dc = CMFEM_ParaViewDataCollection_New("Example5", mesh);
      CMFEM_ParaViewDataCollection_SetPrefixPath(paraview_dc, "ParaView");
      CMFEM_ParaViewDataCollection_SetLevelsOfDetail(paraview_dc, order);
      CMFEM_ParaViewDataCollection_SetCycle(paraview_dc, 0);
      CMFEM_ParaViewDataCollection_SetDataFormatBinary(paraview_dc);
      CMFEM_ParaViewDataCollection_SetHighOrderOutput(paraview_dc, 1);
      CMFEM_ParaViewDataCollection_SetTime(paraview_dc, 0.0);
      CMFEM_ParaViewDataCollection_RegisterFieldGf(paraview_dc, "velocity", u);
      CMFEM_ParaViewDataCollection_RegisterFieldGf(paraview_dc, "pressure", p);
      CMFEM_ParaViewDataCollection_Save(paraview_dc);

      // 13. Send the solution by socket to a GLVis server.
      if (visualization)
      {
         CMFEM_SendSolutionToGLVis(mesh, u, "localhost", 19916);
         CMFEM_SendSolutionToGLVis(mesh, p, "localhost", 19916);
      }

      CMFEM_ParaViewDataCollection_Delete(paraview_dc);
      CMFEM_VisItDataCollection_Delete(visit_dc);
      CMFEM_GridFunction_Delete(p);
      CMFEM_GridFunction_Delete(u);
      CMFEM_BlockDiagonalPreconditioner_Delete(darcy_prec);
      CMFEM_GSSmoother_Delete(inv_s);
      CMFEM_DSmoother_Delete(inv_m);
      CMFEM_BlockOperator_Delete(darcy_op);
      CMFEM_SparseMatrix_Delete(s_matrix);
      CMFEM_SparseMatrix_Delete(minv_bt);
      CMFEM_SparseMatrix_Delete(bt_matrix);
      CMFEM_MixedBilinearForm_Delete(b_varf);
      CMFEM_BilinearForm_Delete(m_varf);
      CMFEM_LinearForm_Delete(gform);
      CMFEM_LinearForm_Delete(fform);
      CMFEM_FunctionCoefficient_Delete(pcoeff);
      CMFEM_VectorFunctionCoefficient_Delete(ucoeff);
      CMFEM_FunctionCoefficient_Delete(gcoeff);
      CMFEM_FunctionCoefficient_Delete(fnatcoeff);
      CMFEM_VectorFunctionCoefficient_Delete(fcoeff);
      CMFEM_ConstantCoefficient_Delete(k);
      CMFEM_FiniteElementSpace_Delete(w_space);
      CMFEM_FiniteElementSpace_Delete(r_space);
      CMFEM_L2FeCollection_Delete(l2_coll);
      CMFEM_RtFeCollection_Delete(hdiv_coll);
      CMFEM_Vector_Destroy(&p_vec);
      CMFEM_Vector_Destroy(&u_vec);
      CMFEM_Vector_Destroy(&md);
      CMFEM_Vector_Destroy(&gvec);
      CMFEM_Vector_Destroy(&fvec);
      CMFEM_Vector_Destroy(&x);
      CMFEM_Vector_Destroy(&rhs);
      CMFEM_ArrayInt_Destroy(&block_offsets);
      CMFEM_Mesh_Delete(mesh);
   }

   return 0;
}
