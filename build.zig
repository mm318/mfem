const std = @import("std");

// Supported floating-point configurations for the generated config header.
const Precision = enum {
    double,
    single,
};

// Common compiler flags for MFEM C++ translation units.
const cxx_flags = &.{
    "-std=c++17",
};

// Common compiler flags for C translation units.
const c_flags = &.{
    "-std=c17",
};

// Extra flags for C++ executables that need runtime path helpers pre-included.
const executable_cxx_flags = &.{
    "-std=c++17",
    "-include",
    "config/runtime_paths.hpp",
};

// Runtime path resolution is shared by the library and installed executables.
const runtime_c_sources: []const []const u8 = &.{
    "config/runtime_paths.c",
};

// Core MFEM implementation sources grouped by major subsystem.
const general_sources: []const []const u8 = &.{
    "core/general/array.cpp",
    "core/general/binaryio.cpp",
    "core/general/cuda.cpp",
    "core/general/device.cpp",
    "core/general/error.cpp",
    "core/general/gecko.cpp",
    "core/general/globals.cpp",
    "core/general/hash.cpp",
    "core/general/hash_util.cpp",
    "core/general/isockstream.cpp",
    "core/general/mem_manager.cpp",
    "core/general/occa.cpp",
    "core/general/optparser.cpp",
    "core/general/osockstream.cpp",
    "core/general/sets.cpp",
    "core/general/socketstream.cpp",
    "core/general/stable3d.cpp",
    "core/general/table.cpp",
    "core/general/tic_toc.cpp",
    "core/general/tinyxml2.cpp",
    "core/general/version.cpp",
    "core/general/hip.cpp",
};

const linalg_sources: []const []const u8 = &.{
    "core/linalg/auxiliary.cpp",
    "core/linalg/batched/batched.cpp",
    "core/linalg/batched/gpu_blas.cpp",
    "core/linalg/batched/magma.cpp",
    "core/linalg/batched/native.cpp",
    "core/linalg/batched/solver.cpp",
    "core/linalg/blockmatrix.cpp",
    "core/linalg/blockoperator.cpp",
    "core/linalg/blockvector.cpp",
    "core/linalg/complex_densemat.cpp",
    "core/linalg/complex_operator.cpp",
    "core/linalg/constraints.cpp",
    "core/linalg/densemat.cpp",
    "core/linalg/filteredsolver.cpp",
    "core/linalg/handle.cpp",
    "core/linalg/matrix.cpp",
    "core/linalg/mma.cpp",
    "core/linalg/ode.cpp",
    "core/linalg/operator.cpp",
    "core/linalg/ordering.cpp",
    "core/linalg/particlevector.cpp",
    "core/linalg/solvers.cpp",
    "core/linalg/sundials.cpp",
    "core/linalg/sparsemat.cpp",
    "core/linalg/sparsesmoothers.cpp",
    "core/linalg/symmat.cpp",
    "core/linalg/vector.cpp",
};

const mesh_sources: []const []const u8 = &.{
    "core/mesh/attribute_sets.cpp",
    "core/mesh/element.cpp",
    "core/mesh/exodus_writer.cpp",
    "core/mesh/face_nbr_geom.cpp",
    "core/mesh/gmsh.cpp",
    "core/mesh/hexahedron.cpp",
    "core/mesh/mesh.cpp",
    "core/mesh/mesh_operators.cpp",
    "core/mesh/mesh_readers.cpp",
    "core/mesh/ncmesh.cpp",
    "core/mesh/ncnurbs.cpp",
    "core/mesh/nurbs.cpp",
    "core/mesh/point.cpp",
    "core/mesh/pyramid.cpp",
    "core/mesh/quadrilateral.cpp",
    "core/mesh/segment.cpp",
    "core/mesh/spacing.cpp",
    "core/mesh/tetrahedron.cpp",
    "core/mesh/triangle.cpp",
    "core/mesh/vertex.cpp",
    "core/mesh/vtk.cpp",
    "core/mesh/vtkhdf.cpp",
    "core/mesh/wedge.cpp",
    "core/mesh/submesh/submesh.cpp",
    "core/mesh/submesh/ncsubmesh.cpp",
    "core/mesh/submesh/submesh_utils.cpp",
    "core/mesh/submesh/transfermap.cpp",
};

const fem_sources: []const []const u8 = &.{
    "core/fem/bilinearform.cpp",
    "core/fem/bilinearform_ext.cpp",
    "core/fem/bilininteg.cpp",
    "core/fem/integ/bilininteg_br2.cpp",
    "core/fem/integ/bilininteg_convection_mf.cpp",
    "core/fem/integ/bilininteg_convection_pa.cpp",
    "core/fem/integ/bilininteg_convection_ea.cpp",
    "core/fem/integ/bilininteg_curlcurl_pa.cpp",
    "core/fem/integ/bilininteg_dgdiffusion_pa.cpp",
    "core/fem/integ/bilininteg_dgtrace_pa.cpp",
    "core/fem/integ/bilininteg_dgtrace_ea.cpp",
    "core/fem/integ/bilininteg_diffusion_mf.cpp",
    "core/fem/integ/bilininteg_diffusion_pa.cpp",
    "core/fem/integ/bilininteg_diffusion_ea.cpp",
    "core/fem/integ/bilininteg_diffusion_patch.cpp",
    "core/fem/integ/bilininteg_divdiv_pa.cpp",
    "core/fem/integ/bilininteg_elasticity_ea.cpp",
    "core/fem/integ/bilininteg_elasticity_pa.cpp",
    "core/fem/integ/bilininteg_gradient_pa.cpp",
    "core/fem/integ/bilininteg_interp_pa.cpp",
    "core/fem/integ/bilininteg_mass_mf.cpp",
    "core/fem/integ/bilininteg_mass_pa.cpp",
    "core/fem/integ/bilininteg_mass_ea.cpp",
    "core/fem/integ/bilininteg_mixedcurl_pa.cpp",
    "core/fem/integ/bilininteg_mixedvecgrad_pa.cpp",
    "core/fem/integ/bilininteg_trace_jump_ea.cpp",
    "core/fem/integ/bilininteg_transpose_ea.cpp",
    "core/fem/integ/bilininteg_vecdiffusion_mf.cpp",
    "core/fem/integ/bilininteg_vecdiffusion_pa.cpp",
    "core/fem/integ/bilininteg_vecdiv_pa.cpp",
    "core/fem/integ/bilininteg_vecmass_mf.cpp",
    "core/fem/integ/bilininteg_vecmass_pa.cpp",
    "core/fem/integ/bilininteg_vectorfediv_pa.cpp",
    "core/fem/integ/bilininteg_vectorfemass_pa.cpp",
    "core/fem/integ/bilininteg_diffusion_kernels.cpp",
    "core/fem/integ/bilininteg_elasticity_kernels.cpp",
    "core/fem/integ/bilininteg_hcurl_kernels.cpp",
    "core/fem/integ/bilininteg_hdiv_ea.cpp",
    "core/fem/integ/bilininteg_hdiv_kernels.cpp",
    "core/fem/integ/bilininteg_hcurlhdiv_kernels.cpp",
    "core/fem/integ/bilininteg_mass_kernels.cpp",
    "core/fem/integ/lininteg_boundary.cpp",
    "core/fem/integ/lininteg_boundary_flux.cpp",
    "core/fem/integ/lininteg_domain.cpp",
    "core/fem/integ/lininteg_domain_grad.cpp",
    "core/fem/integ/lininteg_domain_vectorfe.cpp",
    "core/fem/integ/nonlininteg_vecconvection_pa.cpp",
    "core/fem/integ/nonlininteg_vecconvection_mf.cpp",
    "core/fem/coefficient.cpp",
    "core/fem/complex_fem.cpp",
    "core/fem/convergence.cpp",
    "core/fem/datacollection.cpp",
    "core/fem/dgmassinv.cpp",
    "core/fem/doftrans.cpp",
    "core/fem/dfem/doperator.cpp",
    "core/fem/eltrans.cpp",
    "core/fem/batchitrans.cpp",
    "core/fem/estimators.cpp",
    "core/fem/fe.cpp",
    "core/fem/fe/face_map_utils.cpp",
    "core/fem/fe/fe_base.cpp",
    "core/fem/fe/fe_fixed_order.cpp",
    "core/fem/fe/fe_h1.cpp",
    "core/fem/fe/fe_l2.cpp",
    "core/fem/fe/fe_nd.cpp",
    "core/fem/fe/fe_nurbs.cpp",
    "core/fem/fe/fe_pos.cpp",
    "core/fem/fe/fe_pyramid.cpp",
    "core/fem/fe/fe_rt.cpp",
    "core/fem/fe/fe_ser.cpp",
    "core/fem/fe_coll.cpp",
    "core/fem/fespace.cpp",
    "core/fem/derefmat_op.cpp",
    "core/fem/pderefmat_op.cpp",
    "core/fem/geom.cpp",
    "core/fem/gridfunc.cpp",
    "core/fem/hybridization.cpp",
    "core/fem/hybridization_ext.cpp",
    "core/fem/intrules.cpp",
    "core/fem/intrules_cut.cpp",
    "core/fem/ceed/interface/basis.cpp",
    "core/fem/ceed/interface/restriction.cpp",
    "core/fem/ceed/interface/operator.cpp",
    "core/fem/ceed/interface/util.cpp",
    "core/fem/ceed/integrators/convection/convection.cpp",
    "core/fem/ceed/integrators/diffusion/diffusion.cpp",
    "core/fem/ceed/integrators/nlconvection/nlconvection.cpp",
    "core/fem/ceed/integrators/mass/mass.cpp",
    "core/fem/ceed/solvers/algebraic.cpp",
    "core/fem/ceed/solvers/full-assembly.cpp",
    "core/fem/ceed/solvers/solvers-atpmg.cpp",
    "core/fem/kdtree.cpp",
    "core/fem/linearform.cpp",
    "core/fem/linearform_ext.cpp",
    "core/fem/lininteg.cpp",
    "core/fem/lor/lor.cpp",
    "core/fem/lor/lor_ads.cpp",
    "core/fem/lor/lor_ams.cpp",
    "core/fem/lor/lor_batched.cpp",
    "core/fem/multigrid.cpp",
    "core/fem/nonlinearform.cpp",
    "core/fem/nonlinearform_ext.cpp",
    "core/fem/nonlininteg.cpp",
    "core/fem/fespacehierarchy.cpp",
    "core/fem/qfunction.cpp",
    "core/fem/qinterp/det.cpp",
    "core/fem/qinterp/eval_by_nodes.cpp",
    "core/fem/qinterp/eval_by_vdim.cpp",
    "core/fem/qinterp/eval_hdiv.cpp",
    "core/fem/qinterp/grad_by_nodes.cpp",
    "core/fem/qinterp/grad_by_vdim.cpp",
    "core/fem/qspace.cpp",
    "core/fem/quadinterpolator.cpp",
    "core/fem/quadinterpolator_face.cpp",
    "core/fem/restriction.cpp",
    "core/fem/normal_deriv_restriction.cpp",
    "core/fem/staticcond.cpp",
    "core/fem/tmop.cpp",
    "core/fem/tmop/pa.cpp",
    "core/fem/tmop/assemble/diag2_limit.cpp",
    "core/fem/tmop/assemble/diag2.cpp",
    "core/fem/tmop/assemble/grad2_limit.cpp",
    "core/fem/tmop/assemble/grad2.cpp",
    "core/fem/tmop/assemble/diag3_limit.cpp",
    "core/fem/tmop/assemble/diag3.cpp",
    "core/fem/tmop/assemble/grad3_limit.cpp",
    "core/fem/tmop/assemble/grad3.cpp",
    "core/fem/tmop/metrics/001.cpp",
    "core/fem/tmop/metrics/002.cpp",
    "core/fem/tmop/metrics/007.cpp",
    "core/fem/tmop/metrics/056.cpp",
    "core/fem/tmop/metrics/077.cpp",
    "core/fem/tmop/metrics/080.cpp",
    "core/fem/tmop/metrics/094.cpp",
    "core/fem/tmop/metrics/302.cpp",
    "core/fem/tmop/metrics/303.cpp",
    "core/fem/tmop/metrics/315.cpp",
    "core/fem/tmop/metrics/318.cpp",
    "core/fem/tmop/metrics/321.cpp",
    "core/fem/tmop/metrics/332.cpp",
    "core/fem/tmop/metrics/338.cpp",
    "core/fem/tmop/mult/grad2_limit.cpp",
    "core/fem/tmop/mult/grad2.cpp",
    "core/fem/tmop/mult/mult2_limit.cpp",
    "core/fem/tmop/mult/mult2.cpp",
    "core/fem/tmop/mult/grad3_limit.cpp",
    "core/fem/tmop/mult/grad3.cpp",
    "core/fem/tmop/mult/mult3_limit.cpp",
    "core/fem/tmop/mult/mult3.cpp",
    "core/fem/tmop/tools/det2_jpr.cpp",
    "core/fem/tmop/tools/det3_jpr.cpp",
    "core/fem/tmop/tools/discrete.cpp",
    "core/fem/tmop/tools/energy2_limit.cpp",
    "core/fem/tmop/tools/energy2.cpp",
    "core/fem/tmop/tools/energy3_limit.cpp",
    "core/fem/tmop/tools/energy3.cpp",
    "core/fem/tmop/tools/target2.cpp",
    "core/fem/tmop/tools/target3.cpp",
    "core/fem/tmop_tools.cpp",
    "core/fem/tmop_amr.cpp",
    "core/fem/gslib.cpp",
    "core/fem/gslib/findpts_local_2.cpp",
    "core/fem/gslib/findpts_local_3.cpp",
    "core/fem/gslib/interpolate_local_2.cpp",
    "core/fem/gslib/interpolate_local_3.cpp",
    "core/fem/transfer.cpp",
    "core/fem/hyperbolic.cpp",
    "core/fem/integrator.cpp",
    "core/fem/bounds.cpp",
    "core/fem/particleset.cpp",
};

const example_main_sources: []const []const u8 = &.{
    "examples/cpp/ex0.cpp",
    // "examples/cpp/ex0p.cpp", // MPI
    "examples/cpp/ex1.cpp",
    // "examples/cpp/ex1p.cpp", // MPI
    "examples/cpp/ex2.cpp",
    // "examples/cpp/ex2p.cpp", // MPI
    "examples/cpp/ex3.cpp",
    // "examples/cpp/ex3p.cpp", // MPI
    "examples/cpp/ex4.cpp",
    // "examples/cpp/ex4p.cpp", // MPI
    "examples/cpp/ex5.cpp",
    // "examples/cpp/ex5p.cpp", // MPI
    "examples/cpp/ex6.cpp",
    // "examples/cpp/ex6p.cpp", // MPI
    "examples/cpp/ex7.cpp",
    // "examples/cpp/ex7p.cpp", // MPI
    "examples/cpp/ex8.cpp",
    // "examples/cpp/ex8p.cpp", // MPI
    "examples/cpp/ex9.cpp",
    // "examples/cpp/ex9p.cpp", // MPI
    "examples/cpp/ex10.cpp",
    // "examples/cpp/ex10p.cpp", // MPI
    // "examples/cpp/ex11p.cpp", // MPI
    // "examples/cpp/ex12p.cpp", // MPI
    // "examples/cpp/ex13p.cpp", // MPI
    "examples/cpp/ex14.cpp",
    // "examples/cpp/ex14p.cpp", // MPI
    "examples/cpp/ex15.cpp",
    // "examples/cpp/ex15p.cpp", // MPI
    "examples/cpp/ex16.cpp",
    // "examples/cpp/ex16p.cpp", // MPI
    "examples/cpp/ex17.cpp",
    // "examples/cpp/ex17p.cpp", // MPI
    "examples/cpp/ex18.cpp",
    // "examples/cpp/ex18p.cpp", // MPI
    "examples/cpp/ex19.cpp",
    // "examples/cpp/ex19p.cpp", // MPI
    "examples/cpp/ex20.cpp",
    // "examples/cpp/ex20p.cpp", // MPI
    "examples/cpp/ex21.cpp",
    // "examples/cpp/ex21p.cpp", // MPI
    "examples/cpp/ex22.cpp",
    // "examples/cpp/ex22p.cpp", // MPI
    "examples/cpp/ex23.cpp",
    "examples/cpp/ex24.cpp",
    // "examples/cpp/ex24p.cpp", // MPI
    "examples/cpp/ex25.cpp",
    // "examples/cpp/ex25p.cpp", // MPI
    "examples/cpp/ex26.cpp",
    // "examples/cpp/ex26p.cpp", // MPI
    "examples/cpp/ex27.cpp",
    // "examples/cpp/ex27p.cpp", // MPI
    "examples/cpp/ex28.cpp",
    // "examples/cpp/ex28p.cpp", // MPI
    "examples/cpp/ex29.cpp",
    // "examples/cpp/ex29p.cpp", // MPI
    "examples/cpp/ex30.cpp",
    // "examples/cpp/ex30p.cpp", // MPI
    "examples/cpp/ex31.cpp",
    // "examples/cpp/ex31p.cpp", // MPI
    // "examples/cpp/ex32p.cpp", // MPI
    "examples/cpp/ex33.cpp",
    // "examples/cpp/ex33p.cpp", // MPI
    "examples/cpp/ex34.cpp",
    // "examples/cpp/ex34p.cpp", // MPI
    // "examples/cpp/ex35p.cpp", // MPI
    "examples/cpp/ex36.cpp",
    // "examples/cpp/ex36p.cpp", // MPI
    "examples/cpp/ex37.cpp",
    // "examples/cpp/ex37p.cpp", // MPI
    "examples/cpp/ex38.cpp",
    "examples/cpp/ex39.cpp",
    // "examples/cpp/ex39p.cpp", // MPI
    "examples/cpp/ex40.cpp",
    // "examples/cpp/ex40p.cpp", // MPI
    "examples/cpp/ex41.cpp",
    // "examples/cpp/ex41p.cpp", // MPI
    // "examples/cpp/amgx/ex1.cpp", // AmgX
    // "examples/cpp/amgx/ex1p.cpp", // AmgX + MPI
    // "examples/cpp/caliper/ex1.cpp", // Caliper
    // "examples/cpp/caliper/ex1p.cpp", // Caliper + MPI
    // "examples/cpp/ginkgo/ex1.cpp", // Ginkgo
    // "examples/cpp/ginkgo/ex1p.cpp", // Ginkgo + MPI
    // "examples/cpp/hiop/ex9.cpp", // HiOp
    // "examples/cpp/hiop/ex9p.cpp", // HiOp + MPI
    // "examples/cpp/moonolith/ex1.cpp", // Moonolith
    // "examples/cpp/moonolith/ex1p.cpp", // Moonolith + MPI
    // "examples/cpp/moonolith/ex2p.cpp", // Moonolith + MPI
    // "examples/cpp/petsc/ex1p.cpp", // PETSc + MPI
    // "examples/cpp/petsc/ex2p.cpp", // PETSc + MPI
    // "examples/cpp/petsc/ex3p.cpp", // PETSc + MPI
    // "examples/cpp/petsc/ex4p.cpp", // PETSc + MPI
    // "examples/cpp/petsc/ex5p.cpp", // PETSc + MPI
    // "examples/cpp/petsc/ex6p.cpp", // PETSc + MPI
    // "examples/cpp/petsc/ex9p.cpp", // PETSc + MPI
    // "examples/cpp/petsc/ex10p.cpp", // PETSc + MPI
    // "examples/cpp/petsc/ex11p.cpp", // PETSc + MPI
    // "examples/cpp/pumi/ex1.cpp", // PUMI
    // "examples/cpp/pumi/ex2.cpp", // PUMI
    // "examples/cpp/pumi/ex1p.cpp", // PUMI + MPI
    // "examples/cpp/pumi/ex6p.cpp", // PUMI + MPI
    "examples/cpp/sundials/ex9.cpp", // SUNDIALS
    "examples/cpp/sundials/ex10.cpp", // SUNDIALS
    "examples/cpp/sundials/ex16.cpp", // SUNDIALS
    // "examples/cpp/sundials/ex9p.cpp", // SUNDIALS + MPI
    // "examples/cpp/sundials/ex10p.cpp", // SUNDIALS + MPI
    // "examples/cpp/sundials/ex16p.cpp", // SUNDIALS + MPI
    // "examples/cpp/superlu/ex1p.cpp", // SuperLU_DIST + MPI
};

const c_example_main_sources: []const []const u8 = &.{
    "examples/c/ex0.c",
    "examples/c/ex1.c",
    "examples/c/ex2.c",
    "examples/c/ex3.c",
    "examples/c/ex4.c",
    "examples/c/ex5.c",
    "examples/c/ex6.c",
    "examples/c/ex7.c",
    "examples/c/ex8.c",
    "examples/c/ex9.c",
    "examples/c/ex10.c",
    "examples/c/ex14.c",
    "examples/c/ex15.c",
    "examples/c/ex16.c",
    "examples/c/ex17.c",
    "examples/c/ex18.c",
    "examples/c/ex19.c",
    "examples/c/ex20.c",
    "examples/c/ex21.c",
    "examples/c/ex22.c",
    "examples/c/ex23.c",
    "examples/c/ex24.c",
    "examples/c/ex25.c",
    "examples/c/ex26.c",
    "examples/c/ex27.c",
    "examples/c/ex28.c",
    "examples/c/ex29.c",
    "examples/c/ex30.c",
    "examples/c/ex31.c",
    "examples/c/ex33.c",
    "examples/c/ex34.c",
    "examples/c/ex36.c",
    "examples/c/ex37.c",
    "examples/c/ex38.c",
    "examples/c/ex39.c",
    "examples/c/ex40.c",
    "examples/c/ex41.c",
};

const miniapp_main_sources: []const []const u8 = &.{
    // "miniapps/adjoint/adjoint_advection_diffusion.cpp", // SUNDIALS
    // "miniapps/adjoint/cvsRoberts_ASAi_dns.cpp", // SUNDIALS
    // "miniapps/autodiff/par_example.cpp", // MPI
    "miniapps/autodiff/seq_example.cpp",
    "miniapps/autodiff/seq_test.cpp",
    // "miniapps/contact/contact.cpp", // Axom/Tribol + solver TPLs
    // "miniapps/dfem/dfem-minimal-surface.cpp", // MPI + mfem::future
    // "miniapps/diag-smoothers/abs-l1-jacobi.cpp", // MPI/hypre
    // "miniapps/diag-smoothers/mg-abs-l1-jacobi.cpp", // MPI/hypre
    "miniapps/dpg/acoustics.cpp",
    "miniapps/dpg/convection-diffusion.cpp",
    "miniapps/dpg/diffusion.cpp",
    "miniapps/dpg/maxwell.cpp",
    // "miniapps/dpg/pacoustics.cpp", // MPI
    // "miniapps/dpg/pconvection-diffusion.cpp", // MPI
    // "miniapps/dpg/pdiffusion.cpp", // MPI
    // "miniapps/dpg/pmaxwell.cpp", // MPI
    // "miniapps/electromagnetics/joule.cpp", // MPI
    // "miniapps/electromagnetics/lorentz.cpp", // MPI + GSLIB
    // "miniapps/electromagnetics/maxwell.cpp", // MPI
    // "miniapps/electromagnetics/tesla.cpp", // MPI
    // "miniapps/electromagnetics/volta.cpp", // MPI
    // "miniapps/fluids/navier/navier_3dfoc.cpp", // MPI
    // "miniapps/fluids/navier/navier_bifurcation.cpp", // MPI
    // "miniapps/fluids/navier/navier_cht.cpp", // MPI
    // "miniapps/fluids/navier/navier_kovasznay.cpp", // MPI
    // "miniapps/fluids/navier/navier_kovasznay_vs.cpp", // MPI
    // "miniapps/fluids/navier/navier_mms.cpp", // MPI
    // "miniapps/fluids/navier/navier_shear.cpp", // MPI
    // "miniapps/fluids/navier/navier_tgv.cpp", // MPI
    // "miniapps/fluids/navier/navier_turbchan.cpp", // MPI
    // "miniapps/fluids/schrodinger-flow/pschrodinger_flow.cpp", // MPI
    "miniapps/fluids/schrodinger-flow/schrodinger_flow.cpp",
    // "miniapps/gslib/field-diff.cpp", // GSLIB
    // "miniapps/gslib/field-interp.cpp", // GSLIB
    // "miniapps/gslib/findpts.cpp", // GSLIB
    // "miniapps/gslib/particles_redist.cpp", // GSLIB
    // "miniapps/gslib/pfindpts.cpp", // MPI
    // "miniapps/gslib/schwarz_ex1.cpp", // GSLIB
    // "miniapps/gslib/schwarz_ex1p.cpp", // MPI
    // "miniapps/hdiv-linear-solver/darcy.cpp", // MPI/hypre
    // "miniapps/hdiv-linear-solver/grad_div.cpp", // MPI/hypre
    // "miniapps/hooke/hooke.cpp", // MPI
    "miniapps/meshing/extruder.cpp",
    // "miniapps/meshing/fit-node-position.cpp", // MPI
    "miniapps/meshing/hpref.cpp",
    "miniapps/meshing/klein-bottle.cpp",
    // "miniapps/meshing/mesh-bounding-boxes.cpp", // MPI
    "miniapps/meshing/mesh-explorer.cpp",
    "miniapps/meshing/mesh-optimizer.cpp",
    "miniapps/meshing/mesh-quality.cpp",
    "miniapps/meshing/minimal-surface.cpp",
    "miniapps/meshing/mobius-strip.cpp",
    // "miniapps/meshing/phpref.cpp", // MPI
    // "miniapps/meshing/pmesh-fitting.cpp", // MPI
    // "miniapps/meshing/pmesh-optimizer.cpp", // MPI
    // "miniapps/meshing/pminimal-surface.cpp", // MPI
    "miniapps/meshing/polar-nc.cpp",
    "miniapps/meshing/ref321.cpp",
    "miniapps/meshing/reflector.cpp",
    "miniapps/meshing/shaper.cpp",
    "miniapps/meshing/toroid.cpp",
    "miniapps/meshing/trimmer.cpp",
    "miniapps/meshing/twist.cpp",
    // "miniapps/mtop/mtop_test_iso_elasticity.cpp", // MPI/hypre/future
    // "miniapps/multidomain/multidomain.cpp", // MPI/hypre
    // "miniapps/multidomain/multidomain_nd.cpp", // MPI/hypre
    // "miniapps/multidomain/multidomain_rt.cpp", // MPI/hypre
    "miniapps/nurbs/nurbs_curveint.cpp",
    "miniapps/nurbs/nurbs_ex1.cpp",
    "miniapps/nurbs/nurbs_ex10.cpp",
    // "miniapps/nurbs/nurbs_ex10p.cpp", // MPI
    // "miniapps/nurbs/nurbs_ex11p.cpp", // MPI
    // "miniapps/nurbs/nurbs_ex1p.cpp", // MPI
    "miniapps/nurbs/nurbs_ex24.cpp",
    "miniapps/nurbs/nurbs_ex3.cpp",
    "miniapps/nurbs/nurbs_ex5.cpp",
    "miniapps/nurbs/nurbs_mesh_info.cpp",
    "miniapps/nurbs/nurbs_naca_cmesh.cpp",
    "miniapps/nurbs/nurbs_patch_ex1.cpp",
    "miniapps/nurbs/nurbs_printfunc.cpp",
    "miniapps/nurbs/nurbs_solenoidal.cpp",
    "miniapps/nurbs/nurbs_surface.cpp",
    // "miniapps/parelag/MultilevelHcurlHdivSolver.cpp", // ParELAG
    "miniapps/performance/ex1.cpp",
    // "miniapps/performance/ex1p.cpp", // MPI
    // "miniapps/shifted/diffusion.cpp", // MPI/hypre
    // "miniapps/shifted/distance.cpp", // MPI/hypre
    // "miniapps/shifted/extrapolate.cpp", // MPI/hypre
    // "miniapps/shifted/lsf_integral.cpp", // MPI/hypre
    // "miniapps/solvers/block-solvers.cpp", // MPI/hypre
    // "miniapps/solvers/lor_elast.cpp", // MPI/hypre
    // "miniapps/solvers/lor_solvers.cpp", // MPI/hypre
    // "miniapps/solvers/plor_solvers.cpp", // MPI
    // "miniapps/spde/generate_random_field.cpp", // MPI/hypre
    "miniapps/tools/compare-dc.cpp",
    "miniapps/tools/convert-dc.cpp",
    "miniapps/tools/display-basis.cpp",
    "miniapps/tools/get-values.cpp",
    // "miniapps/tools/gridfunction-bounds.cpp", // MPI
    "miniapps/tools/load-dc.cpp",
    "miniapps/tools/lor-transfer.cpp",
    // "miniapps/tools/nodal-transfer.cpp", // MPI
    // "miniapps/tools/plor-transfer.cpp", // MPI
    "miniapps/tools/tmop-check-metric.cpp",
    "miniapps/tools/tmop-metric-magnitude.cpp",
    "miniapps/toys/automata.cpp",
    "miniapps/toys/life.cpp",
    "miniapps/toys/lissajous.cpp",
    "miniapps/toys/mandel.cpp",
    "miniapps/toys/mondrian.cpp",
    "miniapps/toys/rubik.cpp",
    "miniapps/toys/snake.cpp",
    "miniapps/toys/spiral.cpp",
    // "miniapps/tribol/contact-patch-test.cpp", // Tribol
};

const unit_test_skip_prefixes: []const []const u8 = &.{
    "tests/unit/ceed/", // libCEED
};

const TestSpec = struct {
    artifact_name: []const u8,
    install_name: []const u8,
    main_source: []const u8,
    extra_sources: []const []const u8 = &.{},
};

const InstalledExecutable = struct {
    step: *std.Build.Step.InstallArtifact,
    path: []const u8,
};

pub fn build(b: *std.Build) void {
    // Resolve the user-selected target triple and optimization mode.
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Expose the supported MFEM build toggles as Zig build options.
    const shared = b.option(bool, "shared", "Build MFEM as a shared library") orelse false;
    const precision = b.option(Precision, "precision", "Floating-point precision") orelse .double;
    const thread_safe = b.option(bool, "thread-safe", "Enable MFEM thread safety") orelse false;
    const use_memalloc = b.option(bool, "memalloc", "Enable MFEM internal MEMALLOC") orelse true;
    const use_simd = b.option(bool, "simd", "Enable MFEM SIMD code paths") orelse false;
    const arkode_dep = b.dependency("arkode_zig", .{
        .target = target,
        .optimize = optimize,
    });

    // Generate MFEM's config header from the selected Zig build options.
    const write_files = b.addWriteFiles();
    const io = b.graph.io;
    _ = write_files.add("config/_config.hpp", mfemConfigHeader(
        b,
        target,
        optimize,
        shared,
        precision,
        thread_safe,
        use_memalloc,
        use_simd,
    ));
    const generated_config_dir = write_files.getDirectory();
    const c_api_sources = collectCppSourcesUnder(b, io, "c_api") catch @panic("OOM");
    const c_adapter_sources = collectCppSourcesUnder(b, io, "examples/c/adapters") catch @panic("OOM");

    // Install the repository data files and test fixtures used by the
    // installed examples, miniapps, and tests.
    const install_repo_data = b.addInstallDirectory(.{
        .source_dir = b.path("data"),
        .install_dir = .prefix,
        .install_subdir = "share/data",
    });
    b.getInstallStep().dependOn(&install_repo_data.step);
    const install_test_fixture_data = b.addInstallDirectory(.{
        .source_dir = b.path("tests/unit/data"),
        .install_dir = .bin,
        .install_subdir = "tests/data",
    });
    b.getInstallStep().dependOn(&install_test_fixture_data.step);
    const install_miniapps_meshing_fixtures = b.addInstallDirectory(.{
        .source_dir = b.path("miniapps/meshing"),
        .install_dir = .prefix,
        .install_subdir = "share/miniapps/meshing",
        .include_extensions = &.{ ".mesh", ".vtk" },
    });
    b.getInstallStep().dependOn(&install_miniapps_meshing_fixtures.step);
    const install_miniapps_nurbs_fixtures = b.addInstallDirectory(.{
        .source_dir = b.path("miniapps/nurbs/meshes"),
        .install_dir = .prefix,
        .install_subdir = "share/miniapps/nurbs/meshes",
        .include_extensions = &.{".mesh"},
    });
    b.getInstallStep().dependOn(&install_miniapps_nurbs_fixtures.step);
    const install_miniapps_multidomain_fixtures = b.addInstallDirectory(.{
        .source_dir = b.path("miniapps/multidomain"),
        .install_dir = .prefix,
        .install_subdir = "share/miniapps/multidomain",
        .include_extensions = &.{".mesh"},
    });
    b.getInstallStep().dependOn(&install_miniapps_multidomain_fixtures.step);

    // Build the main MFEM library from the selected core source groups.
    const mfem = b.addLibrary(.{
        .name = "mfem",
        .linkage = if (shared) .dynamic else .static,
        .root_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libc = true,
            .link_libcpp = true,
        }),
    });

    configureMfemModule(mfem.root_module, b, generated_config_dir, arkode_dep);

    if (shared and target.result.os.tag == .windows) {
        mfem.root_module.addCMacro("mfem_EXPORTS", "1");
    }

    addSourceGroup(mfem.root_module, b, general_sources, cxx_flags);
    addSourceGroup(mfem.root_module, b, linalg_sources, cxx_flags);
    mfem.root_module.linkLibrary(arkode_dep.artifact("arkode"));
    addSourceGroup(mfem.root_module, b, mesh_sources, cxx_flags);
    addSourceGroup(mfem.root_module, b, fem_sources, cxx_flags);
    addSourceGroup(mfem.root_module, b, runtime_c_sources, c_flags);

    // Install the library so downstream build steps can depend on the final
    // zig-out layout rather than cache-only artifacts.
    const install_mfem = b.addInstallArtifact(mfem, .{});
    b.getInstallStep().dependOn(&install_mfem.step);

    // Build the C wrapper library on top of the C++ MFEM library.
    const mfem_c_api = b.addLibrary(.{
        .name = "mfem_c_api",
        .linkage = .static,
        .root_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libc = true,
            .link_libcpp = true,
        }),
    });
    configureMfemModule(mfem_c_api.root_module, b, generated_config_dir, arkode_dep);
    mfem_c_api.root_module.addIncludePath(b.path("c_api"));
    mfem_c_api.root_module.addIncludePath(b.path("examples/c"));
    mfem_c_api.root_module.linkLibrary(mfem);
    addSourceGroup(mfem_c_api.root_module, b, c_api_sources, cxx_flags);
    addSourceGroup(mfem_c_api.root_module, b, c_adapter_sources, cxx_flags);

    // Export the public C wrapper surface as a translated Zig module so
    // downstream packages can `@import("cmfem")` instead of relying on the
    // deprecated `@cImport` builtin directly.
    const cmfem_translate = b.addTranslateC(.{
        .root_source_file = b.path("c_api/cmfem.h"),
        .target = target,
        .optimize = optimize,
    });
    const cmfem_module = cmfem_translate.addModule("cmfem");
    cmfem_module.linkLibrary(mfem);
    cmfem_module.linkLibrary(mfem_c_api);
    cmfem_module.link_libcpp = true;

    // Build and install the serial C++ examples.
    const examples_step = b.step("examples", "Build the serial example executables into zig-out/bin/examples/cpp");
    examples_step.dependOn(&install_mfem.step);
    examples_step.dependOn(&install_repo_data.step);
    for (example_main_sources) |main_source| {
        const install_name = exampleInstallName(b, main_source, "");
        const artifact_name = sanitizeArtifactName(b, "example", main_source) catch @panic("OOM");
        const exe = addCppExecutable(
            b,
            target,
            optimize,
            generated_config_dir,
            arkode_dep,
            mfem,
            artifact_name,
            &.{main_source},
            &.{},
        );
        const installed = installExecutableStem(b, exe, install_name);
        examples_step.dependOn(&installed.step.step);
    }

    // Build and install the serial C examples that exercise the CMFEM wrapper.
    const c_examples_step = b.step("c_examples", "Build the C example executables into zig-out/bin/examples/c");
    c_examples_step.dependOn(&install_mfem.step);
    c_examples_step.dependOn(&install_repo_data.step);
    for (c_example_main_sources) |main_source| {
        const install_name = exampleInstallName(b, main_source, "");
        const artifact_name = sanitizeArtifactName(b, "c_example", main_source) catch @panic("OOM");
        const exe = addCExecutable(
            b,
            target,
            optimize,
            generated_config_dir,
            arkode_dep,
            mfem,
            mfem_c_api,
            artifact_name,
            &.{main_source},
            &.{ "examples/c", "c_api" },
        );
        const installed = installExecutableStem(b, exe, install_name);
        c_examples_step.dependOn(&installed.step.step);
    }

    // Collect shared support sources that miniapps reuse across executables.
    const miniapps_common_sources = collectCppSourcesUnder(b, io, "miniapps/common") catch @panic("OOM");

    // Build and install the supported serial miniapps.
    const miniapps_step = b.step("miniapps", "Build the serial miniapp executables into zig-out/bin/miniapps");
    miniapps_step.dependOn(&install_mfem.step);
    miniapps_step.dependOn(&install_repo_data.step);
    for (miniapp_main_sources) |main_source| {
        var sources = std.array_list.Managed([]const u8).init(b.allocator);
        appendUniqueSource(&sources, main_source) catch @panic("OOM");
        appendUniqueSources(&sources, collectSupportSourcesForExecutable(b, io, main_source) catch @panic("OOM")) catch @panic("OOM");
        appendUniqueSources(&sources, miniapps_common_sources) catch @panic("OOM");

        const install_name = sourceStem(main_source);
        const artifact_name = sanitizeArtifactName(b, "miniapp", main_source) catch @panic("OOM");
        const exe = addCppExecutable(
            b,
            target,
            optimize,
            generated_config_dir,
            arkode_dep,
            mfem,
            artifact_name,
            sources.items,
            &.{},
        );
        const installed = installExecutable(b, exe, "miniapps", install_name);
        miniapps_step.dependOn(&installed.step.step);
    }

    // Build the supported tests, install them under zig-out/bin/tests, and run
    // each installed executable from that final location.
    const test_step = b.step("test", "Build zig-out/bin/tests and run the supported test executables");
    test_step.dependOn(&install_repo_data.step);
    test_step.dependOn(&install_test_fixture_data.step);
    test_step.dependOn(&install_miniapps_meshing_fixtures.step);
    test_step.dependOn(&install_miniapps_nurbs_fixtures.step);
    test_step.dependOn(&install_miniapps_multidomain_fixtures.step);
    const unit_test_sources = collectNonMainSourcesUnder(b, io, "tests/unit", unit_test_skip_prefixes) catch @panic("OOM");
    const test_specs = [_]TestSpec{
        .{
            .artifact_name = "test_unit_tests",
            .install_name = "unit_tests",
            .main_source = "tests/unit/unit_test_main.cpp",
            .extra_sources = unit_test_sources,
        },
        // "tests/unit/punit_test_main.cpp", // MPI
        // "tests/unit/gpu_unit_test_main.cpp", // GPU
        // "tests/unit/pgpu_unit_test_main.cpp", // GPU + MPI
        // "tests/unit/ceed/test_ceed_main.cpp", // libCEED
        .{
            .artifact_name = "test_debug_device_tests",
            .install_name = "debug_device_tests",
            .main_source = "tests/unit/miniapps/test_debug_device.cpp",
        },
        .{
            .artifact_name = "test_sedov_tests_cpu",
            .install_name = "sedov_tests_cpu",
            .main_source = "tests/unit/miniapps/test_sedov.cpp",
        },
        .{
            .artifact_name = "test_tmop_pa_tests_cpu",
            .install_name = "tmop_pa_tests_cpu",
            .main_source = "tests/unit/miniapps/test_tmop_pa.cpp",
        },
        .{
            .artifact_name = "test_rates",
            .install_name = "rates",
            .main_source = "tests/convergence/rates.cpp",
        },
        // "tests/convergence/prates.cpp", // MPI
        // "tests/mem_manager/dangling-aliases.cpp", // MPI/hypre
        // "tests/par-mesh-format/ex1p.cpp", // MPI
        // "tests/benchmarks/bench_assembly_levels.cpp", // google benchmark
        // "tests/benchmarks/bench_ceed.cpp", // libCEED + google benchmark
        // "tests/benchmarks/bench_dg_amr.cpp", // google benchmark
        // "tests/benchmarks/bench_elasticity.cpp", // google benchmark
        // "tests/benchmarks/bench_tmop.cpp", // google benchmark
        // "tests/benchmarks/bench_vector.cpp", // google benchmark
        // "tests/benchmarks/bench_virtuals.cpp", // google benchmark
    };

    for (test_specs) |spec| {
        var sources = std.array_list.Managed([]const u8).init(b.allocator);
        appendUniqueSource(&sources, spec.main_source) catch @panic("OOM");
        appendUniqueSources(&sources, spec.extra_sources) catch @panic("OOM");

        const exe = addCppExecutable(
            b,
            target,
            optimize,
            generated_config_dir,
            arkode_dep,
            mfem,
            spec.artifact_name,
            sources.items,
            &.{"tests/unit"},
        );
        const installed = installExecutable(b, exe, "tests", spec.install_name);
        const run = b.addSystemCommand(&.{installed.path});
        run.step.dependOn(&install_mfem.step);
        run.step.dependOn(&install_repo_data.step);
        run.step.dependOn(&install_test_fixture_data.step);
        run.step.dependOn(&install_miniapps_meshing_fixtures.step);
        run.step.dependOn(&install_miniapps_nurbs_fixtures.step);
        run.step.dependOn(&install_miniapps_multidomain_fixtures.step);
        run.step.dependOn(&installed.step.step);
        run.setCwd(b.path("zig-out/bin/tests"));
        run.stdio = .inherit;
        test_step.dependOn(&run.step);
    }
}

// Apply the include path and config macro setup shared by every MFEM artifact.
fn configureMfemModule(
    module: *std.Build.Module,
    b: *std.Build,
    generated_config_dir: std.Build.LazyPath,
    arkode_dep: *std.Build.Dependency,
) void {
    module.addIncludePath(generated_config_dir);
    module.addIncludePath(b.path("config"));
    module.addIncludePath(b.path("core"));
    module.addIncludePath(b.path("."));
    module.addCMacro("MFEM_CONFIG_FILE", "\"config/_config.hpp\"");
    module.addIncludePath(arkode_dep.artifact("arkode").getEmittedIncludeTree());
    module.addCMacro("SUNDIALS_STATIC_DEFINE", "1");
}

// Create a C++ executable that links against MFEM and installs runtime path
// helpers via a forced include.
fn addCppExecutable(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    generated_config_dir: std.Build.LazyPath,
    arkode_dep: *std.Build.Dependency,
    mfem: *std.Build.Step.Compile,
    artifact_name: []const u8,
    sources: []const []const u8,
    extra_include_paths: []const []const u8,
) *std.Build.Step.Compile {
    const exe = b.addExecutable(.{
        .name = artifact_name,
        .root_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libc = true,
            .link_libcpp = true,
        }),
    });

    configureMfemModule(exe.root_module, b, generated_config_dir, arkode_dep);
    for (extra_include_paths) |path| {
        exe.root_module.addIncludePath(b.path(path));
    }
    exe.root_module.linkLibrary(mfem);
    addSourceGroup(exe.root_module, b, sources, executable_cxx_flags);

    return exe;
}

// Create a C executable that links against both MFEM and the CMFEM wrapper.
fn addCExecutable(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    generated_config_dir: std.Build.LazyPath,
    arkode_dep: *std.Build.Dependency,
    mfem: *std.Build.Step.Compile,
    mfem_c_api: *std.Build.Step.Compile,
    artifact_name: []const u8,
    sources: []const []const u8,
    extra_include_paths: []const []const u8,
) *std.Build.Step.Compile {
    const exe = b.addExecutable(.{
        .name = artifact_name,
        .root_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
            .link_libc = true,
            .link_libcpp = true,
        }),
    });

    configureMfemModule(exe.root_module, b, generated_config_dir, arkode_dep);
    for (extra_include_paths) |path| {
        exe.root_module.addIncludePath(b.path(path));
    }
    exe.root_module.linkLibrary(mfem);
    exe.root_module.linkLibrary(mfem_c_api);
    addSourceGroup(exe.root_module, b, sources, c_flags);

    return exe;
}

// Add a set of C/C++ translation units to an existing Zig module.
fn addSourceGroup(
    module: *std.Build.Module,
    b: *std.Build,
    files: []const []const u8,
    flags: []const []const u8,
) void {
    if (files.len == 0) {
        return;
    }

    module.addCSourceFiles(.{
        .root = b.path("."),
        .files = files,
        .flags = flags,
    });
}

// Install an executable under zig-out/bin using an extensionless stem path.
fn installExecutableStem(
    b: *std.Build,
    exe: *std.Build.Step.Compile,
    install_stem: []const u8,
) InstalledExecutable {
    const rel_path = b.fmt("{s}{s}", .{
        install_stem,
        std.fs.path.extension(exe.out_filename),
    });
    const step = b.addInstallArtifact(exe, .{
        .dest_sub_path = rel_path,
    });
    return .{
        .step = step,
        .path = b.getInstallPath(.bin, rel_path),
    };
}

// Install an executable under zig-out/bin/<group>/<name>.
fn installExecutable(
    b: *std.Build,
    exe: *std.Build.Step.Compile,
    group: []const u8,
    install_name: []const u8,
) InstalledExecutable {
    const rel_path = executableInstallRelPath(b, exe, group, install_name);
    return installExecutableStem(b, exe, rel_path[0 .. rel_path.len - std.fs.path.extension(rel_path).len]);
}

// Compute the installed relative path for an executable within a group.
fn executableInstallRelPath(
    b: *std.Build,
    exe: *std.Build.Step.Compile,
    group: []const u8,
    install_name: []const u8,
) []const u8 {
    return b.fmt("{s}/{s}{s}", .{
        group,
        install_name,
        std.fs.path.extension(exe.out_filename),
    });
}

// Strip the extension from a source path to derive the installed executable name.
fn sourceStem(source_path: []const u8) []const u8 {
    const basename = std.fs.path.basename(source_path);
    const extension = std.fs.path.extension(basename);
    return basename[0 .. basename.len - extension.len];
}

// Derive an installed example path stem relative to the source prefix, e.g.
// examples/cpp/sundials/ex9.cpp -> examples/cpp/sundials/ex9.
fn exampleInstallName(
    b: *std.Build,
    source_path: []const u8,
    prefix: []const u8,
) []const u8 {
    const relative_path = if (std.mem.startsWith(u8, source_path, prefix))
        source_path[prefix.len..]
    else
        source_path;
    const extension = std.fs.path.extension(relative_path);
    const stem = relative_path[0 .. relative_path.len - extension.len];

    return b.dupe(stem);
}

// Recursively collect C++ sources under a directory.
fn appendCppSourcesUnder(
    b: *std.Build,
    io: std.Io,
    list: *std.array_list.Managed([]const u8),
    dir_path: []const u8,
) !void {
    var dir = try std.Io.Dir.cwd().openDir(io, b.pathFromRoot(dir_path), .{ .iterate = true });
    defer dir.close(io);

    var iterator = dir.iterate();
    while (try iterator.next(io)) |entry| {
        const child_path = try std.fs.path.join(b.allocator, &.{ dir_path, entry.name });
        switch (entry.kind) {
            .directory => try appendCppSourcesUnder(b, io, list, child_path),
            .file => {
                if (std.mem.endsWith(u8, entry.name, ".cpp")) {
                    try list.append(child_path);
                }
            },
            else => {},
        }
    }
}

// Return the sorted list of C++ sources under a directory.
fn collectCppSourcesUnder(
    b: *std.Build,
    io: std.Io,
    dir_path: []const u8,
) ![]const []const u8 {
    var list = std.array_list.Managed([]const u8).init(b.allocator);
    try appendCppSourcesUnder(b, io, &list, dir_path);
    std.mem.sort([]const u8, list.items, {}, lessThanString);
    return list.toOwnedSlice();
}

// Collect support sources for an executable by taking sibling .cpp files that
// do not define their own main function.
fn collectSupportSourcesForExecutable(
    b: *std.Build,
    io: std.Io,
    main_source: []const u8,
) ![]const []const u8 {
    const dir_path = std.fs.path.dirname(main_source).?;
    const dir_sources = try collectCppSourcesUnder(b, io, dir_path);

    var list = std.array_list.Managed([]const u8).init(b.allocator);
    for (dir_sources) |path| {
        if (std.mem.eql(u8, path, main_source)) {
            continue;
        }
        if (try fileHasMain(b, io, path)) continue;
        try list.append(path);
    }
    std.mem.sort([]const u8, list.items, {}, lessThanString);
    return list.toOwnedSlice();
}

// Collect all non-main C++ sources under a directory while skipping unsupported
// subtrees such as libCEED tests.
fn collectNonMainSourcesUnder(
    b: *std.Build,
    io: std.Io,
    dir_path: []const u8,
    skip_prefixes: []const []const u8,
) ![]const []const u8 {
    const dir_sources = try collectCppSourcesUnder(b, io, dir_path);

    var list = std.array_list.Managed([]const u8).init(b.allocator);
    outer: for (dir_sources) |path| {
        for (skip_prefixes) |skip_prefix| {
            if (std.mem.startsWith(u8, path, skip_prefix)) {
                continue :outer;
            }
        }
        if (try fileHasMain(b, io, path)) continue;
        try list.append(path);
    }
    std.mem.sort([]const u8, list.items, {}, lessThanString);
    return list.toOwnedSlice();
}

// Detect whether a translation unit defines a main function and should
// therefore be treated as a standalone executable entry point.
fn fileHasMain(
    b: *std.Build,
    io: std.Io,
    path: []const u8,
) !bool {
    const contents = try std.Io.Dir.cwd().readFileAlloc(
        io,
        b.pathFromRoot(path),
        b.allocator,
        .limited(16 * 1024 * 1024),
    );
    defer b.allocator.free(contents);

    return std.mem.indexOf(u8, contents, "int main(") != null or
        std.mem.indexOf(u8, contents, "int main (") != null;
}

// Append a source path only if it has not already been added to the list.
fn appendUniqueSource(
    list: *std.array_list.Managed([]const u8),
    item: []const u8,
) !void {
    for (list.items) |existing| {
        if (std.mem.eql(u8, existing, item)) {
            return;
        }
    }
    try list.append(item);
}

// Append a list of sources while preserving uniqueness.
fn appendUniqueSources(
    list: *std.array_list.Managed([]const u8),
    items: []const []const u8,
) !void {
    for (items) |item| {
        try appendUniqueSource(list, item);
    }
}

// Convert a source path into a stable artifact name that Zig accepts.
fn sanitizeArtifactName(
    b: *std.Build,
    prefix: []const u8,
    source_path: []const u8,
) ![]const u8 {
    var name = std.array_list.Managed(u8).init(b.allocator);
    try name.appendSlice(prefix);
    try name.append('_');
    for (source_path) |ch| {
        switch (ch) {
            '/', '\\', '-', '.' => try name.append('_'),
            else => try name.append(ch),
        }
    }
    return name.toOwnedSlice();
}

// Lexicographic comparator used when sorting source path lists.
fn lessThanString(_: void, lhs: []const u8, rhs: []const u8) bool {
    return std.mem.lessThan(u8, lhs, rhs);
}

// Generate MFEM's config header so the Zig build produces the same compile-time
// feature macros that the C++ sources expect from other build systems.
fn mfemConfigHeader(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    shared: bool,
    precision: Precision,
    thread_safe: bool,
    use_memalloc: bool,
    use_simd: bool,
) []const u8 {
    const timer_type: u8 = switch (target.result.os.tag) {
        .windows => 3,
        .macos, .ios, .tvos, .watchos, .visionos => 4,
        else => 2,
    };

    return b.fmt(
        \\#ifndef MFEM_CONFIG_HEADER
        \\#define MFEM_CONFIG_HEADER
        \\
        \\#define MFEM_VERSION 40901
        \\#define MFEM_VERSION_STRING "4.9.1"
        \\#define MFEM_VERSION_MAJOR ((MFEM_VERSION) / 10000)
        \\#define MFEM_VERSION_MINOR (((MFEM_VERSION) / 100) % 100)
        \\#define MFEM_VERSION_PATCH ((MFEM_VERSION) % 100)
        \\#define MFEM_SOURCE_DIR "."
        \\#define MFEM_INSTALL_DIR "zig-out"
        \\#define MFEM_GIT_STRING "(unknown)"
        \\{s}
        \\{s}
        \\{s}
        \\{s}
        \\{s}
        \\{s}
        \\{s}
        \\#define MFEM_TIMER_TYPE {d}
        \\
        \\#endif // MFEM_CONFIG_HEADER
        \\
    ,
        .{
            if (shared) "#define MFEM_SHARED_BUILD" else "",
            switch (precision) {
                .double => "#define MFEM_USE_DOUBLE",
                .single => "#define MFEM_USE_SINGLE",
            },
            if (optimize == .Debug) "#define MFEM_DEBUG" else "",
            if (thread_safe) "#define MFEM_THREAD_SAFE" else "",
            if (use_memalloc) "#define MFEM_USE_MEMALLOC" else "",
            if (use_simd) "#define MFEM_USE_SIMD" else "",
            "#define MFEM_USE_SUNDIALS\n#define MFEM_USE_SUNDIALS_ARKODE",
            timer_type,
        },
    );
}
