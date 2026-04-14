const std = @import("std");

const Precision = enum {
    double,
    single,
};

const general_sources = &.{
    "general/array.cpp",
    "general/binaryio.cpp",
    "general/cuda.cpp",
    "general/device.cpp",
    "general/error.cpp",
    "general/gecko.cpp",
    "general/globals.cpp",
    "general/hash.cpp",
    "general/hash_util.cpp",
    "general/isockstream.cpp",
    "general/mem_manager.cpp",
    "general/occa.cpp",
    "general/optparser.cpp",
    "general/osockstream.cpp",
    "general/sets.cpp",
    "general/socketstream.cpp",
    "general/stable3d.cpp",
    "general/table.cpp",
    "general/tic_toc.cpp",
    "general/tinyxml2.cpp",
    "general/version.cpp",
    "general/hip.cpp",
};

const linalg_sources = &.{
    "linalg/auxiliary.cpp",
    "linalg/batched/batched.cpp",
    "linalg/batched/gpu_blas.cpp",
    "linalg/batched/magma.cpp",
    "linalg/batched/native.cpp",
    "linalg/batched/solver.cpp",
    "linalg/blockmatrix.cpp",
    "linalg/blockoperator.cpp",
    "linalg/blockvector.cpp",
    "linalg/complex_densemat.cpp",
    "linalg/complex_operator.cpp",
    "linalg/constraints.cpp",
    "linalg/densemat.cpp",
    "linalg/filteredsolver.cpp",
    "linalg/handle.cpp",
    "linalg/matrix.cpp",
    "linalg/mma.cpp",
    "linalg/ode.cpp",
    "linalg/operator.cpp",
    "linalg/ordering.cpp",
    "linalg/particlevector.cpp",
    "linalg/solvers.cpp",
    "linalg/sparsemat.cpp",
    "linalg/sparsesmoothers.cpp",
    "linalg/symmat.cpp",
    "linalg/vector.cpp",
};

const mesh_sources = &.{
    "mesh/attribute_sets.cpp",
    "mesh/element.cpp",
    "mesh/exodus_writer.cpp",
    "mesh/face_nbr_geom.cpp",
    "mesh/gmsh.cpp",
    "mesh/hexahedron.cpp",
    "mesh/mesh.cpp",
    "mesh/mesh_operators.cpp",
    "mesh/mesh_readers.cpp",
    "mesh/ncmesh.cpp",
    "mesh/ncnurbs.cpp",
    "mesh/nurbs.cpp",
    "mesh/point.cpp",
    "mesh/pyramid.cpp",
    "mesh/quadrilateral.cpp",
    "mesh/segment.cpp",
    "mesh/spacing.cpp",
    "mesh/tetrahedron.cpp",
    "mesh/triangle.cpp",
    "mesh/vertex.cpp",
    "mesh/vtk.cpp",
    "mesh/vtkhdf.cpp",
    "mesh/wedge.cpp",
    "mesh/submesh/submesh.cpp",
    "mesh/submesh/ncsubmesh.cpp",
    "mesh/submesh/submesh_utils.cpp",
    "mesh/submesh/transfermap.cpp",
};

const fem_sources = &.{
    "fem/bilinearform.cpp",
    "fem/bilinearform_ext.cpp",
    "fem/bilininteg.cpp",
    "fem/integ/bilininteg_br2.cpp",
    "fem/integ/bilininteg_convection_mf.cpp",
    "fem/integ/bilininteg_convection_pa.cpp",
    "fem/integ/bilininteg_convection_ea.cpp",
    "fem/integ/bilininteg_curlcurl_pa.cpp",
    "fem/integ/bilininteg_dgdiffusion_pa.cpp",
    "fem/integ/bilininteg_dgtrace_pa.cpp",
    "fem/integ/bilininteg_dgtrace_ea.cpp",
    "fem/integ/bilininteg_diffusion_mf.cpp",
    "fem/integ/bilininteg_diffusion_pa.cpp",
    "fem/integ/bilininteg_diffusion_ea.cpp",
    "fem/integ/bilininteg_diffusion_patch.cpp",
    "fem/integ/bilininteg_divdiv_pa.cpp",
    "fem/integ/bilininteg_elasticity_ea.cpp",
    "fem/integ/bilininteg_elasticity_pa.cpp",
    "fem/integ/bilininteg_gradient_pa.cpp",
    "fem/integ/bilininteg_interp_pa.cpp",
    "fem/integ/bilininteg_mass_mf.cpp",
    "fem/integ/bilininteg_mass_pa.cpp",
    "fem/integ/bilininteg_mass_ea.cpp",
    "fem/integ/bilininteg_mixedcurl_pa.cpp",
    "fem/integ/bilininteg_mixedvecgrad_pa.cpp",
    "fem/integ/bilininteg_trace_jump_ea.cpp",
    "fem/integ/bilininteg_transpose_ea.cpp",
    "fem/integ/bilininteg_vecdiffusion_mf.cpp",
    "fem/integ/bilininteg_vecdiffusion_pa.cpp",
    "fem/integ/bilininteg_vecdiv_pa.cpp",
    "fem/integ/bilininteg_vecmass_mf.cpp",
    "fem/integ/bilininteg_vecmass_pa.cpp",
    "fem/integ/bilininteg_vectorfediv_pa.cpp",
    "fem/integ/bilininteg_vectorfemass_pa.cpp",
    "fem/integ/bilininteg_diffusion_kernels.cpp",
    "fem/integ/bilininteg_elasticity_kernels.cpp",
    "fem/integ/bilininteg_hcurl_kernels.cpp",
    "fem/integ/bilininteg_hdiv_ea.cpp",
    "fem/integ/bilininteg_hdiv_kernels.cpp",
    "fem/integ/bilininteg_hcurlhdiv_kernels.cpp",
    "fem/integ/bilininteg_mass_kernels.cpp",
    "fem/integ/lininteg_boundary.cpp",
    "fem/integ/lininteg_boundary_flux.cpp",
    "fem/integ/lininteg_domain.cpp",
    "fem/integ/lininteg_domain_grad.cpp",
    "fem/integ/lininteg_domain_vectorfe.cpp",
    "fem/integ/nonlininteg_vecconvection_pa.cpp",
    "fem/integ/nonlininteg_vecconvection_mf.cpp",
    "fem/coefficient.cpp",
    "fem/complex_fem.cpp",
    "fem/convergence.cpp",
    "fem/datacollection.cpp",
    "fem/dgmassinv.cpp",
    "fem/doftrans.cpp",
    "fem/dfem/doperator.cpp",
    "fem/eltrans.cpp",
    "fem/batchitrans.cpp",
    "fem/estimators.cpp",
    "fem/fe.cpp",
    "fem/fe/face_map_utils.cpp",
    "fem/fe/fe_base.cpp",
    "fem/fe/fe_fixed_order.cpp",
    "fem/fe/fe_h1.cpp",
    "fem/fe/fe_l2.cpp",
    "fem/fe/fe_nd.cpp",
    "fem/fe/fe_nurbs.cpp",
    "fem/fe/fe_pos.cpp",
    "fem/fe/fe_pyramid.cpp",
    "fem/fe/fe_rt.cpp",
    "fem/fe/fe_ser.cpp",
    "fem/fe_coll.cpp",
    "fem/fespace.cpp",
    "fem/derefmat_op.cpp",
    "fem/pderefmat_op.cpp",
    "fem/geom.cpp",
    "fem/gridfunc.cpp",
    "fem/hybridization.cpp",
    "fem/hybridization_ext.cpp",
    "fem/intrules.cpp",
    "fem/intrules_cut.cpp",
    "fem/ceed/interface/basis.cpp",
    "fem/ceed/interface/restriction.cpp",
    "fem/ceed/interface/operator.cpp",
    "fem/ceed/interface/util.cpp",
    "fem/ceed/integrators/convection/convection.cpp",
    "fem/ceed/integrators/diffusion/diffusion.cpp",
    "fem/ceed/integrators/nlconvection/nlconvection.cpp",
    "fem/ceed/integrators/mass/mass.cpp",
    "fem/ceed/solvers/algebraic.cpp",
    "fem/ceed/solvers/full-assembly.cpp",
    "fem/ceed/solvers/solvers-atpmg.cpp",
    "fem/kdtree.cpp",
    "fem/linearform.cpp",
    "fem/linearform_ext.cpp",
    "fem/lininteg.cpp",
    "fem/lor/lor.cpp",
    "fem/lor/lor_ads.cpp",
    "fem/lor/lor_ams.cpp",
    "fem/lor/lor_batched.cpp",
    "fem/multigrid.cpp",
    "fem/nonlinearform.cpp",
    "fem/nonlinearform_ext.cpp",
    "fem/nonlininteg.cpp",
    "fem/fespacehierarchy.cpp",
    "fem/qfunction.cpp",
    "fem/qinterp/det.cpp",
    "fem/qinterp/eval_by_nodes.cpp",
    "fem/qinterp/eval_by_vdim.cpp",
    "fem/qinterp/eval_hdiv.cpp",
    "fem/qinterp/grad_by_nodes.cpp",
    "fem/qinterp/grad_by_vdim.cpp",
    "fem/qspace.cpp",
    "fem/quadinterpolator.cpp",
    "fem/quadinterpolator_face.cpp",
    "fem/restriction.cpp",
    "fem/normal_deriv_restriction.cpp",
    "fem/staticcond.cpp",
    "fem/tmop.cpp",
    "fem/tmop/pa.cpp",
    "fem/tmop/assemble/diag2_limit.cpp",
    "fem/tmop/assemble/diag2.cpp",
    "fem/tmop/assemble/grad2_limit.cpp",
    "fem/tmop/assemble/grad2.cpp",
    "fem/tmop/assemble/diag3_limit.cpp",
    "fem/tmop/assemble/diag3.cpp",
    "fem/tmop/assemble/grad3_limit.cpp",
    "fem/tmop/assemble/grad3.cpp",
    "fem/tmop/metrics/001.cpp",
    "fem/tmop/metrics/002.cpp",
    "fem/tmop/metrics/007.cpp",
    "fem/tmop/metrics/056.cpp",
    "fem/tmop/metrics/077.cpp",
    "fem/tmop/metrics/080.cpp",
    "fem/tmop/metrics/094.cpp",
    "fem/tmop/metrics/302.cpp",
    "fem/tmop/metrics/303.cpp",
    "fem/tmop/metrics/315.cpp",
    "fem/tmop/metrics/318.cpp",
    "fem/tmop/metrics/321.cpp",
    "fem/tmop/metrics/332.cpp",
    "fem/tmop/metrics/338.cpp",
    "fem/tmop/mult/grad2_limit.cpp",
    "fem/tmop/mult/grad2.cpp",
    "fem/tmop/mult/mult2_limit.cpp",
    "fem/tmop/mult/mult2.cpp",
    "fem/tmop/mult/grad3_limit.cpp",
    "fem/tmop/mult/grad3.cpp",
    "fem/tmop/mult/mult3_limit.cpp",
    "fem/tmop/mult/mult3.cpp",
    "fem/tmop/tools/det2_jpr.cpp",
    "fem/tmop/tools/det3_jpr.cpp",
    "fem/tmop/tools/discrete.cpp",
    "fem/tmop/tools/energy2_limit.cpp",
    "fem/tmop/tools/energy2.cpp",
    "fem/tmop/tools/energy3_limit.cpp",
    "fem/tmop/tools/energy3.cpp",
    "fem/tmop/tools/target2.cpp",
    "fem/tmop/tools/target3.cpp",
    "fem/tmop_tools.cpp",
    "fem/tmop_amr.cpp",
    "fem/gslib.cpp",
    "fem/gslib/findpts_local_2.cpp",
    "fem/gslib/findpts_local_3.cpp",
    "fem/gslib/interpolate_local_2.cpp",
    "fem/gslib/interpolate_local_3.cpp",
    "fem/transfer.cpp",
    "fem/hyperbolic.cpp",
    "fem/integrator.cpp",
    "fem/bounds.cpp",
    "fem/particleset.cpp",
};

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const shared = b.option(bool, "shared", "Build MFEM as a shared library") orelse false;
    const precision = b.option(Precision, "precision", "Floating-point precision") orelse .double;
    const use_exceptions = b.option(bool, "exceptions", "Enable MFEM exceptions") orelse false;
    const thread_safe = b.option(bool, "thread-safe", "Enable MFEM thread safety") orelse false;
    const use_memalloc = b.option(bool, "memalloc", "Enable MFEM internal MEMALLOC") orelse true;
    const use_simd = b.option(bool, "simd", "Enable MFEM SIMD code paths") orelse false;

    const write_files = b.addWriteFiles();
    _ = write_files.add("config/_config.hpp", mfemConfigHeader(
        b,
        target,
        optimize,
        shared,
        precision,
        use_exceptions,
        thread_safe,
        use_memalloc,
        use_simd,
    ));

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

    mfem.root_module.addIncludePath(write_files.getDirectory());
    mfem.root_module.addIncludePath(b.path("."));
    mfem.root_module.addCMacro("MFEM_CONFIG_FILE", "\"config/_config.hpp\"");

    if (shared and target.result.os.tag == .windows) {
        mfem.root_module.addCMacro("mfem_EXPORTS", "1");
    }

    addSourceGroup(mfem.root_module, b, general_sources);
    addSourceGroup(mfem.root_module, b, linalg_sources);
    addSourceGroup(mfem.root_module, b, mesh_sources);
    addSourceGroup(mfem.root_module, b, fem_sources);

    b.installArtifact(mfem);
}

fn addSourceGroup(
    module: *std.Build.Module,
    b: *std.Build,
    files: []const []const u8,
) void {
    module.addCSourceFiles(.{
        .root = b.path("."),
        .files = files,
        .flags = &.{
            "-std=c++17",
        },
    });
}

fn mfemConfigHeader(
    b: *std.Build,
    target: std.Build.ResolvedTarget,
    optimize: std.builtin.OptimizeMode,
    shared: bool,
    precision: Precision,
    use_exceptions: bool,
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
            if (use_exceptions) "#define MFEM_USE_EXCEPTIONS" else "",
            if (thread_safe) "#define MFEM_THREAD_SAFE" else "",
            if (use_memalloc) "#define MFEM_USE_MEMALLOC" else "",
            if (use_simd) "#define MFEM_USE_SIMD" else "",
            timer_type,
        },
    );
}
