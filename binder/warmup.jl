using Pkg

Pkg.instantiate()
Pkg.precompile()

using GenLAProblems

warmup_dir = "/tmp/genla_warmup"
mkpath(warmup_dir)

GenLAProblems.import_sympy()
GenLAProblems.load_la_figures()
GenLAProblems.load_matrixlayout()

A = Rational{Int}[1 2 0; 0 1 3; 0 0 2]
b = Rational{Int}[5; 7; 4]
pb = ShowGE{Rational{Int}}(A, b; output_dir=warmup_dir)
ref!(pb)
show_layout!(pb; fig_scale=1.1, output_dir=warmup_dir)
show_system(pb; b_col=1, output_dir=warmup_dir)
show_backsubstitution!(pb; b_col=1)
show_solution!(pb; b_col=1)

matrices, pivots, desc = reduce_to_ref(A; gj=true)
show_ge_final(matrices, desc, pivots; n_rhs=0, output_dir=warmup_dir)

nM.show_ge_tbl([1 0; 0 1]; tmp_dir=warmup_dir)
nM.show_eig_tbl([1 0; 0 1]; tmp_dir=warmup_dir)
nM.show_svd_tbl([2 0; 0 1]; tmp_dir=warmup_dir)
nM.qr_tbl_svg([1 0; 0 1]; tmp_dir=warmup_dir)
gram_schmidt_qr([1 0; 0 1]; tmp_dir=warmup_dir)
