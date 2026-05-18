using LATeachingSuite, GenLAProblems

warmup_dir = "/tmp/genla_sysimage"
mkpath(warmup_dir)

GenLAProblems.import_sympy()
GenLAProblems.load_LAFigureSpecs()
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

la = GenLAProblems.load_LAFigureSpecs()
la.ge_bundle([1 0; 0 1])
la.eig_bundle([1 0; 0 1])
la.svd_bundle([2 0; 0 1])
la.qr_bundle([1 0; 0 1])
la.qr_figure([1 0; 0 1])
