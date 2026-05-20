using GenLAProblems
using LAlatex
using LaTeXStrings

warmup_dir = "/tmp/genla_warmup"
mkpath(warmup_dir)

A_ge, X_ge, B_ge = gen_gj_pb(3, 4, 3; maxint=3, num_rhs=2)
display(l_show("A X = B, ", A_ge, X_ge, " = ", B_ge))

A, A_inv = gen_inv_pb(3; maxint=2)
display(l_show("A = ", A, ", ", L"A^{-1} = ", A_inv))

pivot_cols_lu, L_lu, U_lu, A_lu = gen_lu_pb(3, 3, 3; maxint=2)
display(l_show("A = LU, ", A_lu, " = ", L_lu, U_lu))

A_qr = gen_qr_problem(4; family=:pythagorean, maxint=2)
display(l_show("A = ", A_qr))
