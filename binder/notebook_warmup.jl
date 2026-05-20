using GenLAProblems, LAlatex, LinearAlgebra, LaTeXStrings, Random

Random.seed!(42)

A_ge, X_ge, B_ge = gen_gj_pb(3, 4, 3; maxint=3, num_rhs=2)
display(LAlatex.l_show(L"A X = B, \qquad = ", A_ge, X_ge, B_ge))

A, A_inv = gen_inv_pb(3; maxint=2)
display(LAlatex.l_show(L"A = ", A, L",\quad A^{-1} = ", A_inv))

pivot_cols_lu, L_lu, U_lu, A_lu = gen_lu_pb(3, 3, 3; maxint=2)
display(LAlatex.l_show(L"A = LU, \qquad = ", A_lu, L_lu, U_lu))

pivot_cols_plu, P_plu, L_plu, U_plu, A_plu = gen_plu_pb(3, 3, 3; maxint=2)
display(LAlatex.l_show(L"A = P L U, \qquad = ", A_plu, P_plu, L_plu, U_plu))

A_qr = gen_qr_problem(4; family=:pythagorean, maxint=2)
display(LAlatex.l_show(L"A = ", A_qr))

S_eig, Lambda_eig, S_inv_eig, A_eig = gen_eigenproblem([3, -1, 2]; maxint=2)
display(LAlatex.l_show(L"A = S \Lambda S^{-1}, \qquad = ", A_eig, S_eig, Lambda_eig, S_inv_eig))

Q_sym, Lambda_sym, A_sym = gen_symmetric_eigenproblem([4, 1, -2]; maxint=2)
display(LAlatex.l_show(L"A = Q \Lambda Q^T, \qquad = ", A_sym, Q_sym, Lambda_sym, transpose(Q_sym)))

U_svd, Sigma_svd, Vt_svd, A_svd = gen_svd_problem([2, 1], [2, 1], [3, 1, 0]; maxint=2)
display(LAlatex.l_show(L"A = U \Sigma V^T, \qquad = ", A_svd, U_svd, Sigma_svd, Vt_svd))
