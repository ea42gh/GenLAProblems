using Test
using GenLAProblems
using PythonCall

function _py_keys_set(d)
    keys_obj = GenLAProblems._pycall(PythonCall.pygetattr(d, "keys"))
    keys_list = GenLAProblems._pycall(keys_obj)
    return Set(String.(collect(PythonCall.pylist(keys_list))))
end

@testset "bundle spec parity with LAFigureSpecs" begin
    GenLAProblems._ensure_pythoncall()
    la = GenLAProblems.load_LAFigureSpecs()
    A = [1 0; 0 1]

    for (bundle_sym, jl_fn) in [
        (:ge_tbl_bundle, GenLAProblems.ge_tbl_bundle),
        (:eig_tbl_bundle, GenLAProblems.eig_tbl_bundle),
        (:svd_tbl_bundle, GenLAProblems.svd_tbl_bundle),
        (:qr_tbl_bundle, GenLAProblems.qr_tbl_bundle),
    ]
        py_bundle = GenLAProblems._pycall(GenLAProblems._pygetattr(la, bundle_sym), A)
        py_spec = py_bundle["spec"]
        _svg, jl_spec = jl_fn(A)
        @test _py_keys_set(py_spec) == _py_keys_set(jl_spec)
    end
end
