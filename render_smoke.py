import os
import subprocess
import sys


def main() -> int:
    outdir = os.environ.get("GENLAPROBLEMS_SMOKE_OUT", "/tmp/la/smoke/genlaproblems")
    os.makedirs(outdir, exist_ok=True)
    print("GenLAProblems smoke ->", outdir)
    code = r"""
    using GenLAProblems
    using LAlatex
    A, X, B = gen_gj_pb(3, 4, 3; maxint=3, num_rhs=2)
    println(size(A), size(X), size(B))
    println(length(String(LAlatex.L_show("A = ", A))))
    """
    env = dict(os.environ)
    env["GENLAPROBLEMS_SMOKE_OUT"] = outdir
    cmd = ["julia", "-e", code]
    return subprocess.call(cmd, env=env)


if __name__ == "__main__":
    raise SystemExit(main())
