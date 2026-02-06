import os
import subprocess
import sys


def main() -> int:
    outdir = os.environ.get("GENLAPROBLEMS_SMOKE_OUT", "/tmp/la/smoke/genlaproblems")
    os.makedirs(outdir, exist_ok=True)
    print("GenLAProblems smoke render ->", outdir)
    code = r"""
    using GenLAProblems
    A = [1 2; 3 4]
    svg = ge_tbl_svg(A; output_dir=ENV["GENLAPROBLEMS_SMOKE_OUT"])
    println(length(String(svg)))
    """
    env = dict(os.environ)
    env["GENLAPROBLEMS_SMOKE_OUT"] = outdir
    cmd = ["julia", "-e", code]
    return subprocess.call(cmd, env=env)


if __name__ == "__main__":
    raise SystemExit(main())
