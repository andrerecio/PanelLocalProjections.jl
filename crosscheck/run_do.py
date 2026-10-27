"""Run a Stata do-file headless through PyStata.

Usage (from the repository root):
    .venv/bin/python crosscheck/run_do.py crosscheck/stata/smoke.do [--log out.log] [--quiet]

The do-file runs with its own directory as Stata's working directory, and the global
macro PLP_ROOT holds the repository root, so do-files can reach the fixtures and the
project-local ado directory (see crosscheck/stata/common.do). The exit code is 0 on
success and the Stata return code r(#) on failure.

Adapted from LocalProjections.jl (andrerecio fork), stata/run_do.py; the change is the
PLP_ROOT global.
"""

import argparse
import contextlib
import io
import os
import re
import sys
from pathlib import Path

STATA_DIR = os.environ.get("STATA_DIR", "/Applications/StataNow")
STATA_EDITION = os.environ.get("STATA_EDITION", "be")
ROOT = Path(__file__).resolve().parent.parent


class _Tee(io.TextIOBase):
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)
        return len(s)

    def flush(self):
        for st in self.streams:
            st.flush()


def _return_code(msg: str) -> int:
    m = re.search(r"r\((\d+)\);", msg)
    return int(m.group(1)) if m else 1


def main() -> int:
    ap = argparse.ArgumentParser(description="Run a Stata do-file through PyStata.")
    ap.add_argument("dofile", type=Path)
    ap.add_argument("--log", type=Path, help="also write the Stata output to this file")
    ap.add_argument("--quiet", action="store_true", help="do not echo Stata output")
    args = ap.parse_args()

    dofile = args.dofile.resolve()
    if not dofile.is_file():
        print(f"no such do-file: {dofile}", file=sys.stderr)
        return 2
    log = args.log.resolve() if args.log else None

    import stata_setup

    stata_setup.config(STATA_DIR, STATA_EDITION, splash=False)
    from pystata import stata

    buf = io.StringIO()
    rc = 0
    sink = buf if args.quiet else _Tee(sys.stdout, buf)
    with contextlib.redirect_stdout(sink):
        try:
            stata.run(f'global PLP_ROOT "{ROOT}"')
            stata.run(f'quietly cd "{dofile.parent}"')
            stata.run(f'do "{dofile}"')
        except SystemError as e:  # PyStata raises SystemError on a non-zero r()
            rc = _return_code(str(e))
            print(f"\n[run_do] Stata error (rc={rc}):\n{e}", file=sys.stderr)

    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(buf.getvalue())
    return rc


if __name__ == "__main__":
    code = main()
    sys.stdout.flush()
    sys.stderr.flush()
    # The Stata engine's shutdown resets the exit status to 0, so bypass it.
    os._exit(code)
