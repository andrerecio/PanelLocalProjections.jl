"""Record or verify SHA-256 checksums of the project-local Stata ado tree.

    .venv/bin/python crosscheck/ado_manifest.py           # write crosscheck/ado_manifest.sha256
    .venv/bin/python crosscheck/ado_manifest.py --check   # fail if any installed file drifted

SSC packages cannot be pinned to a version, so a changed checksum means the upstream
package changed: re-baseline deliberately (rerun the cross-checks, review the diff of the
committed outputs) before committing a new manifest.
"""

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ADO = ROOT / "crosscheck" / "ado" / "plus"
MANIFEST = ROOT / "crosscheck" / "ado_manifest.sha256"
SKIP = {"backup.trk"}


def current() -> dict:
    out = {}
    for p in sorted(ADO.rglob("*")):
        if p.is_file() and p.name not in SKIP and p.name != ".DS_Store":
            rel = p.relative_to(ADO).as_posix()
            if rel.endswith(".trk"):
                continue  # install log with timestamps
            out[rel] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def main() -> int:
    if not ADO.is_dir():
        print(f"no ado tree at {ADO}: run crosscheck/stata/setup_stata.do first", file=sys.stderr)
        return 2
    cur = current()
    if "--check" in sys.argv:
        old = {}
        for line in MANIFEST.read_text().splitlines():
            if line and not line.startswith("#"):
                h, rel = line.split(None, 1)
                old[rel] = h
        drift = sorted(set(old) ^ set(cur)) + sorted(k for k in old.keys() & cur.keys() if old[k] != cur[k])
        for k in drift:
            print("  drift  ", k)
        print("ado tree matches the manifest" if not drift else f"{len(drift)} file(s) drifted")
        return 0 if not drift else 1
    lines = ["# SHA-256 of crosscheck/ado/plus (written by crosscheck/ado_manifest.py)"]
    lines += [f"{h}  {rel}" for rel, h in cur.items()]
    MANIFEST.write_text("\n".join(lines) + "\n")
    print(f"wrote {MANIFEST.relative_to(ROOT)} ({len(cur)} files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
