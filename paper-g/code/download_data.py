"""Download the public GEO series matrices and platform annotation used by Paper G.

Nothing here is derived data: every file comes straight from NCBI GEO and is
left gzipped exactly as served, so the SHA-256 recorded in
``data/checksums.txt`` identifies the upstream release rather than a local
re-encoding.
"""

import argparse
import csv
import hashlib
import sys
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
RAW_DIR = DATA_DIR / "raw"

# GPL6246 (Affymetrix Mouse Gene 1.0 ST) maps the probe IDs used by GSE54650
# and GSE11923 onto gene symbols.
ANNOTATION_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6246/annot/GPL6246.annot.gz"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, dest: Path) -> None:
    if dest.exists():
        print(f"  have    {dest.name}")
        return
    print(f"  fetching {dest.name} ...", flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as response, tmp.open("wb") as handle:
        while chunk := response.read(1 << 20):
            handle.write(chunk)
    tmp.rename(dest)


def read_accessions() -> list[dict[str, str]]:
    with (DATA_DIR / "accessions.tsv").open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--accession",
        action="append",
        help="Download only this accession (repeatable). Default: all.",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Compare SHA-256 of downloaded files against data/checksums.txt.",
    )
    args = parser.parse_args()

    rows = read_accessions()
    if args.accession:
        wanted = set(args.accession)
        rows = [r for r in rows if r["accession"] in wanted]
        if not rows:
            print(f"No such accession in accessions.tsv: {sorted(wanted)}")
            return 1

    print("Downloading from NCBI GEO:")
    fetch(ANNOTATION_URL, RAW_DIR / "GPL6246.annot.gz")
    for row in rows:
        fetch(row["series_matrix_url"], RAW_DIR / f"{row['accession']}_series_matrix.txt.gz")

    if not args.verify:
        return 0

    expected: dict[str, str] = {}
    checksums = DATA_DIR / "checksums.txt"
    if not checksums.exists():
        print(f"No {checksums.name} to verify against.")
        return 1
    for line in checksums.read_text().splitlines():
        if line and not line.startswith("#"):
            digest, name = line.split()
            expected[name] = digest

    failures = 0
    print("\nVerifying checksums:")
    for name, want in sorted(expected.items()):
        path = RAW_DIR / name
        if not path.exists():
            print(f"  MISSING  {name}")
            failures += 1
            continue
        got = sha256(path)
        status = "ok" if got == want else "MISMATCH"
        if got != want:
            failures += 1
        print(f"  {status:8} {name}")
    if failures:
        print(
            f"\n{failures} file(s) differ from the recorded checksums. GEO does "
            "occasionally re-issue a series matrix; if so the analysis may not "
            "reproduce the committed results exactly."
        )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
