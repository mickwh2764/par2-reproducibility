"""Rerun every Paper G analysis and diff it against the committed results.

Usage:
    python paper-g/verify.py            # rerun everything, compare, report
    python paper-g/verify.py --quick    # skip the two GEO-dependent analyses

The fresh run writes to a temporary directory, so ``results/`` is never
overwritten. Numeric fields are compared with a tolerance because the
permutation p-values depend on the RNG draw even with a fixed seed if NumPy's
generator implementation changes.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

PAPER_DIR = Path(__file__).resolve().parent
CODE_DIR = PAPER_DIR / "code"
COMMITTED = PAPER_DIR / "results"

# Fields that are deterministic given the same inputs, and those that are only
# expected to agree approximately.
TOLERANCE = 1e-9
PERMUTATION_TOLERANCE = 0.01

FRESH = COMMITTED


def run(script: str, args: list[str], results: Path) -> None:
    print(f"\n=== running {script} {' '.join(args)}")
    subprocess.run(
        [sys.executable, str(CODE_DIR / script), *args],
        check=True,
        env={**os.environ, "PAR2G_RESULTS_DIR": str(results)},
    )


def compare_json(name: str, tolerance: float) -> list[str]:
    fresh = json.loads((FRESH / name).read_text())
    old = json.loads((COMMITTED / name).read_text())
    return diff_values(name, old, fresh, tolerance)


def diff_values(path: str, old: object, new: object, tolerance: float) -> list[str]:
    if isinstance(old, dict) and isinstance(new, dict):
        problems = []
        for key in sorted(set(old) | set(new)):
            if key not in old or key not in new:
                problems.append(f"{path}.{key}: present in only one run")
                continue
            problems += diff_values(f"{path}.{key}", old[key], new[key], tolerance)
        return problems
    if isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            return [f"{path}: length {len(old)} vs {len(new)}"]
        problems = []
        for index, (a, b) in enumerate(zip(old, new)):
            problems += diff_values(f"{path}[{index}]", a, b, tolerance)
        return problems
    if isinstance(old, (int, float)) and isinstance(new, (int, float)):
        if abs(float(old) - float(new)) > tolerance:
            return [f"{path}: {old} vs {new}"]
        return []
    return [] if old == new else [f"{path}: {old!r} vs {new!r}"]


def compare_csv(name: str, keys: list[str]) -> list[str]:
    old = pd.read_csv(COMMITTED / name).set_index(keys).sort_index()
    new = pd.read_csv(FRESH / name).set_index(keys).sort_index()
    if list(old.index) != list(new.index):
        return [f"{name}: row keys differ"]

    problems = []
    for column in old.columns:
        if pd.api.types.is_numeric_dtype(old[column]):
            worst = (old[column] - new[column]).abs().max()
            if worst > TOLERANCE:
                problems.append(f"{name}.{column}: max abs difference {worst:.3g}")
        elif not old[column].equals(new[column]):
            problems.append(f"{name}.{column}: values differ")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Only rerun the analyses that need no GEO download.",
    )
    args = parser.parse_args()

    global FRESH
    with tempfile.TemporaryDirectory(prefix="par2g-verify-") as tmp:
        FRESH = Path(tmp)

        run("theory_checks.py", [], FRESH)
        problems = compare_json("theory_checks.json", TOLERANCE)

        if args.quick:
            run("fibonacci_null.py", ["--null-only"], FRESH)
            problems += diff_values(
                "fibonacci_null.json.null",
                json.loads((COMMITTED / "fibonacci_null.json").read_text())["null"],
                json.loads((FRESH / "fibonacci_null.json").read_text())["null"],
                PERMUTATION_TOLERANCE,
            )
        else:
            run("hierarchy.py", [], FRESH)
            problems += compare_csv("category_tests.csv", ["tissue", "category"])
            problems += compare_json("hierarchy_summary.json", TOLERANCE)

            run("fibonacci_null.py", [], FRESH)
            problems += compare_json("fibonacci_null.json", PERMUTATION_TOLERANCE)

            run("permutation.py", [], FRESH)
            problems += compare_json("permutation_tests.json", PERMUTATION_TOLERANCE)

    print("\n" + "=" * 60)
    if problems:
        print(f"{len(problems)} difference(s) against the committed results:")
        for problem in problems:
            print(f"  {problem}")
        return 1
    print("Every rerun analysis matches the committed results.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
