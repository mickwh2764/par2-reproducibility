"""The stability-constrained null model for phi-proximity of AR(2) coefficients.

Paper G's central methodological point: earlier reports of phi-proximity used a
null that admitted explosive (non-stationary) AR(2) processes, which inflates
the expected rate above 80%. Restricting both the data and the null to
strictly stationary processes drops the null expectation to a few percent.

The null draws (beta1, beta2) uniformly from beta1 in [-2, 2], beta2 in [-1, 1],
keeps only stationary draws, and asks how often |beta1/beta2| lands within a
window of phi. Observed rates use the same rule on the AR(2) fits from
GSE54650.

Writes results/fibonacci_null.json.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import geo
from hierarchy import TISSUE_NAMES, fit_tissue, load_categories
from paths import results_dir
from theory_checks import PHI, is_stable

WINDOWS = (0.05, 0.02)


def ratio(beta1: np.ndarray, beta2: np.ndarray) -> np.ndarray:
    """|beta1 / beta2|, with a zero denominator treated as no ratio."""
    with np.errstate(divide="ignore", invalid="ignore"):
        value = np.abs(beta1 / beta2)
    return np.where(np.abs(beta2) > 1e-10, value, np.nan)


def null_rates(n_draws: int, rng: np.random.Generator) -> dict[str, object]:
    beta1 = rng.uniform(-2.0, 2.0, n_draws)
    beta2 = rng.uniform(-1.0, 1.0, n_draws)

    stable = np.array([is_stable(b1, b2) for b1, b2 in zip(beta1, beta2)])
    ratios = ratio(beta1, beta2)

    result: dict[str, object] = {"n_draws": n_draws, "n_stable": int(stable.sum())}
    for window in WINDOWS:
        near = np.abs(ratios - PHI) < PHI * window
        result[f"rate_unfiltered_{int(window * 100)}pct"] = round(
            float(np.nanmean(near)), 4
        )
        result[f"rate_stable_{int(window * 100)}pct"] = round(
            float(np.nanmean(near[stable])), 4
        )
    return result


def observed_rates(accession: str) -> pd.DataFrame:
    categories = load_categories()
    matrix, tissue_columns = geo.load_symbol_matrix(accession)

    rows = []
    for code, columns in tissue_columns.items():
        tissue = TISSUE_NAMES.get(code, code.lower())
        fits = fit_tissue(matrix[columns], categories)
        classified = fits[fits["category"].isin(["clock", "target"])]
        stable = classified[classified["eigenvalue_modulus"] < 1.0]
        ratios = ratio(stable["phi1"].to_numpy(), stable["phi2"].to_numpy())

        row = {
            "tissue": tissue,
            "n_clock_target": len(classified),
            "n_stable": len(stable),
        }
        for window in WINDOWS:
            near = np.abs(ratios - PHI) < PHI * window
            row[f"hits_{int(window * 100)}pct"] = int(np.nansum(near))
            row[f"rate_{int(window * 100)}pct"] = round(float(np.nanmean(near)), 4)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", default="GSE54650")
    parser.add_argument("--draws", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--null-only",
        action="store_true",
        help="Skip the observed rates (no GEO download needed).",
    )
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    null = null_rates(args.draws, rng)

    print(f"Null model, {args.draws} uniform draws ({null['n_stable']} stationary):")
    for window in WINDOWS:
        tag = int(window * 100)
        print(
            f"  +/-{tag}% of phi: {null[f'rate_unfiltered_{tag}pct']:.1%} unfiltered, "
            f"{null[f'rate_stable_{tag}pct']:.1%} stationary only"
        )

    output: dict[str, object] = {"seed": args.seed, "null": null}

    if not args.null_only:
        observed = observed_rates(args.accession)
        observed.to_csv(results_dir() / "fibonacci_observed.csv", index=False)
        output["observed_accession"] = args.accession
        output["observed"] = observed.to_dict(orient="records")

        print(f"\nObserved rates in {args.accession} (clock and target genes, |lambda| < 1):")
        for row in observed.itertuples():
            print(
                f"  {row.tissue:12s} stable={row.n_stable:3d}/{row.n_clock_target:3d} "
                f"+/-5%={row.rate_5pct:.1%}  +/-2%={row.rate_2pct:.1%}"
            )

    results_dir().mkdir(exist_ok=True)
    (results_dir() / "fibonacci_null.json").write_text(json.dumps(output, indent=2) + "\n")
    print(f"\nWrote {results_dir()}/fibonacci_null.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
