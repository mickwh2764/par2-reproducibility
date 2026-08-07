"""Eigenvalue hierarchy across the twelve GSE54650 tissues.

For every tissue, each gene's 24-point time series is fitted with an AR(2)
model and reduced to the eigenvalue modulus |lambda| (the ``par2-circadian``
package on PyPI does the fitting). Genes are then grouped into functional
categories and each category's median |lambda| is compared against the
genome-wide background with a Mann-Whitney U test, Benjamini-Hochberg
corrected within each tissue.

Writes:
  results/gene_lambdas_<tissue>.csv.gz   per-gene fits, one file per tissue
  results/category_tests.csv             per-tissue, per-category statistics
  results/hierarchy_summary.json         the headline numbers quoted in the paper
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from par2.core import fit_ar2
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import geo
from paths import results_dir

CODE_DIR = Path(__file__).resolve().parent

TISSUE_NAMES = {
    "Adr": "adrenal", "Aor": "aorta", "Bstm": "brainstem", "BFat": "brown_fat",
    "Cer": "cerebellum", "Hrt": "heart", "Hyp": "hypothalamus", "Kid": "kidney",
    "Liv": "liver", "Lun": "lung", "Mus": "muscle", "WFat": "white_fat",
}

# Categories are mutually exclusive; the first match in this order wins.
CATEGORY_ORDER = [
    "clock", "target", "housekeeping", "immune", "metabolic",
    "chromatin", "signaling", "dna_repair", "stem",
]


def load_categories() -> dict[str, str]:
    """Return symbol -> category, applying the alias table and priority order."""
    config = json.loads((CODE_DIR / "gene_categories.json").read_text())
    aliases: dict[str, list[str]] = config.get("aliases", {})

    assignment: dict[str, str] = {}
    for category in CATEGORY_ORDER:
        for symbol in config[category]:
            names = [symbol, *aliases.get(symbol, [])]
            for name in names:
                assignment.setdefault(name.upper(), category)
    return assignment


def fit_tissue(matrix: pd.DataFrame, categories: dict[str, str]) -> pd.DataFrame:
    """Fit AR(2) to every gene in one tissue."""
    records = []
    for symbol, values in zip(matrix.index, matrix.to_numpy()):
        if np.isnan(values).any():
            continue
        try:
            fit = fit_ar2(values)
        except ValueError:
            continue
        records.append(
            {
                "gene": symbol,
                "phi1": fit["phi1"],
                "phi2": fit["phi2"],
                "eigenvalue_modulus": fit["eigenvalue"],
                "r2": fit["r2"],
                "root_type": fit["root_type"].lower(),
                "category": categories.get(symbol, "background"),
            }
        )
    return pd.DataFrame.from_records(records)


def test_categories(fits: pd.DataFrame, tissue: str) -> pd.DataFrame:
    """Compare each category against the genome-wide background distribution."""
    background = fits["eigenvalue_modulus"].to_numpy()
    background_median = float(np.median(background))

    rows = []
    for category in CATEGORY_ORDER:
        subset = fits.loc[fits["category"] == category, "eigenvalue_modulus"].to_numpy()
        if subset.size < 3:
            continue
        _, p_value = mannwhitneyu(subset, background, alternative="two-sided")
        pooled_sd = np.sqrt((np.var(subset, ddof=1) + np.var(background, ddof=1)) / 2)
        cohens_d = (subset.mean() - background.mean()) / pooled_sd if pooled_sd > 0 else 0.0
        complex_pct = 100.0 * np.mean(
            fits.loc[fits["category"] == category, "root_type"].eq("complex")
        )
        rows.append(
            {
                "tissue": tissue,
                "category": category,
                "gene_count": int(subset.size),
                "median_lambda": round(float(np.median(subset)), 4),
                "mean_lambda": round(float(subset.mean()), 4),
                "std_lambda": round(float(subset.std(ddof=1)), 4),
                "background_median": round(background_median, 4),
                "mann_whitney_p": p_value,
                "cohens_d": round(float(cohens_d), 4),
                "pct_complex_roots": round(float(complex_pct), 1),
            }
        )

    table = pd.DataFrame(rows)
    if table.empty:
        return table

    # Benjamini-Hochberg within the tissue.
    order = table["mann_whitney_p"].to_numpy().argsort()
    ranks = np.empty_like(order)
    ranks[order] = np.arange(1, len(order) + 1)
    n = len(table)
    adjusted = table["mann_whitney_p"].to_numpy() * n / ranks
    adjusted = np.minimum.accumulate(adjusted[order][::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    table["p_adjusted"] = out
    return table


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", default="GSE54650")
    parser.add_argument(
        "--per-gene",
        action="store_true",
        help="Also write the per-gene fits for every tissue (12 gzipped CSVs).",
    )
    args = parser.parse_args()

    results_dir().mkdir(exist_ok=True)
    categories = load_categories()
    matrix, tissue_columns = geo.load_symbol_matrix(args.accession)

    all_tests = []
    summary: dict[str, object] = {"accession": args.accession, "tissues": {}}
    for code, columns in tissue_columns.items():
        tissue = TISSUE_NAMES.get(code, code.lower())
        fits = fit_tissue(matrix[columns], categories)
        tests = test_categories(fits, tissue)
        all_tests.append(tests)

        if args.per_gene:
            fits.to_csv(results_dir() / f"gene_lambdas_{tissue}.csv.gz", index=False)

        indexed = tests.set_index("category")["median_lambda"]
        clock = float(indexed.get("clock", np.nan))
        target = float(indexed.get("target", np.nan))
        background = float(tests["background_median"].iloc[0])
        summary["tissues"][tissue] = {
            "clock_median": clock,
            "target_median": target,
            "background_median": background,
            "gap_clock_minus_target": round(clock - target, 4),
            "clock_above_target": bool(clock > target),
            "full_hierarchy_holds": bool(clock > target > background),
            "n_genes": len(fits),
        }
        print(
            f"{tissue:12s} clock={clock:.4f} target={target:.4f} "
            f"background={background:.4f} gap={clock - target:+.4f}"
        )

    table = pd.concat(all_tests, ignore_index=True)
    table.to_csv(results_dir() / "category_tests.csv", index=False)

    tissues = summary["tissues"].values()
    gaps = [t["gap_clock_minus_target"] for t in tissues]
    summary["tissues_with_positive_gap"] = int(sum(g > 0 for g in gaps))
    summary["tissues_with_full_hierarchy"] = int(
        sum(t["full_hierarchy_holds"] for t in tissues)
    )
    summary["n_tissues"] = len(gaps)
    summary["median_gap"] = round(float(np.median(gaps)), 4)
    (results_dir() / "hierarchy_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(
        f"\nClock > target in {summary['tissues_with_positive_gap']}/{summary['n_tissues']} "
        f"tissues; median gap {summary['median_gap']:+.4f}. "
        f"Full clock > target > background in "
        f"{summary['tissues_with_full_hierarchy']}/{summary['n_tissues']}."
    )
    print(f"Wrote {results_dir()}/category_tests.csv and hierarchy_summary.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
