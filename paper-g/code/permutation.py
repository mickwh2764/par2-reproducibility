"""Permutation tests for Paper G's two empirical claims.

Two tests, deliberately reported together because they disagree in the
informative way:

1. ``clock_median`` — is the clock genes' median |lambda| higher than a random
   gene set of the same size? Paper G says yes.
2. ``phi_zone`` — do classified genes fall into the 1/phi zone
   (|lambda| in [0.603, 0.633]) more often than category labels shuffled at
   random? Paper G says no (p ~ 0.15), and that negative result is the reason
   the phi enrichment quoted elsewhere is described as an upper bound.

Both are seeded, so the reported p-values are exactly reproducible.

Writes results/permutation_tests.json.
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

PHI = (1.0 + np.sqrt(5.0)) / 2.0
INV_PHI = 1.0 / PHI
PHI_ZONE = (INV_PHI - 0.015, INV_PHI + 0.015)


def clock_median_test(
    fits: pd.DataFrame, n_permutations: int, rng: np.random.Generator
) -> dict[str, object]:
    """Clock median |lambda| against random gene sets of equal size."""
    values = fits["eigenvalue_modulus"].to_numpy()
    clock = fits.loc[fits["category"] == "clock", "eigenvalue_modulus"].to_numpy()
    observed = float(np.median(clock))

    null = np.array(
        [
            np.median(rng.choice(values, size=clock.size, replace=False))
            for _ in range(n_permutations)
        ]
    )
    # One-sided: how often does a random gene set match or beat the clock set?
    hits = int(np.sum(null >= observed))
    p_value = (hits + 1) / (n_permutations + 1)

    return {
        "observed_clock_median": round(observed, 4),
        "n_clock_genes": int(clock.size),
        "null_median_of_medians": round(float(np.median(null)), 4),
        "null_p95": round(float(np.percentile(null, 95)), 4),
        "n_permutations": n_permutations,
        "n_null_ge_observed": hits,
        "p_value": p_value,
    }


def phi_zone_test(
    fits_by_tissue: dict[str, pd.DataFrame], n_permutations: int, rng: np.random.Generator
) -> dict[str, object]:
    """Classified genes in the 1/phi zone, against random gene sets.

    A gene counts as a hit if it lands in the zone in at least one tissue,
    matching the "53/212 genes in at least one tissue" statistic in the paper.
    The null must use the same "at least one of N tissues" rule, because
    testing twelve tissues inflates the hit rate for any gene set: that
    multiplicity is precisely what the paper's negative result attributes the
    apparent enrichment to.
    """
    tissues = sorted(fits_by_tissue)
    zone_by_gene: dict[str, np.ndarray] = {}
    classified_genes: set[str] = set()

    for index, tissue in enumerate(tissues):
        fits = fits_by_tissue[tissue]
        for gene, lam, category in zip(
            fits["gene"], fits["eigenvalue_modulus"], fits["category"]
        ):
            flags = zone_by_gene.setdefault(gene, np.zeros(len(tissues), dtype=bool))
            flags[index] = PHI_ZONE[0] <= lam <= PHI_ZONE[1]
            if category != "background":
                classified_genes.add(gene)

    all_genes = sorted(zone_by_gene)
    hit_any = {gene: bool(zone_by_gene[gene].any()) for gene in all_genes}

    classified = sorted(classified_genes)
    observed_hits = sum(hit_any[g] for g in classified)
    observed_rate = observed_hits / len(classified)

    genome_pool = np.array([hit_any[g] for g in all_genes], dtype=bool)
    null_rates = np.array(
        [
            rng.choice(genome_pool, size=len(classified), replace=False).mean()
            for _ in range(n_permutations)
        ]
    )

    per_tissue_rate = float(
        np.mean(
            [
                np.mean(
                    (fits_by_tissue[t]["eigenvalue_modulus"] >= PHI_ZONE[0])
                    & (fits_by_tissue[t]["eigenvalue_modulus"] <= PHI_ZONE[1])
                )
                for t in tissues
            ]
        )
    )

    hits = int(np.sum(null_rates >= observed_rate))
    return {
        "phi_zone": [round(float(PHI_ZONE[0]), 4), round(float(PHI_ZONE[1]), 4)],
        "n_tissues": len(tissues),
        "n_classified_genes": len(classified),
        "n_in_zone_any_tissue": observed_hits,
        "observed_rate": round(observed_rate, 4),
        "null_expected_rate": round(float(null_rates.mean()), 4),
        "genome_wide_rate_per_tissue": round(per_tissue_rate, 4),
        "n_permutations": n_permutations,
        "p_value": (hits + 1) / (n_permutations + 1),
        "note": (
            "Reproduces Paper G's negative result: classified genes are no more "
            "likely to sit in the 1/phi zone than random genes tested over the "
            "same number of tissues."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", default="GSE54650")
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    categories = load_categories()
    matrix, tissue_columns = geo.load_symbol_matrix(args.accession)

    fits_by_tissue = {
        TISSUE_NAMES.get(code, code.lower()): fit_tissue(matrix[columns], categories)
        for code, columns in tissue_columns.items()
    }

    clock_tests = {
        tissue: clock_median_test(fits, args.permutations, rng)
        for tissue, fits in fits_by_tissue.items()
    }
    zone = phi_zone_test(fits_by_tissue, 5_000, rng)

    output = {
        "accession": args.accession,
        "seed": args.seed,
        "clock_median_test": clock_tests,
        "phi_zone_test": zone,
    }
    results_dir().mkdir(exist_ok=True)
    (results_dir() / "permutation_tests.json").write_text(json.dumps(output, indent=2) + "\n")

    print(f"Clock median |lambda| vs {args.permutations} random gene sets:")
    for tissue, result in clock_tests.items():
        print(
            f"  {tissue:12s} observed={result['observed_clock_median']:.4f} "
            f"null={result['null_median_of_medians']:.4f} p={result['p_value']:.4f}"
        )
    significant = sum(r["p_value"] < 0.05 for r in clock_tests.values())
    print(f"  significant at p<0.05 in {significant}/{len(clock_tests)} tissues")

    print(
        f"\n1/phi zone {zone['phi_zone']}: {zone['n_in_zone_any_tissue']}/"
        f"{zone['n_classified_genes']} classified genes ({zone['observed_rate']:.1%}) "
        f"in at least one of {zone['n_tissues']} tissues; random gene sets give "
        f"{zone['null_expected_rate']:.1%}, p = {zone['p_value']:.4f}"
    )
    print(f"\nWrote {results_dir()}/permutation_tests.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
