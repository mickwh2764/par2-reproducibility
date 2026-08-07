"""Tests that need no network access and no downloaded data."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

CODE_DIR = Path(__file__).resolve().parent.parent / "code"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
sys.path.insert(0, str(CODE_DIR))

import fibonacci_null  # noqa: E402
import geo  # noqa: E402
import hierarchy  # noqa: E402
import theory_checks  # noqa: E402


def test_categories_are_mutually_exclusive():
    categories = hierarchy.load_categories()
    assert categories["ARNTL"] == "clock"
    assert categories["BMAL1"] == "clock"  # alias resolves to the same category
    assert set(categories.values()) <= set(hierarchy.CATEGORY_ORDER)


def test_probe_collapse_keeps_the_widest_ranging_probe():
    expression = pd.DataFrame(
        [[1.0, 1.0, 1.0, 1.0], [0.0, 4.0, 0.0, 4.0], [2.0, 2.5, 2.0, 2.5]],
        index=["p1", "p2", "p3"],
        columns=["s1", "s2", "s3", "s4"],
    )
    mapping = pd.Series({"p1": "GENE1", "p2": "GENE1", "p3": "GENE2"})
    collapsed = geo.collapse_to_symbols(expression, mapping)

    assert sorted(collapsed.index) == ["GENE1", "GENE2"]
    # p2 has the larger interquartile range, so it represents Gene1.
    assert collapsed.loc["GENE1"].tolist() == [0.0, 4.0, 0.0, 4.0]


def test_zt_hours_parses_sample_titles():
    hours = geo.zt_hours(["Liv_CT18", "Liv_CT20", "Liv_CT22"])
    assert hours.tolist() == [18.0, 20.0, 22.0]


def test_ar2_fit_recovers_known_coefficients():
    rng = np.random.default_rng(0)
    beta1, beta2 = 0.6, 0.3
    series = [0.0, 0.0]
    for _ in range(400):
        series.append(beta1 * series[-1] + beta2 * series[-2] + rng.normal(0, 0.05))

    fits = hierarchy.fit_tissue(
        pd.DataFrame([series], index=["ARNTL"]), hierarchy.load_categories()
    )
    assert fits.loc[0, "category"] == "clock"
    assert fits.loc[0, "phi1"] == pytest.approx(beta1, abs=0.05)
    assert fits.loc[0, "phi2"] == pytest.approx(beta2, abs=0.05)


def test_stationarity_triangle():
    assert theory_checks.is_stable(0.5, 0.4)
    assert not theory_checks.is_stable(1.0, 1.0)  # the Fibonacci point
    assert not theory_checks.is_stable(0.0, 1.0)


def test_fibonacci_null_is_seed_deterministic():
    first = fibonacci_null.null_rates(2000, np.random.default_rng(42))
    second = fibonacci_null.null_rates(2000, np.random.default_rng(42))
    assert first == second
    # The stationarity constraint must not push the null rate anywhere near the
    # >80% figure that motivated Paper G's correction.
    assert first["rate_stable_5pct"] < 0.2


def test_committed_theory_checks_all_pass():
    checks = json.loads((RESULTS_DIR / "theory_checks.json").read_text())["checks"]
    assert len(checks) == 7
    assert all(check["passed"] for check in checks)


def test_committed_hierarchy_summary_is_internally_consistent():
    summary = json.loads((RESULTS_DIR / "hierarchy_summary.json").read_text())
    tissues = summary["tissues"]
    assert len(tissues) == summary["n_tissues"] == 12
    for values in tissues.values():
        gap = values["clock_median"] - values["target_median"]
        assert values["gap_clock_minus_target"] == pytest.approx(gap, abs=1e-4)
        assert values["clock_above_target"] == (gap > 0)
    assert summary["tissues_with_positive_gap"] == 12
