"""Numerical checks of the analytic statements in Paper G.

These need no data: each proposition is a claim about the AR(2) companion
matrix or the stability triangle, so it can be checked directly. They exist so
a reader can confirm the algebra is not merely asserted.

Run:  python code/theory_checks.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paths import results_dir

PHI = (1.0 + np.sqrt(5.0)) / 2.0

# Boman's p-number families (Table 6 of Boman 2025) as convergent ratios q,
# against the PAR(2) empirical stable band quoted in Paper G.
STABLE_BAND = (0.52, 0.72)
P_NUMBER_Q = {"p=0/1": 0.6180339887498949, "p=2": 0.6823278038280193, "p=3": 0.7244919590005157}


def companion(beta1: float, beta2: float) -> np.ndarray:
    """Companion matrix of x_t = beta1 x_{t-1} + beta2 x_{t-2}."""
    return np.array([[beta1, beta2], [1.0, 0.0]])


def spectral_radius(beta1: float, beta2: float) -> float:
    return float(np.max(np.abs(np.linalg.eigvals(companion(beta1, beta2)))))


def is_stable(beta1: float, beta2: float) -> bool:
    """Stationarity triangle: |beta2| < 1, beta1 + beta2 < 1, beta2 - beta1 < 1."""
    return abs(beta2) < 1.0 and (beta1 + beta2) < 1.0 and (beta2 - beta1) < 1.0


def check(name: str, statement: str, passed: bool, detail: str) -> dict[str, object]:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}: {statement}\n         {detail}")
    return {"name": name, "statement": statement, "passed": bool(passed), "detail": detail}


def main() -> int:
    checks = []

    # Proposition 1 — Boman's spatial recursion matrix M is the AR(2) companion
    # matrix at (beta1, beta2) = (1, 1): the spatial and temporal recursions are
    # the same operator.
    boman_M = np.array([[1, 1], [1, 0]])
    checks.append(
        check(
            "Proposition 1",
            "Boman's M equals the AR(2) companion matrix C(1,1)",
            np.array_equal(companion(1.0, 1.0), boman_M),
            f"C(1,1) = {companion(1.0, 1.0).tolist()}, M = {boman_M.tolist()}",
        )
    )

    # Theorem 1 — the Fibonacci point is not stationary, and its spectral radius
    # is exactly phi, so Fibonacci dynamics sit outside the admissible region.
    radius = spectral_radius(1.0, 1.0)
    checks.append(
        check(
            "Theorem 1",
            "(1,1) lies outside the stationarity triangle with |lambda| = phi > 1",
            (not is_stable(1.0, 1.0)) and abs(radius - PHI) < 1e-12 and radius > 1.0,
            f"|lambda| = {radius:.15f}, phi = {PHI:.15f}, stable = {is_stable(1.0, 1.0)}",
        )
    )

    # Corollary 1 — phi is a boundary landmark: scaling (1,1) towards the origin
    # crosses into the stable region at exactly the boundary beta1 + beta2 = 1.
    crossing = max(s for s in np.linspace(0.0, 1.0, 100_001) if is_stable(s, s))
    checks.append(
        check(
            "Corollary 1",
            "The ray (s,s) enters the stability triangle at s = 1/2",
            abs(crossing - 0.5) < 1e-4,
            f"largest stable s on the ray = {crossing:.5f}",
        )
    )

    # Proposition 2 — on the equal-coefficient ray (beta1 = beta2 = c),
    # stationarity requires c < 1/2; c = 1/2 gives roots 1 and -1/2 (a unit-root
    # process). The Fibonacci point c = 1 therefore sits exactly a factor of two
    # beyond the stationarity boundary in this direction.
    roots = np.linalg.eigvals(companion(0.5, 0.5))
    roots = np.sort_complex(roots)
    unit_root = abs(max(roots.real) - 1.0) < 1e-12
    damped_root = abs(min(roots.real) + 0.5) < 1e-12
    checks.append(
        check(
            "Proposition 2",
            "equal-coefficient stationarity ends at c = 1/2 (roots 1 and -1/2); "
            "the Fibonacci point c = 1 is twice as far",
            unit_root and damped_root and not is_stable(0.5, 0.5) and is_stable(0.4999, 0.4999),
            f"roots at c=1/2: {[round(float(r.real), 12) for r in roots]}, "
            f"ratio c_fib/c_boundary = {1.0 / 0.5:.1f}",
        )
    )

    # Proposition S1 — the conservation identities that make phi self-inverse
    # under the two operations used in the twinning argument.
    checks.append(
        check(
            "Proposition S1",
            "phi * (1/phi) = 1 and phi - 1 = 1/phi",
            abs(PHI * (1 / PHI) - 1.0) < 1e-15 and abs((PHI - 1) - 1 / PHI) < 1e-15,
            f"phi - 1 = {PHI - 1:.15f}, 1/phi = {1 / PHI:.15f}",
        )
    )

    # Proposition S2 — the integrated memory of the phi-scaled process sums to
    # phi^2, and phi^2 = phi + 1 closes the hierarchy (Corollary S2).
    memory_sum = sum(PHI ** -k for k in range(0, 100)) * (PHI - 1) * PHI
    checks.append(
        check(
            "Proposition S2 / Corollary S2",
            "geometric memory sum equals phi^2, and phi^2 = phi + 1",
            abs(memory_sum - PHI**2) < 1e-9 and abs(PHI**2 - (PHI + 1)) < 1e-12,
            f"sum = {memory_sum:.12f}, phi^2 = {PHI ** 2:.12f}, phi + 1 = {PHI + 1:.12f}",
        )
    )

    # Proposition S3 — the empirical PAR(2) stable band brackets exactly the
    # p=0/1 and p=2 families of Boman's p-numbers, and excludes p>=3.
    inside = {k: STABLE_BAND[0] <= q <= STABLE_BAND[1] for k, q in P_NUMBER_Q.items()}
    checks.append(
        check(
            "Proposition S3",
            f"stable band {STABLE_BAND} contains p=0/1 and p=2 but not p=3",
            inside["p=0/1"] and inside["p=2"] and not inside["p=3"],
            ", ".join(
                f"{k}: q={q:.4f} {'in' if inside[k] else 'out'}" for k, q in P_NUMBER_Q.items()
            ),
        )
    )

    results_dir().mkdir(parents=True, exist_ok=True)
    failures = [c for c in checks if not c["passed"]]
    (results_dir() / "theory_checks.json").write_text(
        json.dumps({"phi": PHI, "checks": checks, "n_failed": len(failures)}, indent=2) + "\n"
    )
    print(f"\n{len(checks) - len(failures)}/{len(checks)} analytic checks pass")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
