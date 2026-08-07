"""Parse GEO series matrices into per-tissue expression matrices.

The two datasets Paper G uses are Affymetrix Mouse Gene 1.0 ST (GPL6246), whose
series matrices are probe-level. Probes are collapsed to gene symbols by taking
the probe with the largest interquartile range, which is the usual choice for
time-series work: it keeps the most dynamic probe rather than averaging a
responsive probe with a flat one.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

RAW_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"


@dataclass(frozen=True)
class Series:
    """A parsed GEO series matrix, probe-level."""

    accession: str
    expression: pd.DataFrame  # probes x samples
    sample_titles: list[str]

    def tissues(self) -> dict[str, list[str]]:
        """Group sample columns by tissue prefix (``Adr_CT18`` -> ``Adr``).

        Sample order within a tissue is the order GEO lists them, which for
        both datasets is ascending circadian time.
        """
        groups: dict[str, list[str]] = {}
        for column, title in zip(self.expression.columns, self.sample_titles):
            match = re.match(r"^([A-Za-z]+)[_-]?CT\d+", title)
            tissue = match.group(1) if match else "all"
            groups.setdefault(tissue, []).append(column)
        return groups


def read_series_matrix(accession: str) -> Series:
    path = RAW_DIR / f"{accession}_series_matrix.txt.gz"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `python code/download_data.py` first."
        )

    sample_titles: list[str] = []
    table_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                sample_titles = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                table_lines.append(line)

    if not table_lines:
        raise ValueError(f"No expression table found in {path}")

    from io import StringIO

    frame = pd.read_csv(StringIO("".join(table_lines)), sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame = frame.apply(pd.to_numeric, errors="coerce")
    return Series(accession=accession, expression=frame, sample_titles=sample_titles)


def probe_to_symbol() -> pd.Series:
    """Map GPL6246 probe IDs to gene symbols (upper-cased)."""
    path = RAW_DIR / "GPL6246.annot.gz"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run `python code/download_data.py` first."
        )

    rows: list[tuple[str, str]] = []
    in_table = False
    header: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            fields = line.rstrip("\n").split("\t")
            if not header:
                header = fields
                continue
            record = dict(zip(header, fields))
            symbol = record.get("Gene symbol", "").strip()
            if symbol and "///" not in symbol:
                rows.append((record["ID"], symbol.upper()))

    mapping = pd.DataFrame(rows, columns=["probe", "symbol"]).set_index("probe")["symbol"]
    return mapping


def collapse_to_symbols(expression: pd.DataFrame, mapping: pd.Series) -> pd.DataFrame:
    """Collapse probes to one row per gene symbol, keeping the widest-ranging probe."""
    symbols = mapping.reindex(expression.index)
    annotated = expression[symbols.notna()]
    symbols = symbols[symbols.notna()]

    spread = annotated.quantile(0.75, axis=1) - annotated.quantile(0.25, axis=1)
    order = spread.sort_values(ascending=False).index
    keep = symbols.loc[order].drop_duplicates()

    collapsed = annotated.loc[keep.index]
    collapsed.index = keep.values
    return collapsed.sort_index()


def load_symbol_matrix(accession: str) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Return a genes-by-samples matrix and the tissue -> columns grouping."""
    series = read_series_matrix(accession)
    collapsed = collapse_to_symbols(series.expression, probe_to_symbol())
    return collapsed, series.tissues()


def zt_hours(sample_titles: list[str]) -> np.ndarray:
    """Extract circadian times (``Adr_CT18`` -> 18.0) for a set of sample titles."""
    times = []
    for title in sample_titles:
        match = re.search(r"CT(\d+)", title)
        times.append(float(match.group(1)) if match else np.nan)
    return np.asarray(times)
