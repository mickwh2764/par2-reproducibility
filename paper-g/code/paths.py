"""Where the scripts read data from and write results to.

``verify.py`` sets ``PAR2G_RESULTS_DIR`` so a fresh run can be diffed against
the committed results without overwriting them.
"""

import os
from pathlib import Path

PAPER_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PAPER_DIR / "data"
RAW_DIR = DATA_DIR / "raw"


def results_dir() -> Path:
    override = os.environ.get("PAR2G_RESULTS_DIR")
    return Path(override).resolve() if override else PAPER_DIR / "results"
