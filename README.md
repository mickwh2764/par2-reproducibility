# PAR(2) reproducibility

Analysis code, dataset accessions and expected results for the published PAR(2)
papers. One directory per paper, each self-contained and runnable from a clean
checkout with public data only.

The AR(2) fitting itself is not reimplemented here — it comes from the
[`par2-circadian`](https://pypi.org/project/par2-circadian/) package on PyPI
([source](https://github.com/mickwh2764/par2discovery)), so the published
package is the citable engine and this repository is the analysis around it.

| Bundle | Paper | Status |
| --- | --- | --- |
| [`paper-g/`](paper-g/) | *A Time-Domain Analogue to Fibonacci Structure via Phase-Gated AR(2) Dynamics*, The Fibonacci Quarterly | accepted, in press |

Further bundles will be added as each paper is accepted.

## Quick start

```bash
pip install -r requirements.txt
python paper-g/code/download_data.py --verify   # ~27 MB from NCBI GEO
python paper-g/code/hierarchy.py --per-gene
python paper-g/code/fibonacci_null.py
python paper-g/code/permutation.py
python paper-g/code/theory_checks.py
```

`python paper-g/verify.py` reruns all of the above into a temporary directory
and diffs it against the committed results, so "does this still reproduce?" is
one command. `python paper-g/verify.py --quick` runs only the analyses that
need no download.

Tested on Python 3.10 with `par2-circadian` 1.1.5, NumPy 2.2, pandas 2.3,
SciPy 1.15.

## Licence

Code is MIT (see [LICENSE](LICENSE)). The GEO data is redistributed by NCBI
under its own terms and is downloaded at run time rather than vendored here.
