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

## Related work

Michael Whiteside, independent computational systems researcher
([ORCID 0009-0000-0643-5791](https://orcid.org/0009-0000-0643-5791)).

| | |
| --- | --- |
| Researcher Profile | [par2discovery.com/profile](https://par2discovery.com/profile) |
| Package | [`par2-circadian`](https://pypi.org/project/par2-circadian/) · [source](https://github.com/mickwh2764/par2discovery) |
| Platform | [par2discovery.com](https://par2discovery.com) |
| Method preprint | [AR(2) eigenvalue modulus as a measure of temporal persistence in gene expression](https://doi.org/10.21203/rs.3.rs-9283100/v1) |
| Tissue-specific dynamics | [A phase-gated autoregressive framework](https://doi.org/10.21203/rs.3.rs-9214347/v1) |
| Half-life independence | [Context-dependent expression persistence](https://doi.org/10.21203/rs.3.rs-9385465/v1) |

## Licence

The analysis code in this repository is MIT (see [LICENSE](LICENSE)). The GEO
data is redistributed by NCBI under its own terms and is downloaded at run time
rather than vendored here.

**Running these analyses requires more than the MIT licence.** The AR(2) fitting
is not in this repository: it comes from `par2-circadian`, which is published
under the [PolyForm Noncommercial License 1.0.0](https://github.com/mickwh2764/par2discovery/blob/main/LICENSE),
and the PAR(2) methodology is additionally the subject of a pending UK patent
application. So the MIT grant covers only the code kept here, and the pipeline as
a whole may be run for noncommercial purposes — which PolyForm defines to include
use by educational institutions, charitable organisations, public research
organisations and government bodies. Commercial use requires a separate licence
for the package (contact mickwh@msn.com); the MIT licence on this repository does
not grant it.
