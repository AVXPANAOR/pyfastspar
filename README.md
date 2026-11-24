# FastSpar Python interface

This repository includes `fastspar_wrapper.py`, a Python module for running FastSpar analyses from pandas dataframes. It provides `infer_correlations` to return correlation and covariance matrices, and `infer_pvalues` to run the bootstrap p-value workflow.

This wrapper uses the original [FastSpar project](https://github.com/scwatts/fastspar). For the method and software, see Watts et al., [“FastSpar: rapid and scalable correlation estimation for compositional data”](https://doi.org/10.1093/bioinformatics/bty734), *Bioinformatics* (2018).

## Requirements and import

Install FastSpar and the Python dependencies `pandas` and `numpy` with conda:

```bash
conda install -c bioconda -c conda-forge fastspar pandas numpy
```

The executables `fastspar`, `fastspar_bootstrap`, and `fastspar_pvalues` must be available on `PATH`, unless executable paths are provided as function arguments.

The wrapper is a module in this repository, not a separately published pip package. From Python, run in the cloned repository directory or add that directory to `PYTHONPATH`:

```python
from fastspar_wrapper import infer_correlations, infer_pvalues
```

## Input dataframe

Both functions accept a pandas dataframe with taxa as rows and samples as columns. The dataframe must have at least four taxa and one sample. Values must be numeric, finite, and nonnegative. Taxa and sample labels must be unique, including after conversion to text, and cannot contain tabs or newlines.

```python
import pandas as pd

counts = pd.DataFrame(
    [[10, 14, 8], [3, 5, 2], [7, 4, 9], [2, 6, 3]],
    index=["taxon_1", "taxon_2", "taxon_3", "taxon_4"],
    columns=["sample_1", "sample_2", "sample_3"],
)
```

## Correlation and covariance

`infer_correlations` returns a pair of dataframes: `(correlation, covariance)`. Both are square matrices with taxa as row and column labels.

```python
correlation, covariance = infer_correlations(counts)
```

Call signature:

```python
infer_correlations(
    counts,
    *,
    fastspar="fastspar",
    iterations=50,
    exclude_iterations=10,
    threshold=0.1,
    threads=1,
    seed=1,
)
```

| Parameter | Meaning and accepted values |
| --- | --- |
| `counts` | Input dataframe described above. |
| `fastspar` | FastSpar executable name or path. Defaults to `fastspar` on `PATH`. |
| `iterations` | SparCC estimation iterations; integer of at least 1. Default: `50`. |
| `exclude_iterations` | Rounds for excluding highly correlated pairs; nonnegative integer. Default: `10`. |
| `threshold` | Pair exclusion correlation threshold; real number from `0` to `1`, inclusive. Default: `0.1`. |
| `threads` | Number of FastSpar threads; positive integer supported by the build and machine. Default: `1`. |
| `seed` | Nonnegative integer random seed. Default: `1`. |

## Bootstrap p-values

`infer_pvalues` calculates correlations for the input table and bootstrap tables, then returns a square p-value dataframe with taxa as row and column labels. Set `exact=False` to request pseudo p-values.

```python
pvalues = infer_pvalues(counts, permutations=1000, threads=8)
```

Call signature:

```python
infer_pvalues(
    counts,
    *,
    permutations=1000,
    exact=True,
    iterations=50,
    bootstrap_iterations=5,
    exclude_iterations=10,
    threshold=0.1,
    threads=1,
    seed=1,
    fastspar="fastspar",
    fastspar_bootstrap="fastspar_bootstrap",
    fastspar_pvalues="fastspar_pvalues",
)
```

| Parameter | Meaning and accepted values |
| --- | --- |
| `counts` | Input dataframe described above. |
| `permutations` | Number of bootstrap tables to generate and analyze; positive integer. Default: `1000`. |
| `exact` | Boolean p-value method selector. `True` requests exact p-values (default); `False` requests pseudo p-values. |
| `iterations` | SparCC estimation iterations for the observed input table; positive integer. Default: `50`. |
| `bootstrap_iterations` | SparCC estimation iterations for each bootstrap table; positive integer. Default: `5`. |
| `exclude_iterations` | Rounds for excluding highly correlated pairs in observed and bootstrap analyses; nonnegative integer. Default: `10`. |
| `threshold` | Pair exclusion correlation threshold; real number from `0` to `1`, inclusive. Default: `0.1`. |
| `threads` | Positive integer used for bootstrap generation and p-value calculation. It also sets the maximum number of concurrent bootstrap correlation processes; each process uses one thread. Default: `1`. |
| `seed` | Nonnegative integer random seed. Default: `1`. |
| `fastspar` | FastSpar executable name or path for observed and bootstrap correlation analyses. Default: `fastspar` on `PATH`. |
| `fastspar_bootstrap` | Bootstrap executable name or path. Default: `fastspar_bootstrap` on `PATH`. |
| `fastspar_pvalues` | P-value executable name or path. Default: `fastspar_pvalues` on `PATH`. |

## Temporary files

FastSpar's programs exchange tabular files, so the wrappers use an automatically removed temporary directory for input and intermediate results. They return pandas dataframes and leave no result files in the working directory. The p-value workflow runs its bootstrap correlations concurrently from Python and does not require GNU Parallel.
