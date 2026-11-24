"""Python interface to the FastSpar command-line program.

The wrapper accepts a pandas dataframe with taxa in rows and samples in
columns, and returns the inferred correlation and covariance matrices.
FastSpar's file-based interface is used internally; all intermediate files
are placed in a temporary directory and removed before this function returns.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Tuple, Union

import numpy as np
import pandas as pd


Executable = Union[str, os.PathLike]


def _prepare_counts(counts: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Index]:
    """Validate the dataframe and return numeric counts and original taxa labels."""
    if not isinstance(counts, pd.DataFrame):
        raise TypeError("counts must be a pandas DataFrame")
    if counts.shape[0] < 4:
        raise ValueError("FastSpar requires at least four taxa")
    if counts.shape[1] < 1:
        raise ValueError("counts must contain at least one sample")
    if counts.index.has_duplicates:
        raise ValueError("taxa index values must be unique")
    if counts.columns.has_duplicates:
        raise ValueError("sample column names must be unique")
    if pd.Index(counts.index.map(str)).has_duplicates:
        raise ValueError("taxa labels must remain unique when converted to text")
    if pd.Index(counts.columns.map(str)).has_duplicates:
        raise ValueError("sample labels must remain unique when converted to text")
    if any("\t" in str(label) or "\n" in str(label) or "\r" in str(label)
           for label in counts.index):
        raise ValueError("taxa labels cannot contain tabs or newlines")
    if any("\t" in str(label) or "\n" in str(label) or "\r" in str(label)
           for label in counts.columns):
        raise ValueError("sample labels cannot contain tabs or newlines")

    try:
        numeric = counts.apply(pd.to_numeric, errors="raise")
        values = numeric.to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError("counts must contain only numeric values") from error
    if not np.isfinite(values).all():
        raise ValueError("counts must contain only finite values")
    if (values < 0).any():
        raise ValueError("counts must be nonnegative")
    return numeric, counts.index.copy()


def _resolve_executable(executable: Executable) -> str:
    executable = os.fspath(executable)
    resolved = shutil.which(executable)
    if resolved is None:
        candidate = Path(executable).expanduser()
        if candidate.is_file() and os.access(str(candidate), os.X_OK):
            resolved = str(candidate.resolve())
        else:
            raise FileNotFoundError(
                "Could not find executable {!r}; install it or pass its path".format(executable)
            )
    return resolved


def _run(command, stage: str) -> None:
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            "{} failed with exit code {}{}".format(
                stage, result.returncode, ":\n" + detail if detail else ""
            )
        )


def infer_correlations(
    counts: pd.DataFrame,
    *,
    fastspar: Executable = "fastspar",
    iterations: int = 50,
    exclude_iterations: int = 10,
    threshold: float = 0.1,
    threads: int = 1,
    seed: int = 1,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Infer SparCC correlation and covariance matrices with FastSpar.

    Parameters
    ----------
    counts:
        Absolute count table, with taxa as the dataframe index and samples as
        columns. Values must be finite, nonnegative numbers. At least four
        taxa are required by FastSpar.
    fastspar:
        FastSpar executable name or path. Defaults to ``fastspar`` on ``PATH``.
    iterations, exclude_iterations, threshold, threads, seed:
        Parameters passed through to the FastSpar command-line program.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        ``(correlation, covariance)`` with taxa as both row and column labels.

    Notes
    -----
    FastSpar's command-line interface requires files. The wrapper uses a
    temporary directory for its input and output tables and deletes it before
    returning; it does not leave result files in the working directory.
    """
    numeric, taxa = _prepare_counts(counts)

    for name, value in (("iterations", iterations),
                        ("exclude_iterations", exclude_iterations),
                        ("threads", threads), ("seed", seed)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise TypeError("{} must be an integer".format(name))
        minimum = 1 if name in ("iterations", "threads") else 0
        if value < minimum:
            raise ValueError("{} must be {}".format(
                name, "at least 1" if minimum == 1 else "nonnegative"
            ))
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float, np.number)):
        raise TypeError("threshold must be a number")
    if not np.isrealobj(threshold):
        raise ValueError("threshold must be a real number")
    if not np.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("threshold must be between 0 and 1")

    resolved_executable = _resolve_executable(fastspar)

    with tempfile.TemporaryDirectory(prefix="fastspar-") as temp_dir:
        temp_path = Path(temp_dir)
        input_path = temp_path / "counts.tsv"
        correlation_path = temp_path / "correlation.tsv"
        covariance_path = temp_path / "covariance.tsv"

        otu_table = numeric.copy()
        otu_table.index = otu_table.index.map(str)
        otu_table.columns = otu_table.columns.map(str)
        otu_table.index.name = "#OTU ID"
        otu_table.to_csv(input_path, sep="\t", float_format="%.17g")

        command = [
            resolved_executable,
            "--otu_table", str(input_path),
            "--correlation", str(correlation_path),
            "--covariance", str(covariance_path),
            "--iterations", str(iterations),
            "--exclude_iterations", str(exclude_iterations),
            "--threshold", str(threshold),
            "--threads", str(threads),
            "--seed", str(seed),
            "--yes",
        ]
        _run(command, "FastSpar")
        if not correlation_path.is_file() or not covariance_path.is_file():
            raise RuntimeError(
                "FastSpar did not produce both output matrices"
            )

        correlation = pd.read_csv(correlation_path, sep="\t", index_col=0)
        covariance = pd.read_csv(covariance_path, sep="\t", index_col=0)

    correlation.index = taxa
    correlation.columns = taxa
    covariance.index = taxa
    covariance.columns = taxa
    return correlation, covariance


def infer_pvalues(
    counts: pd.DataFrame,
    *,
    permutations: int = 1000,
    exact: bool = True,
    iterations: int = 50,
    bootstrap_iterations: int = 5,
    exclude_iterations: int = 10,
    threshold: float = 0.1,
    threads: int = 1,
    seed: int = 1,
    fastspar: Executable = "fastspar",
    fastspar_bootstrap: Executable = "fastspar_bootstrap",
    fastspar_pvalues: Executable = "fastspar_pvalues",
) -> pd.DataFrame:
    """Run FastSpar, bootstrap correlations, and p-value estimation.

    Parameters
    ----------
    counts:
        Absolute count table with taxa as the dataframe index and samples as
        columns.
    permutations:
        Number of bootstrap permutations to generate and analyse.
    exact:
        Calculate exact p-values (default). Set false for pseudo p-values.
    iterations, bootstrap_iterations, exclude_iterations, threshold, threads,
    seed:
        FastSpar inference settings. ``threads`` is also used as the maximum
        number of concurrent bootstrap FastSpar processes. Each bootstrap
        process uses one thread to avoid oversubscribing the machine.
    fastspar, fastspar_bootstrap, fastspar_pvalues:
        Executable names or paths for the three FastSpar programs.

    Returns
    -------
    pandas.DataFrame
        Pairwise p-value matrix with taxa as row and column labels.

    Notes
    -----
    FastSpar's command-line tools use files. The complete pipeline runs in a
    temporary directory that is deleted before this function returns; no
    intermediate or result files are left in the working directory.
    """
    numeric, taxa = _prepare_counts(counts)

    for name, value in (("permutations", permutations),
                        ("iterations", iterations),
                        ("bootstrap_iterations", bootstrap_iterations),
                        ("exclude_iterations", exclude_iterations),
                        ("threads", threads), ("seed", seed)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise TypeError("{} must be an integer".format(name))
        minimum = 1 if name in ("permutations", "iterations", "bootstrap_iterations", "threads") else 0
        if value < minimum:
            raise ValueError("{} must be {}".format(
                name, "at least 1" if minimum == 1 else "nonnegative"
            ))
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float, np.number)):
        raise TypeError("threshold must be a number")
    if not np.isrealobj(threshold):
        raise ValueError("threshold must be a real number")
    if not np.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError("threshold must be between 0 and 1")
    if not isinstance(exact, bool):
        raise TypeError("exact must be a bool")

    fastspar_bin = _resolve_executable(fastspar)
    bootstrap_bin = _resolve_executable(fastspar_bootstrap)
    pvalues_bin = _resolve_executable(fastspar_pvalues)

    with tempfile.TemporaryDirectory(prefix="fastspar-pvalues-") as temp_dir:
        temp_path = Path(temp_dir)
        input_path = temp_path / "counts.tsv"
        observed_cor = temp_path / "observed_correlation.tsv"
        observed_cov = temp_path / "observed_covariance.tsv"
        bootstrap_prefix = temp_path / "bootstrap"
        correlation_prefix = temp_path / "bootstrap_correlation_"
        pvalues_path = temp_path / "pvalues.tsv"

        otu_table = numeric.copy()
        otu_table.index = otu_table.index.map(str)
        otu_table.columns = otu_table.columns.map(str)
        otu_table.index.name = "#OTU ID"
        otu_table.to_csv(input_path, sep="\t", float_format="%.17g")

        # Infer the observed matrix and generate bootstrap count tables.
        _run([
            fastspar_bin, "--otu_table", str(input_path),
            "--correlation", str(observed_cor), "--covariance", str(observed_cov),
            "--iterations", str(iterations), "--exclude_iterations", str(exclude_iterations),
            "--threshold", str(threshold), "--threads", str(threads), "--seed", str(seed),
            "--yes",
        ], "FastSpar observed correlation")
        if not observed_cor.is_file() or not observed_cov.is_file():
            raise RuntimeError("FastSpar did not produce the observed correlation and covariance matrices")
        _run([
            bootstrap_bin, "--otu_table", str(input_path),
            "--number", str(permutations), "--prefix", str(bootstrap_prefix),
            "--threads", str(threads), "--seed", str(seed),
        ], "FastSpar bootstrap")

        def run_bootstrap(index: int) -> None:
            bootstrap_table = temp_path / "{}_{}.tsv".format(bootstrap_prefix.name, index)
            correlation_path = temp_path / "{}{}.tsv".format(correlation_prefix.name, index)
            covariance_path = temp_path / "bootstrap_covariance_{}.tsv".format(index)
            _run([
                fastspar_bin, "--otu_table", str(bootstrap_table),
                "--correlation", str(correlation_path), "--covariance", str(covariance_path),
                "--iterations", str(bootstrap_iterations),
                "--exclude_iterations", str(exclude_iterations),
                "--threshold", str(threshold), "--threads", "1", "--seed", str(seed + index),
                "--yes",
            ], "FastSpar bootstrap correlation {}".format(index))
            if not correlation_path.is_file():
                raise RuntimeError("FastSpar produced no bootstrap correlation for permutation {}".format(index))

        with ThreadPoolExecutor(max_workers=int(threads)) as pool:
            futures = [pool.submit(run_bootstrap, index) for index in range(permutations)]
            for future in as_completed(futures):
                future.result()

        pvalue_command = [
            pvalues_bin, "--otu_table", str(input_path),
            "--correlation", str(observed_cor), "--prefix", str(correlation_prefix),
            "--permutations", str(permutations), "--outfile", str(pvalues_path),
            "--threads", str(threads),
        ]
        if not exact:
            pvalue_command.append("--pseudo")
        _run(pvalue_command, "FastSpar p-value calculation")
        if not pvalues_path.is_file():
            raise RuntimeError("FastSpar p-value program did not produce a p-value matrix")
        pvalues = pd.read_csv(pvalues_path, sep="\t", index_col=0)

    pvalues.index = taxa
    pvalues.columns = taxa
    return pvalues

