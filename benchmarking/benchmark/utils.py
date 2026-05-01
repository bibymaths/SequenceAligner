"""
Utility functions for the benchmarking framework.

This module centralizes common functionality such as running external commands
with timing and memory measurement, reading FASTA files and computing
statistics. Having a single place for these helpers avoids duplication across
different runners and parsers.

The timing and memory measurement relies on :mod:`psutil` to track peak
resident set size (RSS) across a process and all of its descendants. If
psutil is unavailable or access is restricted, the memory usage will be
reported as ``None``.
"""

from __future__ import annotations

import shutil
import subprocess
import time
import statistics
from typing import Dict, List, Optional, Tuple

import psutil  # type: ignore


def check_executable(cmd: str) -> bool:
    """Return True if ``cmd`` exists on the system PATH.

    Parameters
    ----------
    cmd: str
        Name of the executable to search for.

    Returns
    -------
    bool
        True if the command is found in the PATH, False otherwise.
    """
    return shutil.which(cmd) is not None


def read_fasta_lengths(path: str) -> Dict[str, int]:
    """Read a FASTA file and return the length of each sequence.

    This helper reads a FASTA file in a streaming fashion and returns a
    mapping from sequence identifier (the part of the header following
    the initial '>') to the number of residues.

    Parameters
    ----------
    path: str
        Path to the FASTA file on disk.

    Returns
    -------
    Dict[str, int]
        Dictionary mapping sequence IDs to sequence lengths.
    """
    lengths: Dict[str, int] = {}
    current_id: Optional[str] = None
    current_len = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save previous sequence length
                if current_id is not None:
                    lengths[current_id] = current_len
                # Extract ID up to first whitespace
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
        # Save last sequence
        if current_id is not None:
            lengths[current_id] = current_len
    return lengths


def run_subprocess_with_resource_tracking(
    cmd: List[str],
    timeout: Optional[int] = None,
    capture_output: bool = True,
) -> Tuple[float, Optional[float], int, str, str]:
    """Run an external command measuring wall-clock time and peak memory.

    This function uses :mod:`psutil` to monitor the memory usage of the
    spawned process and all of its children. It reports the peak resident
    memory in megabytes (MB) and the total elapsed wall-clock time in
    seconds. If monitoring fails, the peak memory will be returned as
    ``None``.

    Parameters
    ----------
    cmd: list of str
        The command and its arguments to execute.
    timeout: int, optional
        Maximum number of seconds to allow the process to run. If the
        timeout is reached, the process and its children will be killed.
    capture_output: bool
        Whether to capture stdout and stderr from the command. If false,
        stdout and stderr will be set to empty strings.

    Returns
    -------
    Tuple[float, Optional[float], int, str, str]
        A tuple containing (elapsed_time, peak_memory_mb, exit_code,
        stdout, stderr).
    """
    start_time = time.time()
    # Start the process
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE if capture_output else None,
        stderr=subprocess.PIPE if capture_output else None,
        text=True,
    )
    peak_rss = 0.0
    try:
        ps_proc = psutil.Process(process.pid)
        while True:
            if process.poll() is not None:
                break
            try:
                # Compute memory of process and children
                rss = ps_proc.memory_info().rss
                for child in ps_proc.children(recursive=True):
                    try:
                        rss += child.memory_info().rss
                    except psutil.NoSuchProcess:
                        continue
                if rss > peak_rss:
                    peak_rss = rss
            except psutil.NoSuchProcess:
                # Process terminated between polling and memory collection
                pass
            time.sleep(0.1)
        # Wait for process to finish
        try:
            stdout, stderr = process.communicate(timeout=0)
        except subprocess.TimeoutExpired:
            stdout, stderr = "", ""
    except Exception:
        # Unexpected error, kill the process
        process.kill()
        stdout, stderr = process.communicate()
    end_time = time.time()
    elapsed = end_time - start_time
    if timeout is not None and elapsed > timeout:
        # Kill if over time
        process.kill()
        exit_code = -1
    else:
        exit_code = process.returncode
    # Convert memory to MB
    peak_mb: Optional[float] = peak_rss / (1024 * 1024) if peak_rss > 0 else None
    if not capture_output:
        stdout = ""
        stderr = ""
    return elapsed, peak_mb, exit_code, stdout, stderr


def aggregate_numbers(values: List[float]) -> Dict[str, float]:
    """Aggregate a list of numeric values, computing summary statistics.

    Parameters
    ----------
    values: list of float
        The numbers to summarise.

    Returns
    -------
    Dict[str, float]
        Dictionary containing mean, median, standard deviation (std), minimum
        and maximum of the provided values. If the list is empty, all
        statistics will be ``float('nan')``.
    """
    if not values:
        return {
            "mean": float("nan"),
            "median": float("nan"),
            "std": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
        }
    return {
        "mean": float(statistics.mean(values)),
        "median": float(statistics.median(values)),
        "std": float(statistics.stdev(values)) if len(values) > 1 else 0.0,
        "min": float(min(values)),
        "max": float(max(values)),
    }
