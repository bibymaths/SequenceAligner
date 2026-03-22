"""
CLI entrypoint for the sequence alignment benchmarking framework.

This script reads a YAML configuration file specifying input FASTA files,
number of runs, timeout and thread count. It orchestrates execution of a
set of alignment tools on DNA and protein sequences, measures runtime and
memory usage, parses alignment outputs to compute accuracy metrics and
writes structured results to CSV and JSON files. Logs for each run are
written to the ``logs`` directory.

Usage:

    python benchmark/benchmark.py --config configs/default.yaml

The script must be executed from the repository root to resolve relative
paths correctly.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import platform
import statistics
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional

import psutil  # type: ignore
import yaml  # type: ignore


# Support running this script both as a module (``python -m benchmark.benchmark``)
# and as a plain script (``python benchmark/benchmark.py``).  When executed
# directly the ``__package__`` is empty and relative imports would fail.  We
# therefore add the repository root to ``sys.path`` and import the package
# explicitly.  When imported as part of the ``benchmark`` package the
# relative imports below will work as usual.
if __package__ is None or __package__ == "":
    # Determine repository root (two directories up from this file)
    _repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _repo_root not in sys.path:
        sys.path.insert(0, _repo_root)
    from benchmark import utils  # type: ignore
    from benchmark.runners import (
        run_blast,
        run_bowtie2,
        run_bwa,
        run_mafft,
        run_clustal,
    )  # type: ignore
else:
    # Package context: use relative imports
    from . import utils  # type: ignore
    from .runners import (
        run_blast,
        run_bowtie2,
        run_bwa,
        run_mafft,
        run_clustal,
    )


def setup_logging():
    """Configure root logger to output messages to stdout with timestamps."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )


def load_config(path: str) -> Dict[str, Any]:
    """Load YAML configuration file."""
    with open(path, 'r', encoding='utf-8') as f:
        cfg = yaml.safe_load(f)
    return cfg


def get_environment_info(tool_names: List[str]) -> Dict[str, Any]:
    """Collect information about the runtime environment and tool versions."""
    info: Dict[str, Any] = {}
    info['platform'] = platform.platform()
    info['cpu_cores_physical'] = psutil.cpu_count(logical=False)
    info['cpu_cores_logical'] = psutil.cpu_count(logical=True)
    info['memory_total_gb'] = round(psutil.virtual_memory().total / 1e9, 3)
    info['python_version'] = sys.version
    tool_versions: Dict[str, str] = {}
    for tool in tool_names:
        if not utils.check_executable(tool):
            tool_versions[tool] = 'not found'
            continue
        # Determine how to ask for version: many tools support --version
        try:
            import subprocess
            result = subprocess.run([tool, '--version'], capture_output=True, text=True, timeout=10)
            # Use first line of stdout or stderr
            version_output = result.stdout.strip() or result.stderr.strip()
            version_line = version_output.split('\n')[0] if version_output else ''
            tool_versions[tool] = version_line
        except Exception:
            tool_versions[tool] = 'unknown'
    info['tool_versions'] = tool_versions
    return info


def aggregate_metrics(run_results: List[Dict[str, Any]], metric_key: str) -> Dict[str, Optional[float]]:
    """Aggregate a specific metric across multiple runs.

    Parameters
    ----------
    run_results: list of dict
        List of per-run result dictionaries containing the given metric.
    metric_key: str
        Key within the metrics dict to aggregate (e.g., 'identity',
        'alignment_length').

    Returns
    -------
    Dict[str, Optional[float]]
        Dictionary with aggregated statistics (mean, median, std, min, max).
        If all values are None, returns a dictionary with all entries None.
    """
    values: List[float] = []
    for res in run_results:
        val = res['metrics'].get(metric_key)
        if isinstance(val, (int, float)):
            values.append(float(val))
    if not values:
        return {key: None for key in ['mean', 'median', 'std', 'min', 'max']}
    return {
        'mean': float(statistics.mean(values)),
        'median': float(statistics.median(values)),
        'std': float(statistics.stdev(values)) if len(values) > 1 else 0.0,
        'min': float(min(values)),
        'max': float(max(values)),
    }


def main():
    setup_logging()
    parser = argparse.ArgumentParser(description='Benchmark sequence alignment tools')
    parser.add_argument('--config', required=True, help='Path to YAML configuration file')
    args = parser.parse_args()
    cfg = load_config(args.config)
    # Prepare directories
    results_dir = os.path.join('results')
    logs_dir = os.path.join('logs')
    outputs_dir = os.path.join('outputs')
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    # Extract configuration
    dna_cfg = cfg.get('dna', {})
    protein_cfg = cfg.get('protein', {})
    runs = int(cfg.get('runs', 1))
    timeout = cfg.get('timeout')
    if timeout is not None:
        timeout = int(timeout)
    threads = int(cfg.get('threads', 1))
    # Tools per sequence type
    tool_map = {
        'dna': ['blast', 'bowtie2', 'bwa'],
        'protein': ['blast', 'mafft', 'clustal'],
    }
    # Map names to runner callables
    runner_funcs = {
        'blast': run_blast,
        'bowtie2': run_bowtie2,
        'bwa': run_bwa,
        'mafft': run_mafft,
        'clustal': run_clustal,
    }
    # Collect per-run results
    all_results: Dict[str, Dict[str, List[Dict[str, Any]]]] = {'dna': {}, 'protein': {}}
    # For each run
    for run_idx in range(1, runs + 1):
        logging.info(f"Starting run {run_idx}/{runs}")
        for seq_type, cfg_section in [('dna', dna_cfg), ('protein', protein_cfg)]:
            query_path = cfg_section.get('query')
            target_path = cfg_section.get('target')
            if not query_path or not target_path:
                logging.warning(f"No {seq_type} query/target specified; skipping {seq_type} runs")
                continue
            # Ensure result container exists
            if seq_type not in all_results:
                all_results[seq_type] = {}
            for tool_name in tool_map[seq_type]:
                runner = runner_funcs[tool_name]
                # Prepare directories for this tool
                work_dir = os.path.join(outputs_dir, tool_name)
                os.makedirs(work_dir, exist_ok=True)
                log_file = os.path.join(logs_dir, f"{tool_name}_{seq_type}_run{run_idx}.log")
                res = runner(
                    query_path=query_path,
                    target_path=target_path,
                    sequence_type=seq_type,
                    threads=threads,
                    timeout=timeout,
                    work_dir=work_dir,
                    log_path=log_file,
                )
                if res is None:
                    logging.info(f"{tool_name} not executed for {seq_type}")
                    continue
                # Append results
                all_results.setdefault(seq_type, {}).setdefault(tool_name, []).append(res)
        logging.info(f"Finished run {run_idx}/{runs}")
    # Aggregate and write outputs
    # Prepare CSV files
    runtime_lines = ['sequence_type,tool,mean,median,std,min,max']
    memory_lines = ['sequence_type,tool,mean,median,std,min,max']
    accuracy_lines = ['sequence_type,tool,metric,mean,median,std,min,max']
    # metrics keys to aggregate from metrics dict
    accuracy_metrics = ['identity', 'alignment_length', 'mismatches', 'gap_count', 'query_coverage', 'subject_coverage', 'target_coverage']
    # Flatten results into JSON structure
    full_json: Dict[str, Any] = {'runs': runs, 'results': all_results}
    for seq_type, tool_dict in all_results.items():
        for tool_name, run_results in tool_dict.items():
            # Runtime and memory summarisation across runs
            runtimes = [r['runtime'] for r in run_results if isinstance(r.get('runtime'), (int, float))]
            mems = [r['memory'] for r in run_results if isinstance(r.get('memory'), (int, float))]
            rt_stats = utils.aggregate_numbers(runtimes) if runtimes else {k: None for k in ['mean', 'median', 'std', 'min', 'max']}
            mem_stats = utils.aggregate_numbers(mems) if mems else {k: None for k in ['mean', 'median', 'std', 'min', 'max']}
            runtime_lines.append(
                f"{seq_type},{tool_name},{rt_stats['mean']},{rt_stats['median']},{rt_stats['std']},{rt_stats['min']},{rt_stats['max']}"
            )
            memory_lines.append(
                f"{seq_type},{tool_name},{mem_stats['mean']},{mem_stats['median']},{mem_stats['std']},{mem_stats['min']},{mem_stats['max']}"
            )
            # Accuracy metrics
            for metric in accuracy_metrics:
                stats = aggregate_metrics(run_results, metric)
                accuracy_lines.append(
                    f"{seq_type},{tool_name},{metric},{stats['mean']},{stats['median']},{stats['std']},{stats['min']},{stats['max']}"
                )
    # Write CSV files
    runtime_path = os.path.join(results_dir, 'runtime.csv')
    memory_path = os.path.join(results_dir, 'memory.csv')
    accuracy_path = os.path.join(results_dir, 'accuracy.csv')
    # Always end files with a newline for POSIX compatibility
    with open(runtime_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(runtime_lines) + '\n')
    with open(memory_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(memory_lines) + '\n')
    with open(accuracy_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(accuracy_lines) + '\n')
    # Write full results JSON
    full_results_path = os.path.join(results_dir, 'full_results.json')
    with open(full_results_path, 'w', encoding='utf-8') as jf:
        json.dump(full_json, jf, indent=2)
    # Write environment info
    tool_names = ['blastn', 'blastp', 'bowtie2', 'bowtie2-build', 'bwa', 'mafft', 'clustalo']
    env_info = get_environment_info(tool_names)
    env_path = os.path.join(results_dir, 'environment.json')
    with open(env_path, 'w', encoding='utf-8') as ef:
        json.dump(env_info, ef, indent=2)
    logging.info(f"Benchmarking complete. Results written to {results_dir}")


if __name__ == '__main__':
    main()
