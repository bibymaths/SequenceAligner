"""
Alignment task runner for the sequence alignment platform.

This module encapsulates the logic for deciding which C++ binary to run
based on input sequence sizes, invoking the process under MPI, streaming
its logs to the WebSocket manager, downsampling large DP matrices for UI
visualization, triggering plotting scripts and downstream analysis commands,
and updating session status. In this skeleton implementation, the actual
alignment binaries are not invoked. Instead, we simulate long‑running work
using ``asyncio.sleep`` and generate dummy output files.
"""

import asyncio
import json
import os
import random
from pathlib import Path
from typing import Dict

from .common import manager, update_status


async def run_alignment(session_dir: Path, query_path: Path, target_path: Path, params: Dict[str, float]) -> None:
    """Main entrypoint for executing the alignment and analysis pipeline."""
    session_id = session_dir.name
    queue = manager.get_queue(session_id)

    # Determine sequence sizes
    q_size = os.path.getsize(query_path)
    t_size = os.path.getsize(target_path)
    # Rough heuristic: length ~ file size (FASTA includes headers)
    # We use product threshold 1e7; if either >10k, we switch to seed aligner
    use_seed = (q_size * t_size > 10**7) or (q_size > 10000) or (t_size > 10000)
    binary = "seed_aligner" if use_seed else "aligner"
    await update_status(session_dir, status="running")
    await queue.put(f"[info] Starting {binary} for session {session_id}\n")
    await asyncio.sleep(1)  # simulate startup

    # Simulate streaming logs from the alignment binary
    for i in range(5):
        await asyncio.sleep(0.5)
        await queue.put(f"[aligner] step {i+1}/5 completed\n")

    # Write dummy result files
    results_dir = session_dir / "results"
    results_dir.mkdir(exist_ok=True)
    (results_dir / "global_stats.json").write_text(json.dumps({"score": random.randint(100, 1000), "identity": 0.9, "coverage": 0.95, "time": 1.23}))
    (results_dir / "local_stats.json").write_text(json.dumps({"score": random.randint(50, 500), "identity": 0.85, "coverage": 0.8, "time": 0.8}))
    (results_dir / "dp_matrix.bin").write_bytes(os.urandom(1000))

    await queue.put("[info] Alignment completed. Generating plots...\n")
    await asyncio.sleep(1)
    # Simulate plot output
    summary_png = results_dir / "summary.png"
    summary_png.write_bytes(os.urandom(1024))

    await queue.put("[info] Running downstream analysis...\n")
    await update_status(session_dir, status="analyzing")
    await asyncio.sleep(1)

    # Simulate analysis output structure
    analysis_dir = session_dir / "analysis_out"
    analysis_dir.mkdir(exist_ok=True)
    # Create some dummy TSV and PNG files
    for mode in ["global", "local", "lcs"]:
        base = f"sample_{mode}_alignment_summary.tsv"
        with open(analysis_dir / base, "w", encoding="utf-8") as fh:
            fh.write("field\tvalue\nscore\t123\n")
        (analysis_dir / f"sample_{mode}_dp_heatmap.png").write_bytes(os.urandom(2048))
    # Summary JSON
    (analysis_dir / "sample_summary.json").write_text(json.dumps({"modes": ["global", "local", "lcs"], "notes": "Dummy analysis"}))

    await queue.put("[info] Analysis complete\n")
    await update_status(session_dir, status="completed")