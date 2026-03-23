import asyncio
import os
import sys
import shutil
from pathlib import Path
from typing import Dict, Any, Optional

from .common import manager, update_status


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _resolve_binary(binary_name: str) -> str:
    """
    Resolve aligner/seed_aligner to an absolute executable path.
    Adjust candidate paths if your binaries live elsewhere.
    """
    # 1) If already on PATH, use it
    found = shutil.which(binary_name)
    if found:
        return found

    # 2) Common local build locations under repo root
    candidates = [
        PROJECT_ROOT / binary_name,
        PROJECT_ROOT / "build" / binary_name,
        PROJECT_ROOT / "cpp_build" / binary_name,
        PROJECT_ROOT / "bin" / binary_name,
        ]

    for path in candidates:
        if path.exists() and os.access(path, os.X_OK):
            return str(path)

    searched = "\n".join(str(p) for p in candidates)
    raise FileNotFoundError(
        f"Could not find executable '{binary_name}'. "
        f"Searched PATH and these locations:\n{searched}"
    )


async def _stream_process(
        cmd: list[str],
        queue: asyncio.Queue,
        session_id: str,
        *,
        cwd: Optional[Path] = None,
        env: Optional[dict[str, str]] = None,
        step_name: str = "process",
) -> int:
    await queue.put(f"\n[info] Starting {step_name}\n")
    await queue.put(f"[info] CWD: {str(cwd or Path.cwd())}\n")
    await queue.put(f"[info] CMD: {' '.join(cmd)}\n\n")

    process = await asyncio.create_subprocess_exec(
        *cmd,
        cwd=str(cwd) if cwd else None,
        env=env,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.STDOUT,
    )

    assert process.stdout is not None

    while True:
        line = await process.stdout.readline()
        if not line:
            break
        text = line.decode("utf-8", errors="replace")
        await queue.put(text)
        print(f"[{session_id}][{step_name}] {text.rstrip()}")

    return await process.wait()


async def run_alignment(
        session_dir: Path,
        query_path: Path,
        target_path: Path,
        params: Dict[str, Any],
) -> None:
    session_id = session_dir.name
    queue = manager.get_queue(session_id)

    print(f"\n🚀 REAL BACKGROUND TASK STARTED FOR: {session_id}")

    try:
        await update_status(session_dir, status="running")

        q_size = os.path.getsize(query_path)
        t_size = os.path.getsize(target_path)
        use_seed = (q_size * t_size > 10**7) or (q_size > 10000) or (t_size > 10000)

        binary_name = "seed_aligner" if use_seed else "aligner"
        binary_path = _resolve_binary(binary_name)

        choice_map = {
            "global": "1",
            "local": "2",
            "lcs": "3",
            "all": "4",
        }
        cpp_choice = choice_map.get(params["align_method"], "1")

        align_cmd = [
            binary_path,
            "--query", str(query_path),
            "--target", str(target_path),
            "--outdir", str(session_dir),
            "--mode", params["seq_type"],
            "--choice", cpp_choice,
            "--gap_open", str(params["gap_open"]),
            "--gap_extend", str(params["gap_extend"]),
            "--verbose",
            "--binary",
        ]

        rc = await _stream_process(
            align_cmd,
            queue,
            session_id,
            cwd=PROJECT_ROOT,
            step_name="alignment",
        )

        if rc != 0:
            await queue.put(f"\n[error] Alignment failed with exit code {rc}\n")
            await update_status(session_dir, status="failed")
            return

        await queue.put("\n[info] Alignment phase finished successfully\n")

        analysis_outdir = session_dir / "analysis_out"
        py_env = os.environ.copy()
        py_env["PYTHONUNBUFFERED"] = "1"

        # Only run comparative downstream analysis when "all" was selected.
        if params["align_method"] == "all":
            analysis_cmd = [
                sys.executable, "-u",
                "-m", "alignment_tool.cli",
                "full",
                "--results-dir", str(session_dir),
                "--outdir", str(analysis_outdir),
                "--prefix", session_id,
                "--blosum", "blosum62",
                "--plot-dpi", "200",
            ]

            rc = await _stream_process(
                analysis_cmd,
                queue,
                session_id,
                cwd=PROJECT_ROOT,
                env=py_env,
                step_name="analysis-full",
            )

            if rc != 0:
                await queue.put(f"\n[error] Downstream analysis failed with exit code {rc}\n")
                await update_status(session_dir, status="failed")
                return

            await queue.put("\n[info] Analysis complete\n")
        else:
            await queue.put(
                "\n[info] Single-method run detected; skipping comparative analysis\n"
            )
            await queue.put("[info] Alignment complete\n")

        await update_status(session_dir, status="completed")

    except Exception as e:
        import traceback
        traceback.print_exc()
        try:
            await queue.put(f"\n[error] Python backend crash: {e}\n")
            await update_status(session_dir, status="failed")
        except Exception:
            pass