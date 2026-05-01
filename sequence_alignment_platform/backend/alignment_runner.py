import asyncio
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Optional

from .common import manager, update_status

PROJECT_ROOT = Path(__file__).resolve().parents[2]

logger = logging.getLogger("uvicorn")


def _resolve_binary(binary_name: str) -> str:
    """
    Resolve an executable to an absolute path.
    """
    found = shutil.which(binary_name)
    if found:
        return found

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


def _expected_fmindex_path(target_path: Path) -> Path:
    """
    Match the current C++ fmindex naming behavior:
    output file = <target_path.stem>.fmidx in PROJECT_ROOT
    """
    return PROJECT_ROOT / f"{target_path.stem}.fmidx"


async def run_alignment(
    session_dir: Path,
    query_path: Path,
    target_path: Path,
    params: Dict[str, Any],
) -> None:
    session_id = session_dir.name
    queue = manager.get_queue(session_id)

    logger.info(f"\nSession started: {session_id}")

    try:
        await update_status(session_dir, status="running")

        if not query_path.exists():
            await queue.put(f"\n[error] Query file not found: {query_path}\n")
            await update_status(session_dir, status="failed")
            return

        if not target_path.exists():
            await queue.put(f"\n[error] Target file not found: {target_path}\n")
            await update_status(session_dir, status="failed")
            return

        requested_seed = params.get("use_seeded_alignment", False)
        use_seed = bool(requested_seed)

        fmindex_path: Optional[Path] = None

        if use_seed:
            fmindex_bin = _resolve_binary("fmindex")
            fmindex_path = _expected_fmindex_path(target_path)

            # Remove stale index so we do not accidentally reuse an old file
            if fmindex_path.exists():
                try:
                    fmindex_path.unlink()
                except OSError as e:
                    await queue.put(
                        f"\n[error] Could not remove stale FM-index {fmindex_path}: {e}\n"
                    )
                    await update_status(session_dir, status="failed")
                    return

            fmindex_cmd = [
                fmindex_bin,
                str(target_path),
                "-s",
                "$",
            ]

            rc = await _stream_process(
                fmindex_cmd,
                queue,
                session_id,
                cwd=PROJECT_ROOT,
                step_name="fmindex_build",
            )

            if rc != 0:
                await queue.put(
                    f"\n[error] FM-Index generation failed with exit code {rc}\n"
                )
                await update_status(session_dir, status="failed")
                return

            if not fmindex_path.exists():
                await queue.put(
                    f"\n[error] Expected FM-index not found after build: {fmindex_path}\n"
                )
                await update_status(session_dir, status="failed")
                return

            await queue.put(f"[info] FM-index ready: {fmindex_path}\n")

        binary_name = "seed_aligner" if use_seed else "aligner"
        binary_path = _resolve_binary(binary_name)

        choice_map = {
            "global": "1",
            "local": "2",
            "lcs": "3",
            "all": "4",
        }
        cpp_choice = choice_map.get(params.get("align_method", "global"), "1")

        seq_type = params.get("seq_type", "dna")
        if seq_type not in {"dna", "protein"}:
            await queue.put(f"\n[error] Invalid seq_type: {seq_type}\n")
            await update_status(session_dir, status="failed")
            return

        align_cmd = [
            binary_path,
            "--query",
            str(query_path),
            "--target",
            str(target_path),
            "--outdir",
            str(session_dir),
            "--mode",
            seq_type,
            "--choice",
            cpp_choice,
            # "--verbose",
            "--txt",
            "--binary",
        ]

        if use_seed and fmindex_path is not None:
            align_cmd.extend(["--fmindex", str(fmindex_path)])

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

        if params.get("align_method") == "all":
            analysis_cmd = [
                sys.executable,
                "-u",
                "-m",
                "alignment_tool.cli",
                "full",
                "--results-dir",
                str(session_dir),
                "--outdir",
                str(analysis_outdir),
                "--prefix",
                session_id,
                "--blosum",
                "blosum62",
                "--plot-dpi",
                "200",
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
                await queue.put(
                    f"\n[warning] Downstream analysis failed with exit code {rc}\n"
                )
                await queue.put(
                    "[warning] Alignment completed successfully, but comparative analysis was skipped.\n"
                )
            else:
                await queue.put("\n[info] Analysis complete\n")
        else:
            await queue.put(
                "\n[info] Single-method run detected; skipping comparative analysis\n"
            )
            await queue.put("[info] Alignment complete\n")

        await queue.put("\n[info] Session completed successfully\n")
        await update_status(session_dir, status="completed")

    except Exception as e:
        import traceback

        traceback.print_exc()
        try:
            await queue.put(f"\n[error] Python backend crash: {e}\n")
            await update_status(session_dir, status="failed")
        except Exception:
            pass
