import asyncio
import os
from pathlib import Path
from typing import Dict, Any
from .common import manager, update_status

async def run_alignment(session_dir: Path, query_path: Path, target_path: Path, params: Dict[str, Any]) -> None:
    """Main entrypoint for executing the real C++ alignment and analysis pipeline."""
    print(f"\n🚀 REAL BACKGROUND TASK STARTED FOR: {session_dir.name}")

    try:
        session_id = session_dir.name
        queue = manager.get_queue(session_id)
        await update_status(session_dir, status="running")

        # 1. Determine sequence sizes
        q_size = os.path.getsize(query_path)
        t_size = os.path.getsize(target_path)
        use_seed = (q_size * t_size > 10**7) or (q_size > 10000) or (t_size > 10000)

        # NOTE: If your binaries are in the cpp_build folder, update this path:
        # e.g., binary_name = "/app/cpp_build/seed_aligner" (Docker) or "../cpp_build/aligner" (Local)
        binary_name = "seed_aligner" if use_seed else "aligner"

        # Map the string algorithm to an integer for C++ --choice
        # Adjust these integers based on what your C++ switch(align_choice) statement expects!
        choice_map = {"global": "1", "local": "2", "lcs": "3", "all": "4"}
        cpp_choice = choice_map.get(params["align_method"], "1")

        await queue.put(f"[info] Preparing to run {binary_name} via MPI...\n")

        # 2. Construct the terminal command using named flags
        cmd = [
            "mpirun", "-np", "4",
            binary_name,
            "--query", str(query_path),
            "--target", str(target_path),
            "--outdir", str(session_dir),
            "--mode", params["seq_type"],           # "dna" or "protein"
            "--choice", cpp_choice,                 # mapped integer for the algorithm
            "--gap_open", str(params["gap_open"]),
            "--gap_extend", str(params["gap_extend"]),
            "--verbose",                            # Force C++ to print logs for the WebSocket
            "--binary"                              # Force C++ to output dp_matrix.bin
        ]

        await queue.put(f"[info] Executing: {' '.join(cmd)}\n")

        # 3. Launch the C++ process asynchronously
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.STDOUT
        )

        # 4. Stream real-time logs from C++ directly to the React WebSocket
        while True:
            line = await process.stdout.readline()
            if not line:
                break

            log_message = line.decode('utf-8')
            await queue.put(log_message)
            print(f"[{session_id} LOG]: {log_message.strip()}")

        # 5. Wait for the C++ program to exit
        returncode = await process.wait()

        # 6. Handle Success or Failure
        if returncode == 0:
            await queue.put("[info] Analysis complete\n")
            await update_status(session_dir, status="completed")
            print(f"✅ C++ PROCESS FINISHED SUCCESSFULLY FOR: {session_id}\n")
        else:
            await queue.put(f"[error] C++ process failed with exit code {returncode}\n")
            await update_status(session_dir, status="failed")
            print(f"❌ C++ PROCESS FAILED FOR: {session_id} (Code {returncode})\n")

    except Exception as e:
        print(f"❌ PYTHON CRASH IN BACKGROUND TASK:")
        import traceback
        traceback.print_exc()
        try:
            await queue.put(f"[error] Python backend crash: {str(e)}\n")
            await update_status(session_dir, status="failed")
        except:
            pass