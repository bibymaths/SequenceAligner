import json
import logging
import shutil
import uuid
from datetime import datetime
from pathlib import Path

from fastapi.staticfiles import StaticFiles
from fastapi import BackgroundTasks, FastAPI, File, HTTPException, UploadFile, WebSocket, WebSocketDisconnect, Form
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel
from starlette.middleware.cors import CORSMiddleware

from .common import BASE_DATA_DIR, manager
from .alignment_runner import run_alignment
from .analysis_parser import discover_analysis_outputs, parse_tsv

app = FastAPI(title="Sequence Alignment Platform")

# Allow your local React development server to talk to the backend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

logger = logging.getLogger("uvicorn")

static_dir = Path(__file__).parent / "static"

if static_dir.exists():
    app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

@app.get("/")
async def serve_frontend():
    """Serve the React index.html at the root URL."""
    index_file = static_dir / "index.html"
    if index_file.exists():
        return FileResponse(index_file)
    return {"message": "Backend is running, but React index.html was not found."}


class AlignmentRequest(BaseModel):
    gap_open: float = 10.0
    gap_extend: float = 0.5
    align_method: str = "global" # e.g., global, local, lcs
    seq_type: str = "protein"        # dna or protein


class SessionMetadata(BaseModel):
    session_id: str
    timestamp: str
    query_filename: str
    target_filename: str
    parameters: dict
    status: str


def create_session_dir() -> Path:
    sid = str(uuid.uuid4())
    session_dir = BASE_DATA_DIR / sid
    session_dir.mkdir(parents=True, exist_ok=False)
    return session_dir


def write_metadata(session_dir: Path, query_file: str, target_file: str, params: dict, status: str) -> SessionMetadata:
    metadata = SessionMetadata(
        session_id=session_dir.name,
        timestamp=datetime.utcnow().isoformat(),
        query_filename=query_file,
        target_filename=target_file,
        parameters=params,
        status=status,
    )
    with open(session_dir / "metadata.json", "w", encoding="utf-8") as fh:
        json.dump(metadata.dict(), fh, indent=2)
    return metadata


async def update_status(session_dir: Path, status: str) -> None:
    meta_path = session_dir / "metadata.json"
    if not meta_path.exists():
        return
    with open(meta_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    data["status"] = status
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


@app.post("/align")
async def align(
        background_tasks: BackgroundTasks,
        query: UploadFile = File(...),
        target: UploadFile = File(...),
        gap_open: float = Form(10.0),        # Tell FastAPI to look in the FormData body!
        gap_extend: float = Form(0.5),
        align_method: str = Form("global"),
        seq_type: str = Form("dna")
) -> SessionMetadata:

    if seq_type not in {"dna", "protein"}:
        raise HTTPException(status_code=400, detail="Invalid seq_type. Must be 'dna' or 'protein'")

    session_dir = create_session_dir()

    query_path = session_dir / query.filename
    target_path = session_dir / target.filename
    with open(query_path, "wb") as fq:
        shutil.copyfileobj(query.file, fq)
    with open(target_path, "wb") as ft:
        shutil.copyfileobj(target.file, ft)

    # Save all parameters to pass to the runner
    params = {
        "gap_open": gap_open,
        "gap_extend": gap_extend,
        "align_method": align_method,
        "seq_type": seq_type
    }
    metadata = write_metadata(session_dir, query.filename, target.filename, params, status="queued")

    background_tasks.add_task(run_alignment, session_dir, query_path, target_path, params)

    return metadata


@app.get("/session/{session_id}")
async def get_session_metadata(session_id: str) -> SessionMetadata:
    meta_path = BASE_DATA_DIR / session_id / "metadata.json"
    if not meta_path.exists():
        raise HTTPException(status_code=404, detail="Session not found")
    with open(meta_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    return SessionMetadata(**data)


@app.get("/session/{session_id}/results")
async def list_results(session_id: str) -> JSONResponse:
    session_dir = BASE_DATA_DIR / session_id
    if not session_dir.exists():
        raise HTTPException(status_code=404, detail="Session not found")
    files = [str(path.relative_to(session_dir)) for path in session_dir.glob("**/*") if path.is_file()]
    return JSONResponse({"files": files})


@app.get("/session/{session_id}/analysis")
async def list_analysis(session_id: str) -> JSONResponse:
    session_dir = BASE_DATA_DIR / session_id
    analysis_dir = session_dir / "analysis_out"
    if not analysis_dir.exists():
        raise HTTPException(status_code=404, detail="No analysis outputs found")
    data = discover_analysis_outputs(analysis_dir)
    return JSONResponse(data)


@app.get("/session/{session_id}/analysis/table/{filename}")
async def get_analysis_table(session_id: str, filename: str) -> JSONResponse:
    session_dir = BASE_DATA_DIR / session_id
    file_path = session_dir / "analysis_out" / filename
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="File not found")
    records = parse_tsv(file_path)
    return JSONResponse({"records": records})


@app.get("/session/{session_id}/file/{path:path}")
async def get_file(session_id: str, path: str):
    file_path = BASE_DATA_DIR / session_id / path
    if not file_path.exists() or not file_path.is_file():
        raise HTTPException(status_code=404, detail="File not found")
    return FileResponse(file_path)


@app.websocket("/ws/logs/{session_id}")
async def websocket_endpoint(websocket: WebSocket, session_id: str) -> None:
    await manager.connect(session_id, websocket)
    queue = manager.get_queue(session_id)
    try:
        while True:
            msg = await queue.get()
            await websocket.send_text(msg)
    except WebSocketDisconnect:
        manager.disconnect(session_id, websocket)
        logger.info(f"WebSocket disconnected for session {session_id}")