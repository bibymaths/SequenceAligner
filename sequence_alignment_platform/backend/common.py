# backend/common.py
import json
import asyncio
from pathlib import Path
from typing import List
from fastapi import WebSocket

BASE_DATA_DIR = Path("data/sessions")
BASE_DATA_DIR.mkdir(parents=True, exist_ok=True)


async def update_status(session_dir: Path, status: str) -> None:
    """Update the stored status in metadata.json for the session."""
    meta_path = session_dir / "metadata.json"
    if not meta_path.exists():
        return
    with open(meta_path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    data["status"] = status
    with open(meta_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


class ConnectionManager:
    def __init__(self) -> None:
        self.active_connections: dict[str, List[WebSocket]] = {}
        self.log_queues: dict[str, asyncio.Queue] = {}

    async def connect(self, session_id: str, websocket: WebSocket) -> None:
        await websocket.accept()
        self.active_connections.setdefault(session_id, []).append(websocket)
        self.log_queues.setdefault(session_id, asyncio.Queue())

    def disconnect(self, session_id: str, websocket: WebSocket) -> None:
        if session_id in self.active_connections:
            self.active_connections[session_id].remove(websocket)

    def get_queue(self, session_id: str) -> asyncio.Queue:
        return self.log_queues.setdefault(session_id, asyncio.Queue())


manager = ConnectionManager()
