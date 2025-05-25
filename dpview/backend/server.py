from fastapi import FastAPI, Query
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from data_loader import load_matrix
from matrix_utils import downsample

app = FastAPI()
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

MATRIX_PATH = "../data/matrix.bin"
mat, ROWS, COLS = load_matrix(MATRIX_PATH)

@app.get("/meta")
def meta():
    return {"rows": ROWS, "cols": COLS}

@app.get("/tile")
def tile(x: int, y: int, width: int = 100, height: int = 100):
    xmax = min(COLS, x + width)
    ymax = min(ROWS, y + height)
    return JSONResponse(mat[y:ymax, x:xmax].tolist())

@app.get("/downsampled")
def get_downsampled():
    ds = downsample(mat, 1000, 1000)
    return JSONResponse(ds.tolist())