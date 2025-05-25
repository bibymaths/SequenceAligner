## DP Matrix Viewer

A fast, interactive viewer for large dynamic programming (DP) matrices using Python (FastAPI) and JavaScript (PixiJS).

---

### Directory Structure
```
dpviewer/
├── backend/
│   ├── server.py         # FastAPI tile server
│   ├── data_loader.py    # Load binary matrix (int32 memmap)
│   ├── matrix_utils.py   # Downsampling, PCA, clustering
├── frontend/
│   ├── index.html        # HTML shell
│   ├── app.js            # PixiJS tile viewer
│   └── style.css         # Styling
├── data/

