# Sequence Alignment Web Platform

## Introduction

Here's a comprehensive blueprint and skeleton implementation for the Sequence Alignment Web Platform. It covers backend
architecture, intelligent execution logic, real-time streaming, visualization components, data handling, and deployment.
A detailed report with citations and design rationale is provided below.

## Highlights

- **FastAPI backend**: Generates unique session IDs, saves uploaded FASTA files, and runs alignment jobs asynchronously.
  The auto-selector examines sequence sizes and switches to the seed-based aligner when full DP would be prohibitive.
- **WebSockets**: Live console streams logs from background tasks without client polling. Status transitions (queued →
  running → analyzing → completed) are persisted in session metadata.
- **Matrix parsing**: Binary DP matrices are parsed via NumPy. Large matrices (>1,000×1,000) are downsampled so the
  frontend can render heatmaps efficiently.
- **Post-analysis integration**: The backend scans analysis output directories, groups files by mode (global, local,
  LCS, etc.), and exposes parsed TSV/CSV tables and linked PNG plots via API endpoints.
- **React frontend skeleton**: Provides components for file upload, virtualized alignment viewing, interactive DP matrix
  visualization using Plotly.js, and a placeholder analysis dashboard. React-window is used to render only visible
  alignment segments, avoiding DOM overload.
- **Dockerfile**: Installs MPI, Gnuplot, ImageMagick and other dependencies; compiles the C++ alignment tools via CMake;
  builds the React frontend; and runs the FastAPI server.

---

## Sequence Alignment Web Platform

### Introduction

Classical sequence alignment algorithms build a dynamic programming (DP) matrix whose size is proportional to the
product of the query and target sequence lengths. For two sequences of length and , the time and space complexity
is [1]. When aligning long reads or scanning a whole genome, this cost becomes prohibitive.

A popular workaround is seed and extend: index both sequences using short k-mers, find candidate matches, and run DP
only in those regions[2]. Modern aligners implement this strategy via compressed full-text indexes such as the FM-index,
which allow sub-linear pattern search and are widely used in bioinformatics[3].

This project designs a web-based platform that wraps C++ alignment tools (standard DP aligner, seed-based aligner with
FM-index, and an FM-index builder). The platform must scale to large genomic inputs, manage alignment sessions, provide
real-time feedback, visualize DP matrices and tracebacks, and expose downstream analysis outputs.

---

## System Architecture

### Backend (FastAPI)

The backend is built with FastAPI, chosen for its asynchronous support and built-in WebSocket capabilities. Each
alignment run creates a unique SessionID (UUID) and corresponding directory under `data/sessions/{id}/` to isolate files
and metadata. Key features include:

- **Session management**: upon upload, the server saves the query and target FASTA files, writes a `metadata.json`
  containing the filenames, parameters, timestamp and initial status (queued), and schedules a background task to
  execute the alignment (see `backend/main.py`).

- **Intelligent execution engine**: the background task (`alignment_runner.py`) inspects the approximate sequence
  lengths (using file size as a proxy) and automatically switches to the seed-based aligner if the product of lengths
  exceeds a threshold or either exceeds 10,000. This heuristic avoids the explosion of naive DP on large inputs[1]. The
  FM-index is leveraged internally to anchor seeds[2][3].

- **Background execution**: alignment and analysis processes run via FastAPI’s BackgroundTasks so that the HTTP request
  returns immediately. Logs are streamed through an in-memory `asyncio.Queue` into a WebSocket endpoint (
  `/ws/logs/{session_id}`), allowing the frontend to display progress in real-time. WebSockets provide persistent,
  bidirectional communication without polling[4].

- **Matrix parsing and downsampling**: DP matrices written as binary files are parsed using NumPy when available (
  `matrix_parser.py`). Matrices larger than 1000×1000 are downsampled by regularly sampling rows and columns to keep the
  data volume manageable for the browser. The full matrix remains on disk for high-resolution plotting.

- **Post-alignment analysis**: after generating the alignment, the backend invokes downstream analysis via a CLI (e.g.,
  `python -m alignment_tool.cli full --results-dir results/ …`). Analysis outputs such as TSV/CSV tables, PNG heatmaps,
  residue support plots and JSON summaries are stored under `analysis_out/` within the session directory. The API
  exposes endpoints to list these files and to parse tabular data into JSON (see `analysis_parser.py`).

- **Results API**: endpoints allow clients to fetch session metadata, enumerate result files, retrieve downsampled DP
  matrices, access parsed analysis tables, and download arbitrary artifacts from within a session.

---

## WebSockets and Task Status

Real-time feedback is critical for user engagement. A `ConnectionManager` handles WebSocket connections per session.
Background tasks push log strings into an `asyncio.Queue`, and the WebSocket coroutine pops messages and forwards them
to all connected clients.

Status values are updated in `metadata.json` and include:

- queued
- running
- plotting
- analyzing
- completed
- failed

---

## Directory Structure

```

sequence_alignment_platform/
├── backend/
│   ├── main.py             # FastAPI app, endpoints and WebSocket management
│   ├── alignment_runner.py # Simulates alignment and analysis workflows
│   ├── matrix_parser.py    # Binary matrix parser and downsampler
│   ├── analysis_parser.py  # Parses analysis outputs
│   └── ...
├── frontend/
│   ├── package.json        # React dependencies
│   └── src/
│       ├── App.js
│       ├── components/
│       │   ├── DropZone.jsx
│       │   ├── AlignmentViewer.jsx
│       │   ├── MatrixVisualizer.jsx
│       │   └── AnalysisDashboard.jsx
│       └── index.js
├── Dockerfile
└── report.md

```

---

## Frontend (React + Tailwind CSS)

The user interface is built in React and styled with Tailwind CSS. Key components:

1. **Drop & Align**: a drag-and-drop zone (`DropZone.jsx`) where users upload query and target FASTA files, choose gap
   penalties and alignment mode. On submission, it posts to `/align` and displays the returned session ID.

2. **Live Console**: a terminal-style panel subscribes to `/ws/logs/{session_id}` via WebSocket and streams
   stdout/stderr. WebSockets enable immediate display of progress updates without polling[4].

3. **Alignment Viewer**: extremely long alignments cannot be rendered directly in the DOM. The viewer uses
   `react-window` to virtualize rows so that only visible segments are mounted, reducing memory footprint and
   re-rendering overhead. This approach is essential for genomic-scale alignments where rendering every character would
   exhaust browser resources[5].

   The viewer will:
    - fetch alignment segments on demand
    - color matches (green), mismatches (cyan), and gaps (red)
    - support panning/zooming

4. **Matrix Visualizer**: uses Plotly.js to display downsampled DP heatmaps and overlay the traceback path. Users can
   zoom and inspect heat intensities interactively. The full matrix is downsampled on the server to avoid transferring
   multi-gigabyte arrays.

5. **Analysis Dashboard**: automatically discovers analysis artifacts for a session and renders:
    - interactive tables (global/local/LCS summaries, residue support, conserved blocks, path metrics)
    - linked figures
    - sorting, filtering, and pagination
    - JSON summary overview cards

---

## Data Handling and I/O

The platform reads and stores multiple types of files:

- **Core alignment outputs**:
    - `global_stats.json`, `local_stats.json`: score, identity, coverage, runtime
    - `*_path.txt`: coordinates for DP traceback overlay
    - `dp_matrix.bin`: DP matrix (parsed via NumPy and downsampled if needed)

- **Analysis artifacts**:
    - TSV/CSV → parsed into structured JSON
    - PNG → served as static assets
    - JSON summaries → global context

Files are grouped by:

- mode (global, local, LCS)
- category (method comparison, residue support, conserved blocks, substitution summary)

All artifacts remain within the session directory to preserve provenance.

---

## Build & Deployment

The provided Dockerfile:

- installs system dependencies (GNUPlot, ImageMagick, jq, OpenMPI)
- compiles C++ tools via CMake
- installs Python and Node dependencies
- builds the React app
- serves backend + frontend via Uvicorn

Runtime characteristics:

- single container exposing port 8000
- asynchronous FastAPI backend
- WebSockets for real-time interaction
- non-blocking handling of large uploads

---

## Conclusion and Future Work

This skeleton demonstrates how to architect a web platform around computationally intensive sequence alignment. By
coupling FastAPI’s asynchronous backend with a virtualized React frontend, users can:

- upload large FASTA files
- monitor progress in real time
- explore DP matrices interactively
- access post-alignment analytics

The intelligent auto-selector ensures that seed-and-extend algorithms are used when full DP is infeasible, leveraging
FM-index anchors to reduce complexity[2][3].

Future enhancements:

- integration of actual C++ binaries
- optimized streaming of binary matrices
- progressive rendering in alignment viewer
- authentication and multi-user support
- extension toward precision-medicine workflows

## Set up and run the sequence‑alignment platform on your local machine using Docker or by running the backend and frontend separately.

---

## Using Docker (recommended for easiest setup)

1. **Install Docker** – make sure Docker Engine is installed and running on your machine.

2. **Unzip the project** – extract the ZIP file you downloaded into a working directory. You should see a folder
   containing `Dockerfile`, `backend/`, `frontend/`, and `report.md`.

3. **Build the image** – open a terminal, change into the root of the extracted directory (the one containing the
   Dockerfile) and run:

   ```bash
   docker build -t sequence-align-platform .
   ```

   This will compile the C++ tools, install Python and Node dependencies, build the React frontend, and bake everything
   into a single image.

4. **Run the container** – once the build finishes, start the server and expose it on port 8000:

   ```bash
   docker run --name seq-align -p 8000:8000 sequence-align-platform
   ```

   The container will start the FastAPI backend and serve the compiled frontend.

5. **Open your browser** – navigate to `http://localhost:8000/`. You should see the “Drop & Align” interface. Drop two
   FASTA files, click **Run Alignment**, and watch the logs stream in real time.

### Stopping the Docker deployment

* When you’re done, return to the terminal running the container and press `Ctrl+C` to stop it. Alternatively, in
  another terminal run:

  ```bash
  docker stop seq-align
  docker rm seq-align
  ```

  These commands stop and remove the container cleanly.

---

## Local development without Docker

If you prefer to run the backend and frontend separately:

1. **Prerequisites** – ensure you have:

    * Python 3.10+ installed (preferably 3.11).
    * Node.js and npm installed.

2. **Unzip the project** into a working directory.

3. **Set up the backend**:

   ```bash
   cd path/to/extracted/project/backend
   python -m venv venv
   source venv/bin/activate  # or `venv\Scripts\activate` on Windows
   pip install fastapi uvicorn[standard] python-multipart numpy
   ```

   Once dependencies are installed, start the server:

   ```bash
   uvicorn main:app --reload --port 8000
   ```

   The backend API will now be available at `http://localhost:8000`.

4. **Set up the frontend** (in a new terminal window):

   ```bash
   cd path/to/extracted/project/frontend
   npm install
   npm start
   ```

   This will start the React development server on port 3000. It proxies API requests to `localhost:8000`.

5. **Use the app** – open `http://localhost:3000` in your browser to interact with the frontend.

### Closing the local setup

* Stop the backend by pressing `Ctrl+C` in the terminal running `uvicorn`.
* Stop the frontend by pressing `Ctrl+C` in the terminal running `npm start`.
* Deactivate and remove the Python virtual environment if desired:

  ```bash
  deactivate
  rm -rf venv
  ```

---

These steps will get the sequence‑alignment platform running on your machine. Let me know if you encounter any issues or
need guidance on customizing or extending the system.

---

## References

[1] Lecture 9: Alignment - Dynamic Programming and Indexing  
https://data-science-sequencing.github.io/Win2018/lectures/lecture9/

[2] Lecture 9: Alignment - Dynamic Programming and Indexing  
https://data-science-sequencing.github.io/Win2018/lectures/lecture9/

[3] FM-index - Wikipedia  
https://en.wikipedia.org/wiki/FM-index

[4] Real-Time Features in FastAPI: WebSockets, Event Streaming, and Push Notifications  
https://python.plainenglish.io/real-time-features-in-fastapi-websockets-event-streaming-and-push-notifications-fec79a0a6812

[5] Rendering large lists with React Virtualized  
https://blog.logrocket.com/rendering-large-lists-react-virtualized/