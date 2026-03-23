import React, { useEffect, useMemo, useRef, useState, useCallback } from "react";
import axios from "axios";

const METHODS = ["global", "local", "lcs"];

const api = axios.create({
  baseURL: "http://127.0.0.1:8000",
});

function getMatrixFilename(method) {
  if (method === "lcs") return "lcs_traceback_pointers.bin";
  return `${method}_dp_matrix.bin`;
}

function parseHeader(dataView) {
  if (dataView.byteLength < 8) {
    throw new Error("Binary file is too small to contain matrix dimensions.");
  }
  const rows = dataView.getInt32(0, true);
  const cols = dataView.getInt32(4, true);
  if (rows <= 0 || cols <= 0) {
    throw new Error(`Invalid matrix dimensions: rows=${rows}, cols=${cols}`);
  }
  return { rows, cols, offset: 8 };
}

function parseScoreMatrix(arrayBuffer) {
  const view = new DataView(arrayBuffer);
  const { rows, cols, offset: startOffset } = parseHeader(view);
  const z = [];
  const text = [];
  let offset = startOffset;

  for (let i = 0; i < rows; i++) {
    const zRow = [];
    const textRow = [];
    for (let j = 0; j < cols; j++) {
      const value = view.getInt32(offset, true);
      zRow.push(value);
      textRow.push(`Query: ${i} | Target: ${j} | Score: ${value}`);
      offset += 4;
    }
    z.push(zRow);
    text.push(textRow);
  }
  return { rows, cols, z, text };
}

function pointerCodeToLabel(code) {
  if (code === 68) return "D";
  if (code === 85) return "U";
  if (code === 76) return "L";
  if (code === 48 || code === 0) return "∅";
  return String.fromCharCode(code);
}

function pointerCodeToValue(code) {
  if (code === 68) return 3;
  if (code === 85) return 2;
  if (code === 76) return 1;
  return 0;
}

function parsePointerMatrix(arrayBuffer) {
  const view = new DataView(arrayBuffer);
  const { rows, cols, offset: startOffset } = parseHeader(view);
  const z = [];
  const text = [];
  const labels = [];
  let offset = startOffset;

  for (let i = 0; i < rows; i++) {
    const zRow = [];
    const textRow = [];
    const labelRow = [];
    for (let j = 0; j < cols; j++) {
      const code = view.getUint8(offset);
      const label = pointerCodeToLabel(code);
      const value = pointerCodeToValue(code);
      zRow.push(value);
      labelRow.push(label);
      textRow.push(`Query: ${i} | Target: ${j} | Pointer: ${label}`);
      offset += 1;
    }
    z.push(zRow);
    labels.push(labelRow);
    text.push(textRow);
  }
  return { rows, cols, z, text, labels };
}

function parseMatrixByMethod(method, arrayBuffer) {
  if (method === "lcs") return parsePointerMatrix(arrayBuffer);
  return parseScoreMatrix(arrayBuffer);
}

function MethodTabs({ activeTab, onChange }) {
  return (
      <div className="inline-flex rounded-2xl border border-white/10 bg-white/5 p-1">
        {METHODS.map((method) => (
            <button
                key={method}
                type="button"
                onClick={() => onChange(method)}
                className={`px-4 py-2 text-sm font-semibold rounded-xl transition ${
                    activeTab === method
                        ? "bg-cyan-400/15 text-cyan-200 border border-cyan-300/20"
                        : "text-slate-300 hover:text-white"
                }`}
            >
              {method.toUpperCase()}
            </button>
        ))}
      </div>
  );
}

function getDotColor(val, method, maxAbs) {
  if (method === "lcs") {
    if (val === 1) return "#ec4899";
    if (val === 2) return "#a855f7";
    if (val === 3) return "#fde047";
    return "rgba(255,255,255,0.02)";
  }

  const normalized = Math.min(1, Math.abs(val) / maxAbs);

  if (val > 0) {
    const hue = 180 - (normalized * 120);
    const lightness = 45 + (normalized * 20);
    return `hsl(${hue}, 100%, ${lightness}%)`;
  } else if (val < 0) {
    const hue = 320 + (normalized * 40);
    return `hsla(${hue}, 80%, 55%, 0.4)`;
  }

  return "rgba(255,255,255,0.02)";
}

export default function MatrixVisualizer({ sessionId }) {
  const containerRef = useRef(null);
  const canvasRef = useRef(null);
  const [sessionMeta, setSessionMeta] = useState(null);
  const [activeTab, setActiveTab] = useState("global");
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");
  const [parsedMatrix, setParsedMatrix] = useState(null);
  const [parsedPath, setParsedPath] = useState(null); // NEW: State to hold the traceback path
  const [isFullscreen, setIsFullscreen] = useState(false);

  const [hoverInfo, setHoverInfo] = useState({ show: false, x: 0, y: 0, text: "" });

  const alignMethod = sessionMeta?.parameters?.align_method;
  const isAll = alignMethod === "all";
  const effectiveMethod = isAll ? activeTab : alignMethod;

  // 1. Fetch Meta
  useEffect(() => {
    if (!sessionId) return;
    let cancelled = false;
    async function fetchMetadata() {
      try {
        const res = await api.get(`/session/${sessionId}`);
        if (!cancelled) {
          setSessionMeta(res.data);
          setActiveTab(res.data?.parameters?.align_method === "all" ? "global" : res.data?.parameters?.align_method);
        }
      } catch (err) {
        if (!cancelled) { setError("Could not load session metadata."); setLoading(false); }
      }
    }
    fetchMetadata();
    return () => { cancelled = true; };
  }, [sessionId]);

  // 2. Fetch Matrix & Path Data
  useEffect(() => {
    if (!sessionId || !effectiveMethod) return;
    let cancelled = false;

    async function fetchData() {
      try {
        setLoading(true); setError("");

        // Fetch the binary matrix
        const matrixReq = api.get(`/session/${sessionId}/file/${getMatrixFilename(effectiveMethod)}`, { responseType: "arraybuffer" });

        // Fetch the traceback path text file
        const pathReq = api.get(`/session/${sessionId}/file/${effectiveMethod}_path.txt`, { responseType: "text" })
            .catch(() => ({ data: null })); // Catch gracefully if path file doesn't exist

        const [matrixRes, pathRes] = await Promise.all([matrixReq, pathReq]);

        if (cancelled) return;

        // Parse Matrix
        setParsedMatrix(parseMatrixByMethod(effectiveMethod, matrixRes.data));

        // Parse Path (if available)
        if (pathRes.data) {
          const pathCoords = pathRes.data
              .trim()
              .split('\n')
              .map(line => {
                const [r, c] = line.trim().split(/\s+/).map(Number);
                return { r, c };
              })
              .filter(pt => !isNaN(pt.r) && !isNaN(pt.c)); // Ignore empty/invalid lines

          setParsedPath(pathCoords);
        } else {
          setParsedPath(null);
        }

      } catch (err) {
        if (!cancelled) setError("Could not load the matrix data.");
      } finally {
        if (!cancelled) setLoading(false);
      }
    }

    fetchData();
    return () => { cancelled = true; };
  }, [sessionId, effectiveMethod]);

  // 3. Render Canvas
  const drawMatrix = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas || !parsedMatrix) return;

    const ctx = canvas.getContext("2d");
    const { rows, cols, z } = parsedMatrix;

    const parent = canvas.parentElement;
    canvas.width = parent.clientWidth;
    canvas.height = isFullscreen ? window.innerHeight - 150 : 600;

    ctx.clearRect(0, 0, canvas.width, canvas.height);

    const cellWidth = canvas.width / cols;
    const cellHeight = canvas.height / rows;

    let maxAbs = 1;
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        maxAbs = Math.max(maxAbs, Math.abs(z[i][j]));
      }
    }

    // A. Draw the score dots
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        const val = z[i][j];

        if (effectiveMethod === "global" && val <= 0) continue;
        if (val === 0) continue;

        const xPos = j * cellWidth + cellWidth / 2;
        const yPos = i * cellHeight + cellHeight / 2;

        const radius = Math.max(1.5, Math.min(cellWidth/2, cellHeight/2) * 0.8 * (Math.abs(val) / maxAbs + 0.3));

        ctx.beginPath();
        ctx.arc(xPos, yPos, radius, 0, 2 * Math.PI);
        ctx.fillStyle = getDotColor(val, effectiveMethod, maxAbs);
        ctx.fill();
      }
    }

    // B. OVERLAY THE ALIGNMENT PATH
    if (parsedPath && parsedPath.length > 0) {
      ctx.beginPath();

      // Styling for the path line (Stark glowing white)
      ctx.strokeStyle = "rgba(255, 255, 255, 0.9)";
      ctx.lineWidth = Math.max(2, Math.min(cellWidth, cellHeight) * 0.25); // Scale thickness with zoom
      ctx.lineCap = "round";
      ctx.lineJoin = "round";

      // Add a cool neon glow effect to the line
      ctx.shadowBlur = 10;
      ctx.shadowColor = "rgba(255, 255, 255, 0.6)";

      // Draw lines connecting the coordinate centers
      for (let k = 0; k < parsedPath.length; k++) {
        const pt = parsedPath[k];
        const px = pt.c * cellWidth + cellWidth / 2;
        const py = pt.r * cellHeight + cellHeight / 2;

        if (k === 0) {
          ctx.moveTo(px, py);
        } else {
          ctx.lineTo(px, py);
        }
      }

      ctx.stroke();
      ctx.shadowBlur = 0; // Reset shadow so it doesn't affect future renders
    }

  }, [parsedMatrix, parsedPath, effectiveMethod, isFullscreen]);

  // Handle Resize
  useEffect(() => {
    drawMatrix();
    window.addEventListener("resize", drawMatrix);
    return () => window.removeEventListener("resize", drawMatrix);
  }, [drawMatrix]);

  // Handle Canvas Hover Tooltip
  const handleMouseMove = (e) => {
    if (!canvasRef.current || !parsedMatrix) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    const cellWidth = rect.width / parsedMatrix.cols;
    const cellHeight = rect.height / parsedMatrix.rows;

    const col = Math.floor(x / cellWidth);
    const row = Math.floor(y / cellHeight);

    if (row >= 0 && row < parsedMatrix.rows && col >= 0 && col < parsedMatrix.cols) {
      const val = parsedMatrix.z[row][col];

      // Let users inspect the path even if the score is 0
      const isPathPoint = parsedPath?.some(pt => pt.r === row && pt.c === col);

      if (!isPathPoint && (val === 0 || (effectiveMethod === "global" && val <= 0))) {
        setHoverInfo({ show: false, x: 0, y: 0, text: "" });
        return;
      }

      setHoverInfo({
        show: true,
        x: e.clientX,
        y: e.clientY,
        text: parsedMatrix.text[row][col] + (isPathPoint ? "<br><span style='color: #38bdf8'>★ On Alignment Path</span>" : "")
      });
    } else {
      setHoverInfo({ show: false, x: 0, y: 0, text: "" });
    }
  };

  const toggleFullscreen = () => {
    if (!containerRef.current) return;
    if (document.fullscreenElement) document.exitFullscreen();
    else containerRef.current.requestFullscreen();
  };

  useEffect(() => {
    const handler = () => setIsFullscreen(!!document.fullscreenElement);
    document.addEventListener("fullscreenchange", handler);
    return () => document.removeEventListener("fullscreenchange", handler);
  }, []);

  if (!sessionId || !sessionMeta) return null;

  return (
      <div
          ref={containerRef}
          className="mt-6 rounded-3xl border border-white/10 bg-slate-950 shadow-2xl overflow-hidden relative"
      >
        <div className="flex flex-col gap-4 border-b border-white/10 bg-slate-900/80 px-6 py-4 md:flex-row md:items-center md:justify-between">
          <div>
            <h3 className="text-xl font-bold text-white">DP Matrix & Path Explorer</h3>
            <p className="mt-1 text-sm text-slate-400">Ultra-fast rendering overlaid with traceback alignment path.</p>
          </div>
          <div className="flex items-center gap-3">
            {isAll && <MethodTabs activeTab={activeTab} onChange={setActiveTab} />}
            <button onClick={toggleFullscreen} className="rounded-xl border border-cyan-300/20 bg-cyan-400/10 px-4 py-2 text-sm text-cyan-200 hover:bg-cyan-400/20 transition">
              {isFullscreen ? "Exit Fullscreen" : "Fullscreen"}
            </button>
          </div>
        </div>

        <div className="p-5 w-full">
          {loading && <div className="text-cyan-200 text-center py-10 font-medium">Parsing high-performance matrix...</div>}
          {error && <div className="text-red-400 text-center py-10 bg-red-950/20 rounded-xl border border-red-500/20">{error}</div>}

          {!loading && !error && parsedMatrix && (
              <div
                  className="w-full relative rounded-xl border border-white/5 overflow-hidden"
                  style={{ backgroundColor: '#0f172a' }}
                  onMouseMove={handleMouseMove}
                  onMouseLeave={() => setHoverInfo({ show: false, x:0, y:0, text:""})}
              >
                <canvas ref={canvasRef} className="w-full block" style={{ cursor: 'crosshair' }} />
              </div>
          )}
        </div>

        {hoverInfo.show && (
            <div
                className="fixed pointer-events-none bg-slate-800 text-white text-xs p-3 rounded-lg shadow-2xl z-50 border border-slate-600 font-mono tracking-tight"
                style={{ left: hoverInfo.x + 15, top: hoverInfo.y + 15 }}
                dangerouslySetInnerHTML={{ __html: hoverInfo.text }}
            />
        )}
      </div>
  );
}