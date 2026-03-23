import React, { useEffect, useMemo, useRef, useState } from "react";
import axios from "axios";
import Plotly from "plotly.js-dist-min";

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

  const expectedBytes = 8 + rows * cols * 4;
  if (arrayBuffer.byteLength < expectedBytes) {
    throw new Error(
        `Incomplete score matrix file. Expected at least ${expectedBytes} bytes, got ${arrayBuffer.byteLength}.`
    );
  }

  const z = [];
  const text = [];
  let offset = startOffset;

  for (let i = 0; i < rows; i++) {
    const zRow = [];
    const textRow = [];

    for (let j = 0; j < cols; j++) {
      const value = view.getInt32(offset, true);
      zRow.push(value);
      textRow.push(`Row: ${i}<br>Col: ${j}<br>Score: ${value}`);
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
  if (code === 48) return "0";
  if (code === 0) return "∅";
  return String.fromCharCode(code);
}

function pointerCodeToValue(code) {
  if (code === 68) return 3;
  if (code === 85) return 2;
  if (code === 76) return 1;
  if (code === 48 || code === 0) return 0;
  return 0;
}

function parsePointerMatrix(arrayBuffer) {
  const view = new DataView(arrayBuffer);
  const { rows, cols, offset: startOffset } = parseHeader(view);

  const expectedBytes = 8 + rows * cols;
  if (arrayBuffer.byteLength < expectedBytes) {
    throw new Error(
        `Incomplete pointer matrix file. Expected at least ${expectedBytes} bytes, got ${arrayBuffer.byteLength}.`
    );
  }

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
      textRow.push(`Row: ${i}<br>Col: ${j}<br>Pointer: ${label}`);

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

function buildPlotConfig(method, parsed, isFullscreen = false) {
  const commonAxis = {
    zeroline: false,
    showgrid: false,
    ticks: "outside",
    tickfont: { size: isFullscreen ? 12 : 10 },
    titlefont: { size: isFullscreen ? 14 : 12 },
    fixedrange: false,
  };

  const baseLayout = {
    autosize: true,
    margin: isFullscreen
        ? { t: 70, r: 24, b: 70, l: 80 }
        : { t: 56, r: 20, b: 56, l: 68 },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    font: {
      family: "Inter, sans-serif",
      color: "#e5eefc",
      size: isFullscreen ? 13 : 11,
    },
    xaxis: {
      ...commonAxis,
      title: "Target sequence index",
    },
    yaxis: {
      ...commonAxis,
      title: "Query sequence index",
      autorange: "reversed",
      scaleanchor: "x",
      scaleratio: 1,
    },
  };

  if (method === "lcs") {
    return {
      data: [
        {
          z: parsed.z,
          text: parsed.text,
          customdata: parsed.labels,
          type: "heatmap",
          colorscale: [
            [0.0, "#f8fafc"],
            [0.33, "#94a3b8"],
            [0.66, "#38bdf8"],
            [1.0, "#1d4ed8"],
          ],
          zsmooth: false,
          showscale: true,
          colorbar: {
            title: "Pointer",
            tickmode: "array",
            tickvals: [0, 1, 2, 3],
            ticktext: ["∅", "L", "U", "D"],
            thickness: isFullscreen ? 18 : 12,
          },
          hovertemplate: "%{text}<extra></extra>",
        },
      ],
      layout: {
        ...baseLayout,
        title: {
          text: "LCS Traceback Pointer Matrix",
          font: { size: isFullscreen ? 22 : 16, family: "Space Grotesk, Inter, sans-serif" },
        },
        annotations:
            parsed.rows <= 60 && parsed.cols <= 60
                ? parsed.labels.flatMap((row, i) =>
                    row.map((label, j) => ({
                      x: j,
                      y: i,
                      text: label,
                      showarrow: false,
                      font: {
                        size: isFullscreen ? 11 : 9,
                        color: "#111827",
                      },
                    }))
                )
                : [],
      },
    };
  }

  return {
    data: [
      {
        z: parsed.z,
        text: parsed.text,
        type: "heatmap",
        colorscale: "Viridis",
        zsmooth: false,
        showscale: true,
        colorbar: {
          title: "Score",
          thickness: isFullscreen ? 18 : 12,
        },
        hovertemplate: "%{text}<extra></extra>",
      },
    ],
    layout: {
      ...baseLayout,
      title: {
        text: `${method.toUpperCase()} DP Score Matrix`,
        font: { size: isFullscreen ? 22 : 16, family: "Space Grotesk, Inter, sans-serif" },
      },
    },
  };
}

function MatrixVisualizer({ sessionId }) {
  const plotRef = useRef(null);
  const containerRef = useRef(null);
  const [sessionMeta, setSessionMeta] = useState(null);
  const [activeTab, setActiveTab] = useState("global");
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState("");
  const [parsedMatrix, setParsedMatrix] = useState(null);
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [plotHeight, setPlotHeight] = useState(680);

  const alignMethod = sessionMeta?.parameters?.align_method;
  const isAll = alignMethod === "all";
  const effectiveMethod = isAll ? activeTab : alignMethod;

  useEffect(() => {
    if (!sessionId) return;

    let cancelled = false;

    async function fetchMetadata() {
      try {
        const res = await api.get(`/session/${sessionId}`);
        if (cancelled) return;

        const method = res.data?.parameters?.align_method;
        setSessionMeta(res.data);
        setActiveTab(method === "all" ? "global" : method);
      } catch (err) {
        if (!cancelled) {
          setError("Could not load session metadata.");
          setLoading(false);
        }
      }
    }

    fetchMetadata();

    return () => {
      cancelled = true;
    };
  }, [sessionId]);

  useEffect(() => {
    if (!sessionId || !effectiveMethod) return;

    let cancelled = false;

    async function fetchMatrix() {
      try {
        setLoading(true);
        setError("");

        const filename = getMatrixFilename(effectiveMethod);
        const response = await api.get(`/session/${sessionId}/file/${filename}`, {
          responseType: "arraybuffer",
        });

        if (cancelled) return;

        const parsed = parseMatrixByMethod(effectiveMethod, response.data);
        setParsedMatrix(parsed);
      } catch (err) {
        const expectedFile = getMatrixFilename(effectiveMethod);
        if (!cancelled) {
          setError(
              err?.message ||
              `Could not load the matrix for ${effectiveMethod.toUpperCase()}. Expected file: ${expectedFile}`
          );
          setParsedMatrix(null);
          if (plotRef.current) {
            Plotly.purge(plotRef.current);
          }
        }
      } finally {
        if (!cancelled) {
          setLoading(false);
        }
      }
    }

    fetchMatrix();

    return () => {
      cancelled = true;
    };
  }, [sessionId, effectiveMethod]);

  useEffect(() => {
    if (!parsedMatrix || !containerRef.current) return;

    const updateSize = () => {
      const containerWidth = containerRef.current?.clientWidth || 900;
      const viewportCap = isFullscreen ? window.innerHeight - 180 : window.innerHeight * 0.72;
      const idealByAspect = containerWidth * (parsedMatrix.rows / parsedMatrix.cols);
      const squareCap = containerWidth;
      const nextHeight = Math.max(420, Math.min(idealByAspect, viewportCap, squareCap));
      setPlotHeight(nextHeight);
    };

    updateSize();

    const resizeObserver = new ResizeObserver(updateSize);
    resizeObserver.observe(containerRef.current);
    window.addEventListener("resize", updateSize);

    return () => {
      resizeObserver.disconnect();
      window.removeEventListener("resize", updateSize);
    };
  }, [parsedMatrix, isFullscreen]);

  useEffect(() => {
    if (!parsedMatrix || !plotRef.current) return;

    const { data, layout } = buildPlotConfig(effectiveMethod, parsedMatrix, isFullscreen);

    Plotly.react(
        plotRef.current,
        data,
        {
          ...layout,
          height: plotHeight,
        },
        {
          responsive: true,
          displaylogo: false,
          scrollZoom: true,
          doubleClick: "reset+autosize",
          plotGlPixelRatio: 2,
          toImageButtonOptions: {
            format: "svg",
            filename: `${sessionId}_${effectiveMethod}_matrix`,
            height: Math.round(plotHeight * 2),
            width: Math.round((containerRef.current?.clientWidth || 1000) * 2),
            scale: 4,
          },
          modeBarButtonsToAdd: [],
        }
    );
  }, [parsedMatrix, effectiveMethod, plotHeight, isFullscreen, sessionId]);

  useEffect(() => {
    const handleFullscreenChange = () => {
      const fs = document.fullscreenElement === containerRef.current;
      setIsFullscreen(fs);

      setTimeout(() => {
        if (plotRef.current) {
          Plotly.Plots.resize(plotRef.current);
        }
      }, 80);
    };

    document.addEventListener("fullscreenchange", handleFullscreenChange);
    return () => document.removeEventListener("fullscreenchange", handleFullscreenChange);
  }, []);

  useEffect(() => {
    return () => {
      if (plotRef.current) {
        Plotly.purge(plotRef.current);
      }
    };
  }, []);

  const title = useMemo(() => {
    if (!effectiveMethod) return "DP Matrix Visualizer";
    if (effectiveMethod === "lcs") return "LCS Pointer Matrix Visualizer";
    return `${effectiveMethod.toUpperCase()} DP Matrix Visualizer`;
  }, [effectiveMethod]);

  const toggleFullscreen = async () => {
    if (!containerRef.current) return;

    try {
      if (document.fullscreenElement === containerRef.current) {
        await document.exitFullscreen();
      } else {
        await containerRef.current.requestFullscreen();
      }
    } catch (e) {
      console.error("Fullscreen toggle failed:", e);
    }
  };

  if (!sessionId || !sessionMeta) return null;

  return (
      <div className="mt-6 rounded-3xl border border-white/10 bg-white/5 shadow-2xl backdrop-blur-xl overflow-hidden">
        <div className="flex flex-col gap-4 border-b border-white/10 bg-white/5 px-6 py-4 md:flex-row md:items-center md:justify-between">
          <div>
            <h3 className="text-xl font-bold text-white">{title}</h3>
            <p className="mt-1 text-sm text-slate-300">
              Interactive dynamic-programming matrix inspection with zoom and export.
            </p>
          </div>

          <div className="flex items-center gap-3">
            {isAll && <MethodTabs activeTab={activeTab} onChange={setActiveTab} />}

            <button
                type="button"
                onClick={toggleFullscreen}
                className="rounded-xl border border-cyan-300/20 bg-cyan-400/10 px-4 py-2 text-sm font-semibold text-cyan-200 hover:bg-cyan-400/20 transition"
            >
              {isFullscreen ? "Exit Fullscreen" : "Fullscreen"}
            </button>
          </div>
        </div>

        <div ref={containerRef} className="p-5">
          {loading && (
              <div className="rounded-2xl border border-cyan-300/10 bg-cyan-400/5 px-4 py-10 text-center text-cyan-200">
                Downloading and parsing matrix...
              </div>
          )}

          {error && (
              <div className="rounded-2xl border border-red-400/20 bg-red-400/10 px-4 py-6 text-center text-red-200">
                {error}
              </div>
          )}

          {!loading && !error && parsedMatrix && (
              <>
                <div className="mb-4 grid grid-cols-1 gap-3 md:grid-cols-3">
                  <div className="rounded-2xl border border-white/10 bg-white/5 px-4 py-3">
                    <p className="text-xs uppercase tracking-wide text-slate-400">Query Length</p>
                    <p className="mt-1 text-lg font-semibold text-white">{parsedMatrix.rows}</p>
                  </div>
                  <div className="rounded-2xl border border-white/10 bg-white/5 px-4 py-3">
                    <p className="text-xs uppercase tracking-wide text-slate-400">Target Length</p>
                    <p className="mt-1 text-lg font-semibold text-white">{parsedMatrix.cols}</p>
                  </div>
                </div>

                <div
                    className={`rounded-2xl border border-white/10 bg-slate-950/40 p-3 ${
                        isFullscreen ? "h-[calc(100vh-180px)]" : ""
                    }`}
                >
                  <div ref={plotRef} className="w-full rounded-xl" />
                </div>
              </>
          )}
        </div>
      </div>
  );
}

export default MatrixVisualizer;