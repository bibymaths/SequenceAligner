import React, { useEffect, useRef, useState } from 'react';
import Plotly from 'plotly.js-dist-min';
import axios from 'axios';

function MatrixVisualizer({ sessionId }) {
  const ref = useRef(null);
  const [sessionMeta, setSessionMeta] = useState(null);
  const [activeTab, setActiveTab] = useState('global'); // Default tab

  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');

  // 1. Fetch Session Metadata to know what algorithm was run
  useEffect(() => {
    if (!sessionId) return;
    const fetchMetadata = async () => {
      try {
        const res = await axios.get(`http://127.0.0.1:8000/session/${sessionId}`);
        const method = res.data.parameters.align_method;
        setSessionMeta(res.data);

        // If they chose 'all', default to 'global' first. Otherwise, use their choice.
        setActiveTab(method === 'all' ? 'global' : method);
      } catch (err) {
        console.error("Failed to fetch metadata for matrix:", err);
      }
    };
    fetchMetadata();
  }, [sessionId]);

  // 2. Fetch and Plot the Matrix based on the Active Tab
  useEffect(() => {
    if (!sessionId || !activeTab) return;

    const fetchAndPlotMatrix = async () => {
      try {
        setLoading(true);
        setError('');

        // Map the active tab to your specific C++ output filenames
        let filename = `${activeTab}_dp_matrix.bin`;
        if (activeTab === 'lcs') {
          filename = 'lcs_traceback_pointers.bin';
        }

        const url = `http://127.0.0.1:8000/session/${sessionId}/file/${filename}`;
        const response = await axios.get(url, { responseType: 'arraybuffer' });

        // 1. Create a basic DataView so we can read different types from the same buffer
        const dataView = new DataView(response.data);

        // 2. The first two values are ALWAYS 4-byte integers (num_rows, num_cols)
        // 'true' indicates Little Endian (standard for x86/C++)
        const numRows = dataView.getInt32(0, true);
        const numCols = dataView.getInt32(4, true);

        const zData = [];
        let offset = 8; // Start after the two 4-byte ints

        if (activeTab === 'lcs') {
          // --- LCS PARSING (1-byte chars) ---
          for (let i = 0; i < numRows; i++) {
            const row = [];
            for (let j = 0; j < numCols; j++) {
              // Read 1 byte as an unsigned integer (ASCII code)
              const charCode = dataView.getUint8(offset);

              // Convert ASCII to a numeric value for the heatmap
              // e.g., mapping specific direction characters to numbers
              // If your C++ uses 'U', 'L', 'D' for Up, Left, Diagonal:
              let val = 0;
              if (charCode === 68) val = 3; // 'D' -> 3
              if (charCode === 85) val = 2; // 'U' -> 2
              if (charCode === 76) val = 1; // 'L' -> 1

              row.push(val || charCode); // fallback to raw ASCII if no match
              offset += 1;
            }
            zData.push(row);
          }
        } else {
          // --- GLOBAL/LOCAL PARSING (4-byte ints) ---
          for (let i = 0; i < numRows; i++) {
            const row = [];
            for (let j = 0; j < numCols; j++) {
              row.push(dataView.getInt32(offset, true));
              offset += 4;
            }
            zData.push(row);
          }
        }

        // Configure Plotly
        const trace = {
          z: zData,
          type: 'heatmap',
          colorscale: activeTab === 'lcs' ? 'Cividis' : 'Viridis',
          colorbar: { title: 'Score/Pointer' }
        };

        const layout = {
          title: `${activeTab.toUpperCase()} Dynamic Programming Heatmap`,
          margin: { t: 40, b: 40, l: 40, r: 40 },
          xaxis: { title: 'Target Sequence Index' },
          yaxis: { title: 'Query Sequence Index', autorange: 'reversed' },
          hovermode: 'closest'
        };

        Plotly.newPlot(ref.current, [trace], layout);

      } catch (err) {
        console.error("Matrix Visualizer Error:", err);
        setError(`Could not load the matrix for ${activeTab.toUpperCase()}. Expected file: ${activeTab === 'lcs' ? 'lcs_traceback_pointers.bin' : `${activeTab}_dp_matrix.bin`}`);
        if (ref.current) Plotly.purge(ref.current); // Clear old plot on error
      } finally {
        setLoading(false);
      }
    };

    fetchAndPlotMatrix();

    return () => {
      if (ref.current) Plotly.purge(ref.current);
    };
  }, [sessionId, activeTab]);

  if (!sessionId || !sessionMeta) return null;

  const isAll = sessionMeta.parameters.align_method === 'all';

  return (
      <div className="bg-white rounded-lg shadow-md border border-gray-200 mt-6 overflow-hidden">
        {/* Header & Tabs */}
        <div className="bg-gray-50 px-6 py-3 border-b border-gray-200 flex flex-col md:flex-row justify-between items-center gap-4">
          <h3 className="text-lg font-bold text-gray-800">DP Matrix Visualizer</h3>

          {/* Render Tabs ONLY if 'All' was selected */}
          {isAll && (
              <div className="flex bg-gray-200 p-1 rounded-lg">
                {['global', 'local', 'lcs'].map((method) => (
                    <button
                        key={method}
                        onClick={() => setActiveTab(method)}
                        className={`px-4 py-1 text-sm font-semibold rounded-md transition-colors ${
                            activeTab === method ? 'bg-white text-blue-600 shadow' : 'text-gray-600 hover:text-gray-900'
                        }`}
                    >
                      {method.toUpperCase()}
                    </button>
                ))}
              </div>
          )}
        </div>

        {/* Plot Container */}
        <div className="p-4">
          {loading && <p className="text-blue-500 animate-pulse mb-4 font-semibold text-center">Downloading and parsing binary matrix...</p>}
          {error && <p className="text-red-500 bg-red-50 p-4 rounded text-center font-semibold">{error}</p>}

          <div
              ref={ref}
              className="w-full h-[500px] border border-gray-100 rounded"
              style={{ display: error || loading ? 'none' : 'block' }}
          />
        </div>
      </div>
  );
}

export default MatrixVisualizer;