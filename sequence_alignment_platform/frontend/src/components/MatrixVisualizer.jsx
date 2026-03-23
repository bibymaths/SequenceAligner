import React, { useEffect, useRef } from 'react';
import Plotly from 'plotly.js-dist-min';

// This component renders a DP heatmap using Plotly.js. It accepts a 2D
// array of scores and overlays a hypothetical traceback path. For brevity,
// the matrix is generated randomly here. In production, you would fetch
// downsampled matrix data from the backend and the traceback path from
// path files. Plotly handles zooming and panning out of the box.

function MatrixVisualizer() {
  const ref = useRef(null);
  useEffect(() => {
    // Generate dummy heatmap
    const size = 50;
    const z = [];
    for (let i = 0; i < size; i++) {
      const row = [];
      for (let j = 0; j < size; j++) {
        row.push(Math.sin(i / 3) * Math.cos(j / 5));
      }
      z.push(row);
    }
    const trace = {
      x: Array.from({ length: size }, (_, i) => i),
      y: Array.from({ length: size }, (_, i) => i),
      z: z,
      type: 'heatmap',
      colorscale: 'Viridis',
    };
    const pathX = Array.from({ length: size }, (_, i) => i);
    const pathY = Array.from({ length: size }, (_, i) => i);
    const trace2 = {
      x: pathX,
      y: pathY,
      mode: 'lines',
      line: { color: 'red', width: 2 },
      type: 'scatter',
    };
    Plotly.newPlot(ref.current, [trace, trace2], { margin: { t: 30 }, hovermode: false });
  }, []);
  return (
    <div className="border p-2 mt-4 bg-white">
      <h3 className="font-semibold mb-2">DP Matrix Visualizer</h3>
      <div ref={ref} style={{ width: '100%', height: '300px' }} />
    </div>
  );
}

export default MatrixVisualizer;