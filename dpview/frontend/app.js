const canvas = document.getElementById("heatmap");
const ctx = canvas.getContext("2d");

async function drawDownsampled() {
  try {
    const res = await fetch("http://localhost:8000/downsampled");
    const data = await res.json();

    if (!data || !data.length) {
      console.error("No matrix data received");
      return;
    }

    console.log(`Drawing matrix: ${data.length} rows × ${data[0].length} cols`);

    const rows = data.length;
    const cols = data[0].length;

    const cellWidth = canvas.width / cols;
    const cellHeight = canvas.height / rows;

    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        const v = data[i][j];
        const norm = Math.min(1, Math.max(0, (v + 10000) / 20000));
        const color = d3.interpolateInferno(norm);
        ctx.fillStyle = color;
        ctx.fillRect(j * cellWidth, i * cellHeight, cellWidth, cellHeight);
      }
    }
  } catch (e) {
    console.error("Draw failed:", e);
  }
}

drawDownsampled();
