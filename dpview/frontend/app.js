const app = new PIXI.Application({ width: 800, height: 800, backgroundColor: 0x111111, resolution: 1 });
document.body.appendChild(app.view);

let offsetX = 0, offsetY = 0;
const tileSize = 100, cellSize = 8;

async function fetchTile(x, y) {
  const res = await fetch(`http://localhost:8000/tile?x=${x}&y=${y}&width=${tileSize}&height=${tileSize}`);
  return res.json();
}

function valueToColor(v) {
  const norm = Math.min(255, Math.max(0, Math.floor((v + 5000) / 10000 * 255)));
  return PIXI.utils.rgb2hex([norm / 255, norm / 255, norm / 255]);
}

async function drawTile(x, y) {
  const data = await fetchTile(x, y);
  const g = new PIXI.Graphics();
  for (let i = 0; i < data.length; i++) {
    for (let j = 0; j < data[i].length; j++) {
      g.beginFill(valueToColor(data[i][j]));
      g.drawRect(j * cellSize, i * cellSize, cellSize, cellSize);
      g.endFill();
    }
  }
  app.stage.removeChildren();
  app.stage.addChild(g);
}

drawTile(offsetX, offsetY);

document.addEventListener('keydown', e => {
  if (e.key === 'ArrowRight') offsetX += tileSize;
  if (e.key === 'ArrowLeft') offsetX = Math.max(0, offsetX - tileSize);
  if (e.key === 'ArrowDown') offsetY += tileSize;
  if (e.key === 'ArrowUp') offsetY = Math.max(0, offsetY - tileSize);
  drawTile(offsetX, offsetY);
});