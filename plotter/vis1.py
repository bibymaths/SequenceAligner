#!/usr/bin/env python3
import os
# Choose a safe Numba threading backend; adjust to "tbb" if you have it installed
os.environ["NUMBA_THREADING_LAYER"] = "omp"

import sys
import xarray as xr
import numpy as np
from ome_zarr.writer import write_image
from ome_zarr.io import parse_url
import datashader as ds
import datashader.transfer_functions as tf

def convert_to_ome_zarr(input_zarr, ome_zarr_dir):
    ds = xr.open_zarr(input_zarr, consolidated=True)
    arr = ds["matrix"].values
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D matrix, got shape {arr.shape}")
    # build 5D array: (t=1, c=1, z=1, y, x)
    data = arr[np.newaxis, np.newaxis, np.newaxis, :, :].astype(np.uint8)
    store = parse_url(ome_zarr_dir, mode="w").store
    write_image(store, data,
                dimension_separator="/",
                scale_factors=[[2, 2]],
                axes="tczyx")

def plot_multiscale_zarr(ome_zarr_dir, output_png, px=1000, py=1000):
    ds = xr.open_zarr(ome_zarr_dir, consolidated=True)
    da = ds["0"]  if "0" in ds.data_vars else ds["matrix"]
    # depending on how write_image named your pyramid; adjust above if needed
    canvas = ds.Canvas(plot_width=px, plot_height=py)
    agg = canvas.raster(da)
    img = tf.shade(agg, cmap=["black", "blue", "cyan", "white"])
    img.to_pil().save(output_png)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python process_and_plot_zarr.py <input.zarr> <ome.zarr> <output.png> <pixels>")
        print("  pixels: WIDTHxHEIGHT, e.g. 1000x1000")
        sys.exit(1)

    inp, ome, out_png, pxpy = sys.argv[1:]
    w, h = map(int, pxpy.lower().split("x"))
    print(f"→ Converting {inp} → multiscale Zarr at {ome}")
    convert_to_ome_zarr(inp, ome)
    print(f"→ Rendering heatmap → {out_png} ({w}×{h})")
    plot_multiscale_zarr(ome, out_png, px=w, py=h)
    print("Done.")
