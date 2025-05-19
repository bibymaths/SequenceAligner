import dask.array as da
import dask.dataframe as dd
import datashader as ds
import datashader.transfer_functions as tf
import argparse

def visualize_zarr_matrix(zarr_path, output_png="matrix_preview.png", cmap='inferno', width=1000, height=1000):
    """
    Efficiently visualize a 2D matrix from Zarr as a heatmap using Datashader.
    """
    arr = da.from_zarr(zarr_path)
    if arr.ndim != 2:
        raise ValueError("Expected a 2D matrix")

    h, w = arr.shape

    x = da.tile(da.arange(w), h)
    y = da.repeat(da.arange(h), w)
    val = arr.ravel()

    ddf = dd.from_dask_array(da.stack([x, y, val], axis=1), columns=["x", "y", "val"])

    cvs = ds.Canvas(plot_width=width, plot_height=height,
                    x_range=(0, w), y_range=(0, h))
    agg = cvs.points(ddf, 'x', 'y', agg=ds.mean('val'))
    img = tf.shade(agg, cmap=cmap)
    img.to_pil().save(output_png)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--zarrfile", required=True, help="Path to the .zarr file")
    parser.add_argument("--output", default="matrix_preview.png", help="Output image name")
    parser.add_argument("--cmap", default="inferno", help="Colormap")
    parser.add_argument("--width", type=int, default=1000, help="Canvas width")
    parser.add_argument("--height", type=int, default=1000, help="Canvas height")
    args = parser.parse_args()

    visualize_zarr_matrix(args.zarrfile, args.output, args.cmap, args.width, args.height)
