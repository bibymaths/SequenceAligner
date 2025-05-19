import zarr
import xarray as xr
import datashader as ds
import datashader.transfer_functions as tf

# Path to the converted Zarr file
zarr_path = "your_converted_matrix.zarr"  # replace with actual output from bin_to_zarr

# Step 1: Open the Zarr file using xarray
ds_zarr = xr.open_zarr(zarr_path)

# Step 2: Extract the array variable (Zarr saves with default variable name 'arr_0')
data = ds_zarr["arr_0"]

# Step 3: Plot with Datashader (downsample to 1000x1000 pixels)
canvas = ds.Canvas(plot_width=1000, plot_height=1000)
agg = canvas.raster(data)

# Step 4: Apply shading and save
img = tf.shade(agg, cmap=["black", "red", "yellow", "white"])
img.to_pil().save("heatmap_from_zarr.png")
