import zarr
import numpy as np
import napari

# Load Zarr array
z = zarr.open('../results/protein/lcs_traceback.zarr', mode='r')
char_arr = np.array(z)  # shape: (H, W), dtype: object or str

# Map characters to integers safely
char_to_int = {'D': 1, 'U': 2, 'L': 3, ' ': 0}
# Use np.vectorize with fallback for unexpected chars
vectorized = np.vectorize(lambda c: char_to_int.get(c, 0), otypes=[np.uint8])
numeric = vectorized(char_arr)

# Visualize with napari
viewer = napari.Viewer()
viewer.add_image(numeric, colormap='viridis', name='LCS Traceback')
napari.run()
