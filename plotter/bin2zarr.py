#!/usr/bin/env python3
import numpy as np
import xarray as xr
import argparse

def bin_to_zarr_int32(binfile, zarrfile):
    """
    Convert an int32 binary file to a Zarr-backed xarray Dataset.
    """
    with open(binfile, 'rb') as f:
        rows = np.frombuffer(f.read(4), dtype=np.int32)[0]
        cols = np.frombuffer(f.read(4), dtype=np.int32)[0]
        data = np.frombuffer(f.read(), dtype=np.int32).reshape((rows, cols))

    # Wrap in xarray and write Zarr
    da = xr.DataArray(data, dims=("y","x"), name="matrix")
    ds = da.to_dataset()
    ds.to_zarr(zarrfile, mode="w", consolidated=True)
    print(f"Written int32 matrix to {zarrfile}")

def bin_to_zarr_char(binfile, zarrfile):
    """
    Convert a uint8 (char) binary file to a Zarr-backed xarray Dataset.
    """
    with open(binfile, 'rb') as f:
        rows = np.frombuffer(f.read(4), dtype=np.int32)[0]
        cols = np.frombuffer(f.read(4), dtype=np.int32)[0]
        data = np.frombuffer(f.read(), dtype=np.uint8).reshape((rows, cols))

    da = xr.DataArray(data, dims=("y","x"), name="matrix")
    ds = da.to_dataset()
    ds.to_zarr(zarrfile, mode="w", consolidated=True)
    print(f"Written char matrix to {zarrfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert binary file to zarr-backed xarray dataset.')
    parser.add_argument('--binfile', type=str,  required=True, help='Input binary file')
    parser.add_argument('--zarrfile', type=str, required=True, help='Output Zarr directory')
    parser.add_argument('--type', choices=['int32','char'], default='int32',
                        help="Data type: 'int32' or 'char'")
    args = parser.parse_args()

    if args.type == 'int32':
        bin_to_zarr_int32(args.binfile, args.zarrfile)
    else:
        bin_to_zarr_char(args.binfile, args.zarrfile)