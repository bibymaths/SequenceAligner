import numpy as np
import zarr
import argparse

def bin_to_zarr_int32(binfile, zarrfile):
    """
    Convert a binary file to zarr format.
    Args:
        binfile (str): Path to the input binary file.
        zarrfile (str): Path to the output zarr file.
    """
    with open(binfile, 'rb') as f:
        rows = np.frombuffer(f.read(4), dtype=np.int32)[0]
        cols = np.frombuffer(f.read(4), dtype=np.int32)[0]
        data = np.frombuffer(f.read(), dtype=np.int32).reshape((rows, cols))

    zarr.save(zarrfile, data)

def bin_to_zarr_char(binfile, zarrfile):
    """
    Convert a binary file to zarr format.
    Args:
        binfile (str): Path to the input binary file.
        zarrfile (str): Path to the output zarr file.

    """
    with open(binfile, 'rb') as f:
        rows = np.frombuffer(f.read(4), dtype=np.int32)[0]
        cols = np.frombuffer(f.read(4), dtype=np.int32)[0]
        data = np.frombuffer(f.read(), dtype=np.uint8).reshape((rows, cols))

    zarr.save(zarrfile, data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert binary file to zarr format.')
    parser.add_argument('--binfile', type=str, help='Input binary file')
    parser.add_argument('--zarrfile', type=str, help='Output zarr file')
    parser.add_argument('--type', type=str, choices=['int32', 'char'], default='int32', help='Data type of the binary file')
    args = parser.parse_args()

    if args.type == 'int32':
        bin_to_zarr_int32(args.binfile, args.zarrfile)
    elif args.type == 'char':
        bin_to_zarr_char(args.binfile, args.zarrfile)
    else:
        raise ValueError("Unsupported data type. Use 'int32' or 'char'.")

