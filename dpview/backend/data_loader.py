import numpy as np

HEADER_SIZE = 8  # bytes (2 int32s)

def load_matrix(path):
    with open(path, 'rb') as f:
        header = np.frombuffer(f.read(8), dtype='<i4')
    rows, cols = header
    mat = np.memmap(path, dtype='<i4', mode='r', offset=HEADER_SIZE, shape=(rows, cols))
    return mat, rows, cols