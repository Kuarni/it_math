import math
from io import BufferedReader
from struct import pack, unpack, calcsize

from PIL import Image
import numpy as np

from .Packing import Packing
from .Transformer import Transformer


class SVD(Transformer):
    name = "SVD"

    def __init__(self, compression_degree=2):
        self._u = None
        self._v = None
        self._s = None
        self.compression_degree = compression_degree

    def _saving_part(self, h, w):
        return math.ceil(h * w // (4 * self.compression_degree * (h + w + 1)))

    def _standard_svd(self, matrix: np.ndarray, svd_algo):
        h, w = matrix.shape
        saving_part = self._saving_part(h, w)
        u, s, v = svd_algo(matrix, full_matrices=False)
        self._u = u[:, :saving_part].astype(np.float32)
        self._s = s[:saving_part].astype(np.float32)
        self._v = v[:saving_part, :].astype(np.float32)
        self._matrices_sizes.value = (h, saving_part, w)
        return self

    def encode(self, image: Image):
        return self._standard_svd(np.array(image), np.linalg.svd)

    def size(self):
        if self._u is None or self._s is None or self._v is None:
            raise ValueError('Cannot calculate SVD size because it\'s not initialized')
        return calcsize(self._matrices_sizes.packing) + self._u.nbytes + self._s.nbytes + self._v.nbytes

    def decode(self):
        return np.dot(np.dot(self._u, np.diag(self._s)), self._v).astype(np.uint8)

    _matrices_sizes: Packing = Packing(None, "<III")

    def to_bytes(self):
        if self._u is None or self._s is None or self._v is None:
            raise ValueError('Cannot convert to bytes an empty values')
        if self._matrices_sizes.value is None:
            raise ValueError('Matrices\' sizes is not initialized')
        return pack(self._matrices_sizes.packing,
                    *self._matrices_sizes.value) + self._u.tobytes() + self._s.tobytes() + self._v.tobytes()

    def read(self, f: BufferedReader):
        n, k, m = unpack(self._matrices_sizes.packing, f.read(calcsize(self._matrices_sizes.packing)))
        self._u = np.fromfile(f, dtype=np.float32, count=n * k).reshape(n, k)
        self._s = np.fromfile(f, dtype=np.float32, count=k)
        self._v = np.fromfile(f, dtype=np.float32, count=k * m).reshape(k, m)
        return self
