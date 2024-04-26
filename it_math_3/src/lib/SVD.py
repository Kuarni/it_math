import math
from io import BufferedReader
from struct import pack, unpack, calcsize
from typing import NamedTuple

import PIL
import numpy as np

from src.lib.Packing import Packing
from src.lib.Transformer import Transformer


class SVD(Transformer):
    name = "SVD"

    def __init__(self, compression_degree=2):
        self.__u = None
        self.__v = None
        self.__s = None
        self.compression_degree = compression_degree

    def encode(self, image: PIL.Image):
        h, w = image.height, image.width
        saving_part = math.ceil(h * w // (4 * self.compression_degree * (h + w + 1)))
        u, s, v = np.linalg.svd(image, full_matrices=False)
        self.__u = u[:, :saving_part].astype(np.float32)
        self.__s = s[:saving_part].astype(np.float32)
        self.__v = v[:saving_part, :].astype(np.float32)
        self.__matrices_sizes.value = (h, saving_part, w)
        return self

    def size(self):
        if self.__u is None or self.__s is None or self.__v is None:
            raise ValueError('Cannot calculate SVD size because it\'s not initialized')
        return calcsize(self.__matrices_sizes.packing) + self.__u.nbytes + self.__s.nbytes + self.__v.nbytes

    def decode(self):
        return np.dot(np.dot(self.__u, np.diag(self.__s)), self.__v).astype(np.uint8)

    __matrices_sizes: Packing = Packing(None, "<III")

    def to_bytes(self):
        if self.__u is None or self.__s is None or self.__v is None:
            raise ValueError('Cannot convert to bytes an empty values')
        if self.__matrices_sizes.value is None:
            raise ValueError('Matrices\' sizes is not initialized')
        return pack(self.__matrices_sizes.packing,
                    *self.__matrices_sizes.value) + self.__u.tobytes() + self.__s.tobytes() + self.__v.tobytes()

    def read(self, f: BufferedReader):
        n, k, m = unpack(self.__matrices_sizes.packing, f.read(calcsize(self.__matrices_sizes.packing)))
        self.__u = np.fromfile(f, dtype=np.float32, count=n * k).reshape(n, k)
        self.__s = np.fromfile(f, dtype=np.float32, count=k)
        self.__v = np.fromfile(f, dtype=np.float32, count=k * m).reshape(k, m)
        return self
