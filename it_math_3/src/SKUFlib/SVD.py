import math
from abc import ABC
from struct import pack, unpack, calcsize
from typing import BinaryIO

from PIL import Image
import numpy as np

from .Packing import Packing
from .Transformer import Transformer


class SVD(Transformer, ABC):
    def __init__(self, compression_degree: int = 2):
        self._u = None
        self._v = None
        self._s = None
        self.compression_degree = compression_degree

    def _saving_part(self, h, w):
        return math.ceil(h * w // (4 * self.compression_degree * (h + w + 1)))

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

    def read(self, f: BinaryIO):
        n, k, m = unpack(self._matrices_sizes.packing, f.read(calcsize(self._matrices_sizes.packing)))
        self._u = np.fromfile(f, dtype=np.float32, count=n * k).reshape(n, k)
        self._s = np.fromfile(f, dtype=np.float32, count=k)
        self._v = np.fromfile(f, dtype=np.float32, count=k * m).reshape(k, m)
        self._matrices_sizes.value = (n, k, m)
        return self


class NpSVD(SVD):
    name = "npSVD"

    def encode(self, image: Image):
        matrix = np.array(image)
        h, w = matrix.shape
        saving_part = self._saving_part(h, w)
        u, s, v = np.linalg.svd(matrix, full_matrices=False)
        self._u = u[:, :saving_part].astype(np.float32)
        self._s = s[:saving_part].astype(np.float32)
        self._v = v[:saving_part, :].astype(np.float32)
        self._matrices_sizes.value = (h, saving_part, w)
        return self


class MySVD(SVD, ABC):
    def __init__(self, compression_degree=2, max_iter=0, tolerance=1e-2, start_vector_seed=0xebac0c):
        """
        :param max_iter: maximum number of iterations for 1d vector computation. If <= 0 when there is no limits
        :param tolerance: for 1d vector computation.
        :param start_vector_seed: seed for random number generator for 1d vector init.
        """
        self.max_iter = max_iter
        self.tolerance = tolerance
        np.random.seed(start_vector_seed)
        super().__init__(compression_degree)

    @staticmethod
    def _random_unit_vector(size):
        unnormalized = np.random.uniform(0, 1, size)
        return unnormalized / np.linalg.norm(unnormalized)


class PowerMethodSVD(MySVD):
    name = "pwmSVD"

    def __svd_1d(self, matrix):
        n, m = matrix.shape
        current_vector = super()._random_unit_vector(m)
        b = np.dot(matrix.T, matrix)

        iterations = self.max_iter if self.max_iter > 0 else -1
        while iterations != 0:
            iterations -= 1
            last_vector = current_vector
            current_vector = np.dot(b, last_vector)
            current_vector /= np.linalg.norm(current_vector)

            if (abs(np.dot(current_vector, last_vector))) > 1 - self.tolerance:
                break
        return current_vector

    def __algo(self, matrix, k):
        matrix = matrix.astype(np.float64)
        svd_so_far = []

        for i in range(k):
            matrix_for_1d = matrix.copy()

            for u, singular_value, v in svd_so_far:
                matrix_for_1d -= singular_value * np.outer(u, v)

            v = self.__svd_1d(matrix_for_1d)
            u_unnormalized = np.dot(matrix, v)
            sigma = np.linalg.norm(u_unnormalized)
            u = u_unnormalized / sigma

            svd_so_far.append((u, sigma, v))

        us, singular_value, vs = [np.array(x) for x in zip(*svd_so_far)]
        return us.T, singular_value, vs

    def encode(self, image: Image):
        matrix = np.array(image)
        h, w = matrix.shape
        saving_part = super()._saving_part(h, w)
        u, s, v = self.__algo(matrix, saving_part)
        self._u = u.astype(np.float32)
        self._s = s.astype(np.float32)
        self._v = v.astype(np.float32)
        self._matrices_sizes.value = (h, saving_part, w)
        return self
