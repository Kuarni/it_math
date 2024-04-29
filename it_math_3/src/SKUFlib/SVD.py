import math
from abc import ABC, abstractmethod
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
        return np.clip(np.dot(np.dot(self._u, np.diag(self._s)), self._v), 0, 255).astype(np.uint8)

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

    @abstractmethod
    def _algo(self, matrix):
        pass

    def encode(self, image: Image):
        matrix = np.array(image).astype(np.float64)
        h, w = matrix.shape
        u, s, v = self._algo(matrix)
        self._u = u.astype(np.float32)
        self._s = s.astype(np.float32)
        self._v = v.astype(np.float32)
        self._matrices_sizes.value = (h, len(s), w)
        return self


class NpSVD(SVD):
    name = "npSVD"

    def _algo(self, matrix):
        u, s, v = np.linalg.svd(matrix, full_matrices=False)
        saving_part = self._saving_part(*matrix.shape)
        return u[:, :saving_part], s[:saving_part], v[:saving_part, :]


class PowSVD(SVD, ABC):
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


class PowerMethodSVD(PowSVD):
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

    def _algo(self, matrix):
        s = self._saving_part(matrix.shape[1], matrix.shape[0])
        svd_so_far = []

        for i in range(s):
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


class BlockPowerMethodSVD(PowSVD, ABC):
    name = "bpmSVD"

    def _algo(self, matrix):
        h, w = matrix.shape
        s = self._saving_part(h, w)
        v = np.array([self._random_unit_vector(s) for _ in range(w)])
        u = np.zeros((h, s))
        sigma = np.zeros((s, s))
        err = 10 ** 9
        while err > self.tolerance:
            q, r = np.linalg.qr(np.dot(matrix, v))
            u = q[:, :s]
            q, r = np.linalg.qr(np.dot(matrix.T, u))
            v = q[:, :s]
            sigma = r[:s, :s]
            err = np.linalg.norm(np.dot(matrix, v) - np.dot(u, sigma))
        return u, np.diag(sigma), v.T
