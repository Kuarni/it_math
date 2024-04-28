import time
from random import normalvariate

import numpy as np
from PIL import Image

from .SVD import SVD


class PowerMethodSVD(SVD):
    name = "pwmSVD"

    def __init__(self, compression_degree=2, max_iter=0, tolerance=1e-2, start_vector_seed=0xebac0c):
        self.max_iter = max_iter
        self.tolerance = tolerance
        np.random.seed(start_vector_seed)
        super().__init__(compression_degree)

    def __svd_1d(self, matrix):
        def random_unit_vector(size):
            unnormalized = np.random.uniform(0, 1, size)
            return unnormalized / np.linalg.norm(unnormalized)

        n, m = matrix.shape
        x = random_unit_vector(m)
        current_vector = x
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
            start = time.time()
            matrix_for_1d = matrix.copy()

            for u, singular_value, v in svd_so_far:
                matrix_for_1d -= singular_value * np.outer(u, v)

            v = self.__svd_1d(matrix_for_1d)
            u_unnormalized = np.dot(matrix, v)
            sigma = np.linalg.norm(u_unnormalized)
            u = u_unnormalized / sigma

            svd_so_far.append((u, sigma, v))
            print(time.time() - start, i)

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
