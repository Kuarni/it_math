import numpy as np

path_a = "result_py.txt"
# path_b = "C:/Users/kotyr/CLionProjects/it_math_1/cmake-build-debug/result.txt"
path_b = "C:/Users/kotyr/Downloads/LENA_FINAL_BREATH (4).txt"
delta = 0.5
N = 100

matrix_a = np.loadtxt(path_a)
matrix_b = np.loadtxt(path_b)

max_delta = 0

for i in range(N+2):
    for j in range(N+2):
        cur_delta = abs(matrix_b[i][j] - matrix_a[i][j])
        if cur_delta > delta:
            print(f"pizdec in position {i},{j}: in a = {matrix_a[i][j]}, in b = {matrix_b[i][j]}")
        if cur_delta > max_delta:
            max_delta = cur_delta

print(f"max_delta = {max_delta}")