import numpy as np

path_a = "result_py.txt"
path_b = "C:/Users/kotyr/CLionProjects/it_math/it_math_1/cmake-build-debug/result_algo=6_n=100_task=task_1.txt"
# path_b = "C:/Users/kotyr/Downloads/LENA_FINAL_BREATH (9).txt"
delta = 0.5

matrix_a = np.loadtxt(path_a)
matrix_b = np.loadtxt(path_b)

max_delta = [0, (0, 0), 0, 0]

for i in range(len(matrix_a)):
    for j in range(len(matrix_a[0])):
        cur_delta = abs(matrix_b[i][j] - matrix_a[i][j])
        if cur_delta > delta:
            print(f"pizdec in position {i},{j}: in a = {matrix_a[i][j]}, in b = {matrix_b[i][j]}")
        if cur_delta > max_delta[0]:
            max_delta[0] = cur_delta
            max_delta[1] = (i, j)
            max_delta[2] = matrix_a[i][j]
            max_delta[3] = matrix_b[i][j]

print(
    f"max_delta = {max_delta[0]}, position = ({max_delta[1][0]},{max_delta[1][0]}), value = {max_delta[2]}, {max_delta[3]}")
