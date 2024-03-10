import numpy as np

SEED = 111111
np.random.seed(SEED)

N = 100


def poisson(N, u, f, iterations):
    h = 1 / (N + 1)
    while iterations > 0:
        dmax = 0
        for i in range(1, N + 1):
            for j in range(1, N + 1):
                temp = u[i][j]
                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - h * h * f[i][j])
                dm = abs(temp - u[i][j])
                if dmax < dm:
                    dmax = dm
        iterations -= 1


def create_f_matrix(N: int):
    matrix = np.zeros((N + 2, N + 2))
    return matrix


f_matrix = create_f_matrix(N)


def create_u_matrix(N):
    matrix = np.random.uniform(-100, 100, [N + 2, N + 2])
    for i in range(N + 2):
        for j in range(N + 2):
            x = i / (N + 1)
            y = j / (N + 1)
            if y == 0:
                matrix[i][j] = 100 - 200 * x
            elif x == 0:
                matrix[i][j] = 100 - 200 * y
            elif y == 1:
                matrix[i][j] = -100 + 200 * x
            elif x == 1:
                matrix[i][j] = -100 + 200 * y
    return matrix


u_matrix = create_u_matrix(N)

poisson(N, u_matrix, f_matrix, 100)

np.savetxt("result_py.txt", u_matrix, fmt='%.2f', delimiter=' ', newline='\n', header='', footer='', encoding='ascii')
