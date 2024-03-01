#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include <minmax.h>

#define OUTPUT_FILE "./result.txt"
#define error(...) fprintf(stderr, __VA_ARGS__)

#define N 100
#define EPSILON 0.0001
#define SEED 0xEBAC0C

#define CHUNK_SIZE 51

double randfrom(double min, double max) {
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

enum algo {
    sequential = 1,
    parallel_string = 3,
    wave_chunk = 6
};

typedef struct {
    double epsilon;
    double h;
    uint32_t n;
    uint32_t chunk_size;
    enum algo;
    int seed;
};

int parse_args(int argc, char **argv) {

    for (int i = 0; i < argc; i ++) {

    }
}

int allocate_empty_matrix(double ***new_matrix) {
    double **matrix = calloc(N + 2, sizeof(double *));
    if (!matrix)
        return -ENOMEM;

    for (int i = 0; i < N + 2; i++) {
        matrix[i] = calloc(N + 2, sizeof(double));
        if (!matrix[i])
            return -ENOMEM;
    }
    *new_matrix = matrix;
    return 0;
}

void fill_u_matrix(double **matrix) {
    for (int i = 0; i < N + 2; i++)
        for (int j = 0; j < N + 2; j++) {
            double x = i / (double) (N + 1);
            double y = j / (double) (N + 1);
            double new_val = 0;
            if (y == 0)
                new_val = 100 - 200 * x;
            else if (x == 0)
                new_val = 100 - 200 * y;
            else if (y == 1)
                new_val = -100 + 200 * x;
            else if (x == 1)
                new_val = -100 + 200 * y;
            else new_val = randfrom(-100, 100);
            matrix[i][j] = new_val;
        }
}

int sequential_poisson(double **u, double **f) {
    int iterations = 0;
    const double h = 1 / (N + 1);
    double dmax, temp, dm;
    do {
        dmax = 0;
        for (int i = 1; i < N + 1; i++)
            for (int j = 1; j < N + 1; j++) {
                temp = u[i][j];
                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - h * h * f[i][j]);
                dm = fabs(temp - u[i][j]);
                if (dmax < dm)
                    dmax = dm;
            }
        iterations += 1;
    } while (dmax > EPSILON);
    return iterations;
}

int parallel_string_poisson(double **u, double **f) { //11.3
    omp_lock_t dmax_lock;
    omp_init_lock(&dmax_lock);
    int iterations = 0;
    const double h = 1 / (N + 1);
    int i, j;
    double dmax, temp, dm, d;
    do {
        dmax = 0;
#pragma omp parallel for shared(u, dmax) private(i, j, temp, d, dm)
        for (i = 1; i < N + 1; i++) {
            dm = 0;
            for (j = 1; j < N + 1; j++) {
                temp = u[i][j];
                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - h * h * f[i][j]);
                d = fabs(temp - u[i][j]);
                if (dm < d)
                    dm = d;
            }
            omp_set_lock(&dmax_lock);
            if (dmax < dm) dmax = dm;
            omp_unset_lock(&dmax_lock);
        }
        iterations += 1;
    } while (dmax > EPSILON);
    return iterations;
}

typedef struct {
    int x;
    int y;
} chunk_t;

double process_chunk(double **u, double **f, chunk_t chunk) {
    const double h = 1 / (N + 1);
    double dmax, temp, dm;
    dmax = 0;
    for (int i = chunk.x * CHUNK_SIZE; i < (chunk.x + 1) * CHUNK_SIZE; i++) {
        if (i == 0 || i == N + 1)
            continue;
        for (int j = chunk.y * CHUNK_SIZE; j < (chunk.y + 1) * CHUNK_SIZE; j++) {
            if (j == 0 || j == N + 1)
                continue;
            temp = u[i][j];
            u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - h * h * f[i][j]);
            dm = fabs(temp - u[i][j]);
            if (dmax < dm)
                dmax = dm;
        }
    }
    return dmax;
}

int wave_chunk_poisson(double **u, double **f) { //11.6
    if ((N + 2) % CHUNK_SIZE) {
        error("Invalid chunk_t size/n");
        return -EINVAL;
    }
#define NB (N+2) / CHUNK_SIZE
    omp_lock_t dmax_lock;
    omp_init_lock(&dmax_lock);
    int nx, i, j, iterations = 0;
    double d, dmax, dm[NB];
    int dmax_chunk = min(100, NB);
    do {
        dmax = 0;
        for (nx = 0; nx < NB; nx++) {
            dm[nx] = 0;
#pragma omp parallel for shared(nx) private(i, j)
            for (i = 0; i < nx + 1; i++) {
                j = nx - i;
                d = process_chunk(u, f, (chunk_t) {i, j});
                if (dm[i] < d) dm[i] = d;
            }
        }
        for (nx = NB - 1; nx > 0; nx--) {
#pragma omp parallel for shared(nx) private(i, j)
            for (i = NB-nx; i < NB; i++) {
                j = 2 * (NB - 1) - nx - i + 1;
                d = process_chunk(u, f, (chunk_t) {i, j});
                if (dm[i] < d) dm[i] = d;
            }
        }
#pragma omp parallel for shared(n, dm, dmax) private(i, d)
        for (i = 0; i < NB; i += dmax_chunk) {
            d = 0;
            for (j = i; j < i + dmax_chunk; j++)
                if (d < dm[j]) d = dm[j];
            omp_set_lock(&dmax_lock);
            if (dmax < d) dmax = d;
            omp_unset_lock(&dmax_lock);
        }
        iterations++;
    } while (dmax > EPSILON);
    return iterations;
}

int print_results(double **matrix) {
    int rc = 0;
    FILE *output_file = fopen(OUTPUT_FILE, "w");
    if (!output_file) {

        return -EIO;
    }

    for (int i = 0; i < N + 2; i++) {
        for (int j = 0; j < N + 2; j++) {
            if (fprintf(output_file, "%.2f ", matrix[i][j]) <= 0) {
                rc = -EIO;
                goto close;
            }
        }
        if (fprintf(output_file, "\n") <= 0) {
            rc = -EIO;
            goto close;
        }
    }
    close:
    fclose(output_file);
    return rc;
}

void free_matrix(double **matrix) {
    for (int i = 0; i < N + 2; i++)
        free(matrix[i]);
    free(matrix);
}

int main(int argc, char **argv) {
    double **f_matrix = NULL;
    double **u_matrix = NULL;
    int rc = 0;
    srand(SEED);

    if (allocate_empty_matrix(&f_matrix)) {
        rc = -ENOMEM;
        goto clean;
    }

    if (allocate_empty_matrix(&u_matrix)) {
        rc = -ENOMEM;
        goto clean;
    }

    fill_u_matrix(u_matrix);

    int iterations = 0;
//    iterations = sequential_poisson(u_matrix, f_matrix);
//    iterations = parallel_string_poisson(u_matrix, f_matrix);
    iterations = wave_chunk_poisson(u_matrix, f_matrix);
    if (iterations < 0) { //if iterations < 0 it's the error code of poisson
        rc = iterations;
        goto clean;
    } else {
        printf("%d", iterations);
    }

    rc = print_results(u_matrix);

    clean:
    free_matrix(f_matrix);
    free_matrix(u_matrix);
    return rc;
}
