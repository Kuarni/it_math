#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <omp.h>
#include <inttypes.h>
#include "argparse.h"

#define OUTPUT_FILE "./result.txt"
#define error(...) fprintf(stderr, __VA_ARGS__)

#define min(a, b) a < b ? a : b

double randfrom(double min, double max) {
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

enum algo {
    Sequential = 1,
    ParallelString = 3,
    WaveChunk = 6
};

typedef struct {
    //size of side of the f function matrix
    uint32_t n;
    //size of side of WaveChunk algo's chunk
    uint32_t chunk_size;
    double epsilon;
    //interval between values
    double h;
    struct borders {
        double left;
        double right;
    } borders;
    enum algo algo;
    int random_seed;
    int random_max;
    int random_min;
} params_t;

params_t parse_args(int argc, const char **argv) {
    params_t params = {100, 51, 0.0001, 0, {0, 1}, 6, 0xEBAC0C, 100, -100};
    int algo = -1;

    const char *const usage[] = {
            "poisson [options]",
            NULL,
    };

    struct argparse_option options[] = {
            OPT_HELP(),
            OPT_GROUP("Basic options"),
            OPT_INTEGER('n', "size", &params.n, "size of side of the f function matrix"),
            OPT_INTEGER('c', "chunk_size", &params.chunk_size, "size of side of WaveChunk algo's chunk"),
            OPT_FLOAT('e', "epsilon", &params.epsilon, "accuracy of calculations"),
            OPT_FLOAT('l', "left_border", &params.borders.left, "left border"),
            OPT_FLOAT('r', "right_border", &params.borders.right, "right border"),
            OPT_INTEGER('a', "algo", &algo, "algorithm: 1-Sequential; 3-Parallel String; 6-Wave Chunk"),
            OPT_INTEGER('s', "seed", &params.random_seed, "random seed"),
            OPT_INTEGER('\0', "max", &algo, "random max"),
            OPT_INTEGER('\0', "min", &algo, "random min"),
            OPT_END()
    };


    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argparse_describe(&argparse, "\nIt works and I'm so tired to write more\n", "");
    argparse_parse(&argparse, argc, argv);

    if (algo >= 0) {
        params.algo = algo;
    }
    params.h = fabs(params.borders.right - params.borders.left) / (params.n + 1);

    return params;
}

int allocate_empty_matrix(double ***new_matrix, const params_t *params) {
    double **matrix = calloc(params->n + 2, sizeof(double *));
    if (!matrix)
        return -ENOMEM;

    for (int i = 0; i < params->n + 2; i++) {
        matrix[i] = calloc(params->n + 2, sizeof(double));
        if (!matrix[i])
            return -ENOMEM;
    }
    *new_matrix = matrix;
    return 0;
}

void fill_u_matrix(double **matrix, const params_t *params) {
    for (int i = 0; i < params->n + 2; i++)
        for (int j = 0; j < params->n + 2; j++) {
            double x =  params->borders.left + (i / (double) (params->n + 1));
            double y = params->borders.left + (j / (double) (params->n + 1));
            double new_val = 0;
            if (y == params->borders.left)
                new_val = 100 - 200 * x;
            else if (x == params->borders.left)
                new_val = 100 - 200 * y;
            else if (y == params->borders.right)
                new_val = -100 + 200 * x;
            else if (x == params->borders.right)
                new_val = -100 + 200 * y;
            else new_val = randfrom(params->random_min, params->random_max);
            matrix[i][j] = new_val;
        }
}

int sequential_poisson(double **u, double **f, const params_t *params) {
    int iterations = 0;
    double dmax, temp, dm, h = params->h;
    do {
        dmax = 0;
        for (int i = 1; i < params->n + 1; i++)
            for (int j = 1; j < params->n + 1; j++) {
                temp = u[i][j];
                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - h * h * f[i][j]);
                dm = fabs(temp - u[i][j]);
                if (dmax < dm)
                    dmax = dm;
            }
        iterations += 1;
    } while (dmax > params->epsilon);
    return iterations;
}

int parallel_string_poisson(double **u, double **f, const params_t *params) { //11.3
    omp_lock_t dmax_lock;
    omp_init_lock(&dmax_lock);
    int iterations = 0;
    int i, j;
    double dmax, temp, dm, d, h = params->h;
    do {
        dmax = 0;
#pragma omp parallel for shared(u, dmax) private(i, j, temp, d, dm)
        for (i = 1; i < params->n + 1; i++) {
            dm = 0;
            for (j = 1; j < params->n + 1; j++) {
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
    } while (dmax > params->epsilon);
    return iterations;
}

typedef struct {
    uint32_t x;
    uint32_t y;
} chunk_t;

double process_chunk(double **u, double **f, chunk_t chunk, const params_t *params) {
    double dmax, temp, dm, h = params->h;
    dmax = 0;
    for (uint32_t i = chunk.x * params->chunk_size; i < (chunk.x + 1) * params->chunk_size; i++) {
        if (i == 0 || i == params->n + 1)
            continue;
        for (uint32_t j = chunk.y * params->chunk_size; j < (chunk.y + 1) * params->chunk_size; j++) {
            if (j == 0 || j == params->n + 1)
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

int wave_chunk_poisson(double **u, double **f, const params_t *params) { //11.6
    if ((params->n + 2) % params->chunk_size) {
        error("Invalid chunk size: %d\n", params->chunk_size);
        return -EINVAL;
    }
    const uint32_t nb = (params->n + 2) / params->chunk_size;
    omp_lock_t dmax_lock;
    omp_init_lock(&dmax_lock);
    uint32_t nx, i, j, iterations = 0;
    double d, dmax, dm[nb];
    int dmax_chunk = min(100, nb);
    do {
        dmax = 0;
        for (nx = 0; nx < nb; nx++) {
            dm[nx] = 0;
#pragma omp parallel for shared(nx) private(i, j)
            for (i = 0; i < nx + 1; i++) {
                j = nx - i;
                d = process_chunk(u, f, (chunk_t) {i, j}, params);
                if (dm[i] < d) dm[i] = d;
            }
        }
        for (nx = nb - 1; nx > 0; nx--) {
#pragma omp parallel for shared(nx) private(i, j)
            for (i = nb - nx; i < nb; i++) {
                j = 2 * (nb - 1) - nx - i + 1;
                d = process_chunk(u, f, (chunk_t) {i, j}, params);
                if (dm[i] < d) dm[i] = d;
            }
        }
#pragma omp parallel for shared(n, dm, dmax) private(i, d)
        for (i = 0; i < nb; i += dmax_chunk) {
            d = 0;
            for (j = i; j < i + dmax_chunk; j++)
                if (d < dm[j]) d = dm[j];
            omp_set_lock(&dmax_lock);
            if (dmax < d) dmax = d;
            omp_unset_lock(&dmax_lock);
        }
        iterations++;
    } while (dmax > params->epsilon);
    return (int) iterations;
}

int print_results(double **matrix, const params_t *params) {
    int rc = 0;
    FILE *output_file = fopen(OUTPUT_FILE, "w");
    if (!output_file) {

        return -EIO;
    }

    for (int i = 0; i < params->n + 2; i++) {
        for (int j = 0; j < params->n + 2; j++) {
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

void free_matrix(double **matrix, const params_t *params) {
    for (int i = 0; i < params->n + 2; i++)
        free(matrix[i]);
    free(matrix);
}

int main(int argc, const char **argv) {
    const params_t params = parse_args(argc, argv);
    double **f_matrix = NULL;
    double **u_matrix = NULL;
    int rc = 0;
    srand(params.random_seed);

    if (allocate_empty_matrix(&f_matrix, &params)) {
        rc = -ENOMEM;
        goto clean;
    }

    if (allocate_empty_matrix(&u_matrix, &params)) {
        rc = -ENOMEM;
        goto clean;
    }

    fill_u_matrix(u_matrix, &params);

    double t1, t2;
    t1 = omp_get_wtime();

    int iterations;
    switch (params.algo) {
        case Sequential:
            printf("run sequential algorithm\n");
            iterations = sequential_poisson(u_matrix, f_matrix, &params);
            break;
        case ParallelString:
            printf("run parallel string algorithm\n");
            iterations = parallel_string_poisson(u_matrix, f_matrix, &params);
            break;
        case WaveChunk:
            printf("run wave chunk algorithm\n");
            iterations = wave_chunk_poisson(u_matrix, f_matrix, &params);
            break;
        default:
            error("Invalid algo was given\n");
            rc = -EINVAL;
            goto clean;
    }
    if (iterations < 0) { //if iterations < 0 it's the error code of poisson
        rc = iterations;
        goto clean;
    } else {
        t2 = omp_get_wtime();
        printf("iterations: %d\n", iterations);
        printf("time: %f\n", t2-t1);
    }

    rc = print_results(u_matrix, &params);

    clean:
    free_matrix(f_matrix, &params);
    free_matrix(u_matrix, &params);
    return rc;
}
