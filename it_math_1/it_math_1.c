#include "it_math_1.h"

double randfrom(double min, double max) {
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int parse_args(int argc, const char **argv, params_t **params) {
    char *name = NULL;
    *params = calloc(1, sizeof(params_t));
    if (!(*params)) {
        return -ENOMEM;
    }
    params_t *params_init = *params;

    //init values
    name = "task_1";
    *params_init = (params_t) {NULL, 100, 51, 0.0001, 0, {0, 1}, 6, 0xEBAC0C, 100, -100, false};
    int algo = -1;

    const char *const usage[] = {
            "poisson [options]",
            NULL,
    };

    struct argparse_option options[] = {
            OPT_HELP(),
            OPT_GROUP("Basic options"),
            OPT_STRING('t', "task", &name, "name of task"),
            OPT_INTEGER('n', "size", &params_init->n, "size of side of the f function matrix"),
            OPT_INTEGER('c', "chunk_size", &params_init->chunk_size, "size of side of WaveChunk algo's chunk"),
            OPT_FLOAT('e', "epsilon", &params_init->epsilon, "accuracy of calculations"),
            OPT_FLOAT('l', "left_border", &params_init->borders.left, "left border"),
            OPT_FLOAT('r', "right_border", &params_init->borders.right, "right border"),
            OPT_INTEGER('a', "algo", &algo, "algorithm: 1-Sequential; 3-Parallel String; 6-Wave Chunk"),
            OPT_INTEGER('s', "seed", &params_init->random_seed, "random seed"),
            OPT_INTEGER('\0', "max", &algo, "random max"),
            OPT_INTEGER('\0', "min", &algo, "random min"),
            OPT_BOOLEAN('\0', "time-only", &params_init->time_only, "print only time of algo execution in seconds"),
            OPT_END()
    };


    struct argparse argparse;
    argparse_init(&argparse, options, usage, 0);
    argparse_describe(&argparse, "\nIt works and I'm so tired to write more\n", "");
    argparse_parse(&argparse, argc, argv);

    if (algo >= 0) {
        params_init->algo = algo;
    }
    params_init->h = fabs(params_init->borders.right - params_init->borders.left) / (params_init->n + 1);

    if (find_task(name, &(params_init->task))) {
        error("Cannot find given task\n");
        return -EINVAL;
    }

    return 0;
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
            double x = params->borders.left + (i / (double) (params->n + 1));
            double y = params->borders.left + (j / (double) (params->n + 1));
            double new_val;
            if (y == params->borders.left)
                new_val = params->task->bottom(x);
            else if (x == params->borders.left)
                new_val = params->task->left(y);
            else if (y == params->borders.right)
                new_val = params->task->upper(x);
            else if (x == params->borders.right)
                new_val = params->task->right(y);
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
    char output_path[OUTPUT_FILE_MAX_SIZE + 1];
    if (snprintf(output_path, OUTPUT_FILE_MAX_SIZE + 1, "%s_algo=%d_n=%d_task=%s%s", OUTPUT_FILE_PREFIX,
                 params->algo, params->n, params->task->name, OUTPUT_FILE_POSTFIX) <= 0) {
        error("Error in generating the output file path\n");
        return -ENAMETOOLONG;
    }
    FILE *output_file = fopen(output_path, "w");
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
    int rc = 0;
    params_t *params = NULL;
    double **f_matrix = NULL;
    double **u_matrix = NULL;

    rc = parse_args(argc, argv, &params);
    if (rc) {
        goto clean;
    }

    srand(params->random_seed);

    if (allocate_empty_matrix(&f_matrix, params)) {
        rc = -ENOMEM;
        goto clean;
    }

    rc = params->task->f(f_matrix, (int) params->n);
    if (rc) {
        error("Error in the f function of task: \"%s\"", params->task->name);
        goto clean;
    }

    if (allocate_empty_matrix(&u_matrix, params)) {
        rc = -ENOMEM;
        goto clean;
    }

    fill_u_matrix(u_matrix, params);

    double t1, t2;
    t1 = omp_get_wtime();

    int iterations;
    switch (params->algo) {
        case Sequential:
            if (!params->time_only)
                printf("run sequential algorithm\n");
            iterations = sequential_poisson(u_matrix, f_matrix, params);
            break;
        case ParallelString:
            if (!params->time_only)
                printf("run parallel string algorithm\n");
            iterations = parallel_string_poisson(u_matrix, f_matrix, params);
            break;
        case WaveChunk:
            if (!params->time_only)
                printf("run wave chunk algorithm\n");
            iterations = wave_chunk_poisson(u_matrix, f_matrix, params);
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
        if (!params->time_only) {
            printf("iterations: %d\n", iterations);
            printf("time: %fs\n", t2 - t1);
        } else {
            printf("%f\n", t2 - t1);
        }
    }

    rc = print_results(u_matrix, params);

    clean:
    free_matrix(f_matrix, params);
    free_matrix(u_matrix, params);
    free(params);
    return rc;
}
