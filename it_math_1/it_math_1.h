#ifndef POISSON_IT_MATH_1_H
#define POISSON_IT_MATH_1_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <omp.h>
#include <inttypes.h>
#include "argparse.h"
#include "tasks.h"

#define OUTPUT_FILE_PREFIX "./result"
#define OUTPUT_FILE_POSTFIX ".txt"
#define OUTPUT_FILE_MAX_SIZE 1000
#define DEBUG false
#define error(...) fprintf(stderr, __VA_ARGS__)
#define debug_print(...) if (DEBUG) printf(__VA_ARGS__)

#define min(a, b) a < b ? a : b

enum algo {
    Sequential = 1,
    ParallelString = 3,
    WaveChunk = 6
};

typedef struct {
    task_t *task;
    //size of side of the f function matrix
    uint32_t n;
    //size of side of WaveChunk algo's chunk
    uint32_t chunk_size;
    float epsilon;
    //interval between values
    double h;
    struct borders {
        float left;
        float right;
    } borders;
    enum algo algo;
    int random_seed;
    int random_max;
    int random_min;
    bool time_only;
} params_t;

#endif //POISSON_IT_MATH_1_H
