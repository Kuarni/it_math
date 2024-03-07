#include "../tasks.h"

static int f_task_1(double **matrix, int n) {
    return 0;
}

static double bottom_1(double x, double y) {
    return 100 - 200 * x;
}

static double left_1(double x, double y) {
    return 100 - 200 * y;
}

static double upper_1(double x, double y) {
    return -100 + 200 * x;
}

static double right_1(double x, double y) {
    return -100 + 200 * y;
}

const task_t task_1 = {"task_1", f_task_1, bottom_1, left_1, upper_1, right_1};