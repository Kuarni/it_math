#include "../tasks.h"
#include "math.h"

static double f_task_2(double x, double y) {
    return 1000 * pow(x, 3) + 2000 * pow(y, 3);
}

static double border_2(double x, double y) {
    return 6000 * x + 12000 * y;
}

const task_t task_2 = {"task_2", f_task_2, border_2, border_2, border_2, border_2};