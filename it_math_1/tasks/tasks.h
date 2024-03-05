#ifndef POISSON_TASKS_H
#define POISSON_TASKS_H

#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

typedef struct {
    //size = TASK_NAME_SIZE
    char *name;

    //given function
    int (*f)(double **matrix, int n);

    //boundaries for sought function (u)
    double (*bottom)(double x);

    double (*left)(double y);

    double (*upper)(double x);

    double (*right)(double y);
} task_t;

int find_task(const char *name, task_t **task_pointer);

#endif //POISSON_TASKS_H
