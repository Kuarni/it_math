#include "tasks.h"
#include "task_1/task_1.c"

//add new task here
const task_t tasks[] = {task_1};

int find_task(const char *name, task_t **task_pointer) {
    for (int i = 0; i < sizeof(tasks); i++) {
        if (!strcmp(name, tasks[i].name)) {
            *task_pointer = (task_t*) &(tasks[i]);
            return 0;
        }
    }
    return -EINVAL;
}
