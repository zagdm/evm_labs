#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

#define Total_task_count 1440
#define Total_weight_unit 2000000
#define Weight_unit 20
#define Global_iteration 10

typedef struct task_data {
    int task_weight;
} task_data;

typedef struct task_pull{
    task_data* task_list;
    int list_size;
    int worked_weight;
} task_pull;

typedef struct proc_data {
    int size;
    int rank;
    task_pull pull;
} proc_data;

void simple_work(task_data* task) {
    usleep(task->task_weight);
}


int init_proc_data(proc_data* proc) {
    MPI_Comm_size(MPI_COMM_WORLD, &(proc->size));
    MPI_Comm_rank(MPI_COMM_WORLD, &(proc->rank));
    proc->pull.list_size = 0;
    proc->pull.task_list = (task_data*)malloc(((Total_task_count / proc->size) + 1) * sizeof(task_data));
    if (proc->pull.task_list == NULL) return -1;
    proc->pull.worked_weight = 0;
    return 0;
}

int distribution(double mean_value, double max_dist) {
    double dist = mean_value - max_dist;
    if (dist < 0) dist *= -1;
    if (dist > 0.5) dist = 1 - dist;
    dist *= 2;
    double norm = 2;
    norm = norm * dist - 2.5;
    norm = norm * dist - 0.5;
    norm = norm * dist + 1.5;
    return (int)(norm * 1500);
}

void fill_task_pull(proc_data* proc, int seed) {
    proc->pull.worked_weight = 0;
    int min_count = Total_task_count % proc->size;
    int pull_size = Total_task_count / proc->size + (proc->rank < min_count);
    int pull_shift = (Total_task_count / proc->size) * proc->rank + ((proc->rank < min_count) ? proc->rank : min_count);
    for (int i = 0; i < pull_size; i++) {
        proc->pull.task_list[i].task_weight = 0;
    }
    double mead_value = (double)seed / Global_iteration;
    int weight_left = Total_weight_unit;
    for (int it = 0; (it < Total_task_count) && (weight_left != 0); it++){
        int requested_weight = distribution(mead_value, (double)it / Total_task_count);
        int i = it - pull_shift;
        if (requested_weight > weight_left) requested_weight = weight_left;
        if ((i >= 0) && (i < pull_size)) proc->pull.task_list[i].task_weight += Weight_unit * requested_weight;
        weight_left -= requested_weight;
    }
    srand(seed);
    while (weight_left != 0) {
        int rand_int = rand();
        int requested_weight = (rand_int % 1000) + 1;
        if (requested_weight > weight_left) requested_weight = weight_left;
        int i = (rand_int % Total_task_count) - pull_shift;
        if ((i >= 0) && (i < pull_size)) proc->pull.task_list[i].task_weight += Weight_unit * requested_weight;
        weight_left -= requested_weight;
    }
    proc->pull.list_size = pull_size;
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    proc_data proc;
    if (init_proc_data(&proc) == -1) {
        fprintf(stderr, "Invalid allocating memory\n");
        MPI_Finalize();
        return 1;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double global_start = MPI_Wtime();
    for (int i = 0; i < Global_iteration; i++) {
        fill_task_pull(&proc, i);
        MPI_Barrier(MPI_COMM_WORLD);
        double local_start = MPI_Wtime();

        while (proc.pull.list_size != 0) {
            proc.pull.list_size--;
            task_data task = proc.pull.task_list[proc.pull.list_size];
            proc.pull.worked_weight += task.task_weight;
            simple_work(&task);
        }
        double local_end = MPI_Wtime();

        int worked = proc.pull.worked_weight;
        fprintf(stdout, "Local work time: %lf, worked %d, iter %d, rank %d\n", local_end - local_start, worked, i + 1, proc.rank);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    double global_end = MPI_Wtime();
    if (proc.rank == 0) {
        fprintf(stdout, "Total time: %lf\n", global_end - global_start);
    }

    free(proc.pull.task_list);
    MPI_Finalize();

    return 0;
}
