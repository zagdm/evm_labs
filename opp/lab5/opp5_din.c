#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <pthread.h>
#include <unistd.h>

#define Max_send_count 100
#define Min_task_count 3
#define Total_task_count 1440
#define Total_weight_unit 2000000
#define Weight_unit 20
#define Global_iteration 10

#define task_request_tag 1
#define task_pull_tag 2
#define task_count_tag 3

typedef enum pull_state {
    wait, stop, ready, finish
} pull_state;

typedef struct task_data {
    int task_weight;
} task_data;

typedef struct task_pull{
    pthread_mutex_t mutex;
    pthread_cond_t wait;
    task_data* task_list;
    int list_size;
    pull_state state;
    int worked_weight;
} task_pull;

typedef struct proc_data {
    int size;
    int rank;
    MPI_Datatype task_type;
    task_pull pull;
} proc_data;

void simple_work(task_data* task) {
    usleep(task->task_weight);
}

void* task_worker_thread(void* proc_ptr) {
    proc_data* proc = (proc_data*)proc_ptr;
    task_pull* pull = &(proc->pull);
    int worked = 0;
    while (1) {
        task_data task;
        pthread_mutex_lock(&(pull->mutex));
        if (pull->state == finish) {
            pthread_mutex_unlock(&(pull->mutex));
            break;
        }
        worked = (pull->list_size != 0);
        if (!worked) {
            if (pull->state != ready) pthread_cond_signal(&(pull->wait));
            pull->state = wait;
            pthread_cond_wait(&(pull->wait), &(pull->mutex));
        }
        else {
            pull->list_size--;
            if ((pull->list_size < Min_task_count) && (pull->state == wait)) {
                pull->state = ready;
                pthread_cond_signal(&(pull->wait));
            }
            task = pull->task_list[pull->list_size];
            pull->worked_weight += task.task_weight;
        }
        pthread_mutex_unlock(&(pull->mutex));
        if (worked) simple_work(&task);
    }
    return NULL;
}

int get_sent_count(int size) {
    if (size <= Min_task_count) return 0;
    if (size / 3 > Max_send_count) return Max_send_count;
    return size / 3;
}

void move_task_data(task_data* inp, task_data* out, int count) {
    for (int i = 0; i < count; i++) {
        out[i] = inp[i];
    }
}

void* task_sender_thread(void* proc_ptr) {
    proc_data* proc = (proc_data*)proc_ptr;
    task_pull* pull = &(proc->pull);
    int request_rank;
    task_data tmp_send_list[Max_send_count];
    MPI_Request request;
    MPI_Recv_init(&request_rank, 1, MPI_INT, MPI_ANY_SOURCE, task_request_tag, MPI_COMM_WORLD, &request);
    while (1) {
        MPI_Start(&request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        if (request_rank == proc->rank) break;
        pthread_mutex_lock(&(pull->mutex));
        int send_count = get_sent_count(pull->list_size);
        pull->list_size -= send_count;
        move_task_data(pull->task_list + pull->list_size, tmp_send_list, send_count);
        pthread_mutex_unlock(&(pull->mutex));
        MPI_Send(&send_count, 1, MPI_INT, request_rank, task_count_tag, MPI_COMM_WORLD);
        if (send_count != 0) MPI_Send(tmp_send_list, send_count, proc->task_type, request_rank, task_pull_tag, MPI_COMM_WORLD);
    }
    MPI_Request_free(&request);
    return NULL;
}

void task_receiver_thread(proc_data* proc) {
    task_pull* pull = &(proc->pull);
    int empty_process = 0;
    int current_rank = proc->rank;
    task_data tmp_receive_list[Max_send_count];
    while (empty_process != proc->size) {
        pthread_mutex_lock(&(pull->mutex));
        if (pull->list_size >= Min_task_count) {
            pull->state = wait;
            pthread_cond_wait(&(pull->wait), &(pull->mutex));
        }
        pthread_mutex_unlock(&(pull->mutex));
        int receive_count = 0;
        if (current_rank != proc->rank) {
            MPI_Send(&(proc->rank), 1, MPI_INT, current_rank, task_request_tag, MPI_COMM_WORLD);
            MPI_Recv(&receive_count, 1, MPI_INT, current_rank, task_count_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (receive_count == 0) {
            empty_process++;
            current_rank = (current_rank + 1) % proc->size;
        }
        else {
            empty_process = 0;
            MPI_Recv(tmp_receive_list, receive_count, proc->task_type, current_rank, task_pull_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            pthread_mutex_lock(&(pull->mutex));
            move_task_data(tmp_receive_list, pull->task_list + pull->list_size, receive_count);
            pull->list_size += receive_count;
            if (pull->state == wait) {
                pull->state = ready;
                pthread_cond_signal(&(pull->wait));
            }
            pthread_mutex_unlock(&(pull->mutex));
        }
    }
}

int init_proc_data(proc_data* proc) {
    MPI_Comm_size(MPI_COMM_WORLD, &(proc->size));
    MPI_Comm_rank(MPI_COMM_WORLD, &(proc->rank));
    proc->pull.state = ready;
    proc->pull.list_size = 0;
    proc->pull.task_list = (task_data*)malloc(((Total_task_count / proc->size) + 1) * sizeof(task_data));
    if (proc->pull.task_list == NULL) return -1;
    proc->pull.worked_weight= 0;
    int td_block_count = 1;
    int td_int_count = 1;
    long int td_int_offset = 0;
    MPI_Datatype int_type = MPI_INT;
    MPI_Type_create_struct(td_block_count, &td_int_count, &td_int_offset, &int_type, &(proc->task_type));
    MPI_Type_commit(&(proc->task_type));
    pthread_mutex_init(&(proc->pull.mutex), NULL);
    pthread_cond_init(&(proc->pull.wait), NULL);
    return 0;
}

void free_proc_data(proc_data* proc) {
    free(proc->pull.task_list);
    MPI_Type_free(&(proc->task_type));
    pthread_mutex_destroy(&(proc->pull.mutex));
    pthread_cond_destroy(&(proc->pull.wait));
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
    pthread_mutex_lock(&(proc->pull.mutex));
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
    pthread_mutex_unlock(&(proc->pull.mutex));
}

int main(int argc, char** argv) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        fprintf(stderr, "Invalid thread support level\n");
        MPI_Finalize();
        return 1;
    }

    proc_data proc;
    if (init_proc_data(&proc) == -1) {
        fprintf(stderr, "Invalid allocating memory\n");
        MPI_Finalize();
        return 1;
    }

    pthread_attr_t attribut;
    pthread_attr_init(&attribut);
    pthread_attr_setdetachstate(&attribut, PTHREAD_CREATE_JOINABLE);

    pthread_t worker_thread, sender_thread;
    pthread_create(&worker_thread, &attribut, &task_worker_thread, &proc);
    pthread_create(&sender_thread, &attribut, &task_sender_thread, &proc);
    pthread_attr_destroy(&attribut);

    pthread_mutex_lock(&(proc.pull.mutex));
    if (proc.pull.state != wait) {
        proc.pull.state = stop;
        pthread_cond_wait(&(proc.pull.wait), &(proc.pull.mutex));
    }
    pthread_mutex_unlock(&(proc.pull.mutex));

    MPI_Barrier(MPI_COMM_WORLD);
    double global_start = MPI_Wtime();
    for (int i = 0; i < Global_iteration; i++) {
        fill_task_pull(&proc, i);
        MPI_Barrier(MPI_COMM_WORLD);
        double local_start = MPI_Wtime();

        pthread_mutex_lock(&(proc.pull.mutex));
        proc.pull.state = ready;
        pthread_cond_signal(&(proc.pull.wait));
        pthread_mutex_unlock(&(proc.pull.mutex));

        task_receiver_thread(&proc);
        double local_end = MPI_Wtime();

        pthread_mutex_lock(&(proc.pull.mutex));
        if (proc.pull.state != wait) {
            proc.pull.state = stop;
            pthread_cond_wait(&(proc.pull.wait), &(proc.pull.mutex));
        }
        int worked = proc.pull.worked_weight;
        pthread_mutex_unlock(&(proc.pull.mutex));
        fprintf(stdout, "Local work time: %lf, worked %d, iter %d, rank %d\n", local_end - local_start, worked, i + 1, proc.rank);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    double global_end = MPI_Wtime();
    if (proc.rank == 0) {
        fprintf(stdout, "Total time: %lf\n", global_end - global_start);
    }

    pthread_mutex_lock(&(proc.pull.mutex));
    if (proc.pull.state == wait) pthread_cond_signal(&(proc.pull.wait));
    proc.pull.state = finish;
    pthread_mutex_unlock(&(proc.pull.mutex));
    pthread_join(worker_thread, NULL);

    MPI_Send(&(proc.rank), 1, MPI_INT, proc.rank, task_request_tag, MPI_COMM_WORLD);
    pthread_join(sender_thread, NULL);

    free_proc_data(&proc);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    
    return 0;
}