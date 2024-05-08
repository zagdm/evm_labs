#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define Min_pos (-1)
#define Max_pos 1

const double eps_param = 1e-8;
const double a_param = 1e5;
const int discret_param = 256;

typedef struct process {
    int rank;
    int size;
} process;

typedef struct proc_mem {
    double** layers;
    double* z_up;
    double* z_down;
    double* tmp_res;
    double* tmp_last;
    double* tmp_up;
    double* tmp_down;
    void* mem;
} proc_mem;

void swap_p(double** a, double** b) {
    double* tmp = *a;
    *a = *b;
    *b = tmp;
}

double phi_func(double x, double y, double z) {
    return x * x + y * y + z * z;
}

double ro(double x, double y, double z) {
    return 6 - a_param * phi_func(x, y, z);
}

double max_d(double a, double b) {
    if (a < b) {
        return b;
    }
    return a;
}

double abs_d(double a) {
    if (a > 0) {
        return a;
    }
    return -a;
}

double calc_step_c(double step_x, double step_y, double step_z) {
    double step_c = 0;
    step_c += 2 / (step_x * step_x);
    step_c += 2 / (step_y * step_y);
    step_c += 2 / (step_z * step_z);
    step_c += a_param;
    step_c = 1 / step_c;
    return step_c;
}

double calculate_layer(double* layer, double z, int x_size, int y_size, double dx, double dy, double dz, double* result, double* d_layer, double* u_layer) {
    static double step_c = 0;
    if (step_c == 0) step_c = calc_step_c(dx, dy, dz);
    double max_diff = 0;
    for (int y = 0; y < y_size; y++) {
        for (int x = 0; x < x_size; x++) {
            if ((x % (x_size - 1)) && (y % (y_size - 1))) {
                double step_phi = (layer[y * x_size + x + 1] + layer[y * x_size + x - 1]) / (dx * dx);
                step_phi += (layer[(y - 1) * x_size + x] + layer[(y - 1) * x_size + x]) / (dy * dy);
                step_phi += (d_layer[y * x_size] + u_layer[y * x_size + x]) / (dz * dz);
                step_phi -= ro(Min_pos + x * dx, Min_pos + y * dy, z);
                step_phi *= step_c;
                max_diff = max_d(max_diff, abs_d(step_phi - layer[y * x_size + x]));
                result[y * x_size + x] = step_phi;
            }
            else {
                result[y * x_size + x] = layer[y * x_size + x];
            }
        }
    }
    return max_diff;
}

int init_proc_mem(proc_mem* memory, int x_size, int y_size, int z_size, const process* proc) {
    size_t layer_size = x_size * y_size;
    size_t memory_size = layer_size * z_size + z_size + layer_size;
    double* last_ptr;
    if (proc->rank != 0) {
        memory_size += 2 * layer_size;
    }
    if (proc->rank != proc->size - 1) {
        memory_size += 3 * layer_size;
    }
    last_ptr = memory->mem = calloc(memory_size, sizeof(double));
    if (memory->mem == NULL) return 1;
    memory->layers = (double**)last_ptr;
    last_ptr += z_size;
    for (int i = 0; i < z_size; i++) {
        memory->layers[i] = last_ptr;
        last_ptr += layer_size;
    }
    memory->tmp_res = last_ptr;
    last_ptr += layer_size;
    if (proc->rank != 0) {
        memory->z_down = last_ptr;
        last_ptr += layer_size;
        memory->tmp_down = last_ptr;
        last_ptr += layer_size;
    }
    else {
        memory->tmp_down = memory->z_down = NULL;
    }
    if (proc->rank != proc->size - 1) {
        memory->z_up = last_ptr;
        last_ptr += layer_size;
        memory->tmp_up = last_ptr;
        last_ptr += layer_size;
        memory->tmp_last = last_ptr;
    }
    else {
        memory->tmp_up = memory->z_up = memory->tmp_last = NULL;
    }
    return 0;
}

void init_bound_layer(double* layer, int x_size, int y_size, double dx, double dy, double z, double z_bound) {
    if (z_bound) {
        for (int i = 0; i < x_size * y_size; i++) {
            layer[i] = phi_func(Min_pos + dx * (i % x_size), Min_pos + dy * (i / x_size), z);
        }
    }
    else {
        for (int i = 0; i < x_size; i++) {
            layer[i] = phi_func(Min_pos + dx * i, Min_pos, z);
            layer[i + x_size * (y_size - 1)] = phi_func(Min_pos + dx * i, Max_pos, z);
        }
        for (int i = 1; i < y_size - 1; i++) {
            layer[i * x_size] = phi_func(Min_pos, Min_pos + dy * i ,z);
            layer[x_size - 1 + i * x_size] = phi_func(Max_pos, Min_pos + dy * i, z);
        }
    }
}

double phi_different(double** layers, int x_size, int y_size, int z_size, double dx, double dy, double dz, double z_min) {
    double max_phi_diff = 0;
    for (int k = 0; k < z_size; k++) {
        for (int i = 0; i < y_size * x_size; i++) {
            double tmp = phi_func(Min_pos + dx * (i % x_size), Min_pos + dy * (i / x_size), z_min + dz * k);
            tmp = abs_d(layers[k][i] - tmp);
            max_phi_diff = max_d(max_phi_diff, tmp);
        }
    }
    return max_phi_diff;
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    process proc;
    MPI_Comm_size(MPI_COMM_WORLD, &(proc.size));
    MPI_Comm_rank(MPI_COMM_WORLD, &(proc.rank));

    if ((discret_param % proc.size) || ((discret_param / proc.size) < 3)) {
        fprintf(stderr, "Invalid number of processes\n");
        MPI_Finalize();
        return 4;
    };

    double start = MPI_Wtime();
    int x_size, y_size, z_size;
    x_size = y_size = z_size = discret_param;
    z_size /= proc.size;

    double dx, dy, dz;
    dx = dy = dz = ((double)Max_pos - Min_pos) / (discret_param - 1);
    double z_min = Min_pos + dz * proc.rank * z_size;

    proc_mem phi;
    int error = init_proc_mem(&phi, x_size, y_size, z_size, &proc);
    if (error) {
        fprintf(stderr, "Memory error %d\n", proc.rank);
        MPI_Finalize();
        return error;
    }
    if (proc.rank != 0) {
        init_bound_layer(phi.z_down, x_size, y_size, dx, dy, z_min - dz, 0);
    }
    if (proc.rank != proc.size - 1) {
        init_bound_layer(phi.z_up, x_size, y_size, dx, dy, z_min + dz * z_size, 0);
    }
    init_bound_layer(phi.layers[0], x_size, y_size, dx, dy, z_min, (proc.rank == 0));
    init_bound_layer(phi.layers[z_size - 1], x_size, y_size, dx, dy, z_min + dz * (z_size - 1), (proc.rank == proc.size - 1));
    for (int i = 1; i < z_size - 1; i++) {
        init_bound_layer(phi.layers[i], x_size, y_size, dx, dy, z_min + dz * i, 0);
    }

    double max_diff;
    double tmp_max = 0;
    do {
        max_diff = 0;
        MPI_Request rq[4];
        if (proc.rank != proc.size - 1) {
            MPI_Irecv(phi.tmp_up, x_size * y_size, MPI_DOUBLE, proc.rank + 1, 80, MPI_COMM_WORLD, &rq[2]);
            tmp_max = calculate_layer(phi.layers[z_size - 1], z_min + dz * (z_size - 1), x_size, y_size, dx, dy, dz, phi.tmp_last, phi.layers[z_size - 2], phi.z_up);
            MPI_Isend(phi.tmp_last, x_size * y_size, MPI_DOUBLE, proc.rank + 1, 80, MPI_COMM_WORLD, &rq[3]);
            max_diff = max_d(max_diff, tmp_max);
        }
        if (proc.rank != 0) {
            MPI_Irecv(phi.tmp_down, x_size * y_size, MPI_DOUBLE, proc.rank - 1, 80, MPI_COMM_WORLD, &rq[0]);
            tmp_max = calculate_layer(phi.layers[0], z_min, x_size, y_size, dx, dy, dz,  phi.tmp_res, phi.z_down, phi.layers[1]);
            MPI_Isend(phi.tmp_res, x_size * y_size, MPI_DOUBLE, proc.rank - 1, 80, MPI_COMM_WORLD, &rq[1]); 
        }
        else {
            tmp_max = calculate_layer(phi.layers[1], z_min + dz, x_size, y_size, dx, dy, dz, phi.tmp_res, phi.layers[0], phi.layers[2]);
        }
        max_diff = max_d(max_diff, tmp_max);
        swap_p(&(phi.tmp_res), &(phi.layers[(proc.rank == 0)]));
        for (int i = 1 + (proc.rank == 0); i < z_size - 1; i++) {
            tmp_max = calculate_layer(phi.layers[i], z_min + dz * i, x_size, y_size, dx, dy, dz, phi.tmp_res, phi.tmp_res, phi.layers[i + 1]);
            max_diff = max_d(max_diff, tmp_max);
            swap_p(&(phi.tmp_res), &(phi.layers[i]));
        }
        if (proc.rank != 0) {
            MPI_Wait(&rq[0], MPI_STATUS_IGNORE);
            MPI_Wait(&rq[1], MPI_STATUS_IGNORE);
            swap_p(&(phi.tmp_down), &(phi.z_down));
        }
        if (proc.rank != proc.size - 1) {
            swap_p(&(phi.tmp_last), &(phi.layers[z_size - 1]));
            MPI_Wait(&rq[2], MPI_STATUS_IGNORE);
            MPI_Wait(&rq[3], MPI_STATUS_IGNORE);
            swap_p(&(phi.tmp_up), &(phi.z_up));
        }
        MPI_Allreduce(&max_diff, &tmp_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    } while (tmp_max >= eps_param);

    double end = MPI_Wtime();

    tmp_max = phi_different(phi.layers, x_size, y_size, z_size, dx, dy, dz, z_min);
    MPI_Allreduce(&tmp_max, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    free(phi.mem);
    MPI_Finalize();

    if (proc.rank == 0) {
        fprintf(stdout, "Time passed: %lf\n", end - start);
        fprintf(stdout, "Maximud phi different %lf\n", max_diff);
    }

    return 0;
}
