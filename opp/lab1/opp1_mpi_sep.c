#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

const double eps = 0.00001;
const int N = 46000;

void v_matrix_inc_mult(const double* matrix, const double* vector, double* result, int v_size, int m_size, int m_rank) {
    for (int i = 0; i < m_size; i++) {
        for (int j = 0; j < v_size; j++) {
            result[i] += matrix[i * m_rank + j] * vector[j];
        }
    }
}

void v_matrix_mult(const double* matrix, const double* vector, double* result, int v_size, int m_size) {
    for (int i = 0; i < m_size; i++) {
        result[i] = 0;
        for (int j = 0; j < v_size; j++) {
            result[i] += matrix[i * v_size + j] * vector[j];
        }
    }
}

double v_scalar_mult(const double* first, const double* second, int size) {
    double res = 0;
    for (int i = 0; i < size; ++i) res += first[i] * second[i];
    return res;
}

void v_subtract(double* vector, const double* sub_vector, int size) {
    for (int i = 0; i < size; i++) vector[i] -= sub_vector[i];
}

void v_linear_comb(double* v_first, double* v_second, double* result, int size, double c_sec) {
    for (int i = 0; i < size; i++) result[i] = v_first[i] + c_sec * v_second[i];
}

int main_exit(void** mem) {
    for (int i = 0; i < 8; i++) free(mem[i]);
    return 0;
}

int main(int argc, char** argv) {
    int error = MPI_Init(&argc, &argv);
    if (error) return error;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > N) {
        fprintf(stderr, "To many process\n");
        MPI_Finalize();
        return 0;
    }
    
    double start = MPI_Wtime();

    void* mem[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    double* y = (double*)(mem[0] = malloc(N * sizeof(double)));
    int* part_table = (int*)(mem[1] = malloc(size * sizeof(int)));
    int* shift_table = (int*)(mem[2] = malloc(size * sizeof(int)));
    if (!(part_table && shift_table && y)) {
        fprintf(stderr, "Memory error\n");
        MPI_Finalize();
        return main_exit(mem);
    }
    for (int i = 0; i < size; i++) part_table[i] = (N / size) + (i < (N% size));
    for (int i = 0; i < size; i++) shift_table[i] = (N / size) * i + min(i, N % size);
    int part_size = part_table[rank];
    int part_shift = shift_table[rank];

    double* part_A = (double*)(mem[3] = malloc(N * part_size * sizeof(double)));
    double* part_x = (double*)(mem[4] = calloc(N, sizeof(double)));
    double* part_y = (double*)(mem[5] = malloc(part_size * sizeof(double)));
    double* tmp = (double*)(mem[6] = malloc(part_table[0] * sizeof(double)));
    double* b = (double*)(mem[7] = malloc(part_size * sizeof(double)));
    if (!(part_A && part_y && tmp && part_x)) {
        fprintf(stderr, "Memory error\n");
        MPI_Finalize();
        return main_exit(mem);
    }

    for (int i = 0; i < N * part_size; i++) A[i] = 1;
    for (int i = 0; i < part_size; i++) b[i] = N + 1;
    for (int i = 0; i < part_size; ++i) part_A[i * N + i + part_shift] += 1;
 
    int iteration = 0;
    double part_b_scalar = v_scalar_mult(b, b, part_size);
    double b_scalar = 0.0;
    MPI_Allreduce(&part_b_scalar, &b_scalar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double crit = b_scalar * eps * eps;

    while (++iteration) {
        for (int i = 0; i < part_size; i++) part_y[i] = 0;
        for (size_t i = 0; i < size; ++i) {
            int load_rank = (rank + i) % size;
            MPI_Request request[2];
            MPI_Status status[2];
            MPI_Isend(part_x, part_size, MPI_DOUBLE, (rank - i + size) % size, 42, MPI_COMM_WORLD, &request[0]);
            MPI_Irecv(tmp, part_table[load_rank], MPI_DOUBLE, load_rank, 42, MPI_COMM_WORLD, &request[1]);
            MPI_Waitall(2, request, status);
            v_matrix_inc_mult(part_A + shift_table[load_rank], tmp, part_y, part_table[load_rank], part_size, N);
        }
        v_subtract(part_y, b, part_size);

        double part_y_scalar = v_scalar_mult(part_y, part_y, part_size);
        double y_scalar = 0.0;
        MPI_Allreduce(&part_y_scalar, &y_scalar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (y_scalar < crit) break;
       
        MPI_Allgatherv(part_y, part_size, MPI_DOUBLE, y, part_table, shift_table, MPI_DOUBLE, MPI_COMM_WORLD);
        v_matrix_mult(part_A, y, tmp, N, part_size);

        double part_y_Ay_scalar = v_scalar_mult(part_y, tmp, part_size);
        double y_Ay_scalar = 0.0;
        MPI_Allreduce(&part_y_Ay_scalar, &y_Ay_scalar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double part_Ay_Ay_scalar = v_scalar_mult(tmp, tmp, part_size);
        double Ay_Ay_scalar = 0.0;
        MPI_Allreduce(&part_Ay_Ay_scalar, &Ay_Ay_scalar, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (Ay_Ay_scalar == 0) Ay_Ay_scalar = 1;
        double t = (-1) * y_Ay_scalar / Ay_Ay_scalar;
        v_linear_comb(part_x, part_y, part_x, part_size, t);
    }

    double end = MPI_Wtime();

    MPI_Finalize();
    if (rank == 0) {
        fprintf(stdout, "Iteration count: %d\n", iteration);
        fprintf(stdout, "Time passed: %lf\n", end - start);
    }
    return main_exit(mem);
}
