#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

const double eps = 0.0001;
const int N = 23700;

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
    for (int i = 0; i < 5; i++) free(mem[i]);
    return 0;
}

int main() {
    auto start = clock();

    void* mem[5] = { 0, 0, 0, 0, 0 };
    double* x = (double*)(mem[0] = calloc(N, sizeof(double)));
    double* b = (double*)(mem[1] = malloc(N * sizeof(double)));
    double* A = (double*)(mem[2] = malloc((size_t)N * N * sizeof(double)));
    double* y = (double*)(mem[3] = malloc(N * sizeof(double)));
    double* tmp = (double*)(mem[4] = malloc(N * sizeof(double)));
    if (!(x && b && A && y && tmp)) {
        fprintf(stderr, "Memory error");
        return main_exit(mem);
    }

    memset(A, 1, N * N);
    memset(b, N + 1, N);
    for (int i = 0; i < N; ++i) A[i * N + i] += 1;

    int iteration = 0;
    double criterial = v_scalar_mult(b, b, N) * eps * eps;
    while (++iteration) {
        v_matrix_mult(A, x, y, N, N);
        v_subtract(y, b, N);
        if (v_scalar_mult(y, y, N) < criterial) break;
        v_matrix_mult(A, y, tmp, N, N);
        double t = v_scalar_mult(tmp, tmp, N);
        t = v_scalar_mult(y, tmp, N) / ((t == 0) ? 1 : t);
        v_linear_comb(x, y, x, N, t * (-1));
    }

    auto end = clock();

    fprintf(stdout, "Iteration count: %d\n", iteration);
    fprintf(stdout, "Time passed: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
    return main_exit(mem);
}
