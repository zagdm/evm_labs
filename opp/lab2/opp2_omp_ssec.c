#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

const double eps = 0.0001;
const size_t N = 55000;

int main_exit(void** mem) {
    for (int i = 0; i < 5; i++) free(mem[i]);
    return 0;
}

int main() {
    double start = omp_get_wtime();

    void* mem[5] = { 0, 0, 0, 0, 0 };
    double* x = (double*)(mem[0] = calloc(N, sizeof(double)));
    double* b = (double*)(mem[1] = malloc(N * sizeof(double)));
    double* A = (double*)(mem[2] = malloc(N * N * sizeof(double)));
    double* y = (double*)(mem[3] = malloc(N * sizeof(double)));
    double* tmp = (double*)(mem[4] = malloc(N * sizeof(double)));
    if (!(x && b && A && y && tmp)) {
        fprintf(stderr, "Memory error\n");
        return main_exit(mem);
    }
    int iteration = 0;
    double y_scal;
    double t, t_m;
    double criterial = 0;
    #pragma omp parallel
    {
        #pragma omp for
        for (size_t i = 0; i < N * N; i++) A[i] = 1;
        #pragma omp for
        for (size_t i = 0; i < N; i++) b[i] = N + 1;
        #pragma omp for
        for (size_t i = 0; i < N; ++i) A[i * N + i] += 1;

        #pragma omp for reduction(+:criterial)
        for (size_t i = 0; i < N; ++i) criterial += b[i] * b[i];
        #pragma omp single nowait
        criterial = criterial * eps * eps;
        while (1) {
            #pragma omp single nowait 
            {
                iteration++;
                t = 0;
                t_m = 0;
                y_scal = 0;
            }
            #pragma omp for
            for (size_t i = 0; i < N; i++) {
                y[i] = 0;
                for (size_t j = 0; j < N; j++) {
                    y[i] += A[i * N + j] * x[j];
                }
            }
            #pragma omp for
            for (size_t i = 0; i < N; i++) y[i] -= b[i];
            #pragma omp for
            for (size_t i = 0; i < N; i++) {
                tmp[i] = 0;
                for (size_t j = 0; j < N; j++) {
                    tmp[i] += A[i * N + j] * y[j];
                }
            }
            #pragma omp for reduction(+:t)
            for (size_t i = 0; i < N; ++i) t += tmp[i] * tmp[i];
            #pragma omp for reduction(+:t_m)
            for (size_t i = 0; i < N; ++i) t_m += y[i] * tmp[i];
            #pragma omp single
            t = ((t == 0) ? t_m : (t_m / t));
            #pragma omp for
            for (size_t i = 0; i < N; i++) x[i] = x[i] - y[i] * t;
            #pragma omp for reduction(+:y_scal)
            for (size_t i = 0; i < N; ++i) y_scal += y[i] * y[i];
            if (y_scal < criterial) break;
        }
    }
    double end = omp_get_wtime();

    fprintf(stdout, "Iteration count: %d\n", iteration);
    fprintf(stdout, "Time passed: %lf\n", end - start);
    return main_exit(mem);
}
