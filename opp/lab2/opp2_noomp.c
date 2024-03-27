#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const double eps = 0.0001;
const int N = 24000;

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
    int start = clock();

    void* mem[5] = { 0, 0, 0, 0, 0 };
    double* x = (double*)(mem[0] = calloc(N, sizeof(double)));
    double* b = (double*)(mem[1] = malloc(N * sizeof(double)));
    double* A = (double*)(mem[2] = malloc((size_t)N * N * sizeof(double)));
    double* y = (double*)(mem[3] = malloc(N * sizeof(double)));
    double* tmp = (double*)(mem[4] = malloc(N * sizeof(double)));
    if (!(x && b && A && y && tmp)) {
        fprintf(stderr, "Memory error\n");
        return main_exit(mem);
    }

    for (int i = 0; i < N * N; i++) A[i] = 1;
    for (int i = 0; i < N; i++) b[i] = N + 1;  
    for (int i = 0; i < N; ++i) A[i * N + i] += 1;

    int iteration = 0;
    double criterial = v_scalar_mult(b, b, N) * eps * eps;
    while (++iteration) {
        v_matrix_mult(A, x, y, N, N);
        v_subtract(y, b, N);
        v_matrix_mult(A, y, tmp, N, N);
        double t = v_scalar_mult(tmp, tmp, N);
        t = v_scalar_mult(y, tmp, N) / ((t == 0) ? 1 : t);
        v_linear_comb(x, y, x, N, t * (-1));
        if (v_scalar_mult(y, y, N) < criterial) break;
    }

    int end = clock();

    fprintf(stdout, "Iteration count: %d\n", iteration);
    fprintf(stdout, "Time passed: %lf\n", (double)(end - start) / CLOCKS_PER_SEC);
    return main_exit(mem);
}



/*
* N = 5000, eps = 0.0001, b = A(rand)*u(sin)
0.002119, .. , 0.587852, .. , 0.069271, .. , 0.948860, .. , 0.224209, .. , 0.001729, .. , -0.259550, .. , -0.740170, .. , -0.152699, .. , -0.254938, .. ,
0.000000, .. , 0.587785, .. , 0.951057, .. , 0.951057, .. , 0.587785, .. , 0.000000, .. , -0.587785, .. , -0.951057, .. , -0.951057, .. , -0.587785, .. ,
*/

/*
* N = 500, eps = 0.0001, b = A(rand)*u(sin)
-0.008027, .. , 0.208514, .. , 0.945603, .. , 0.801689, .. , 0.466958, .. , -0.007307, .. , -0.588712, .. , -0.252286, .. , -0.458136, .. , -0.571164, .. ,
0.000000, .. , 0.587785, .. , 0.951057, .. , 0.951057, .. , 0.587785, .. , 0.000000, .. , -0.587785, .. , -0.951057, .. , -0.951057, .. , -0.587785, .. ,
*/

/*
* N = 500, eps = 0.00000001, b = A(rand)*u(sin)
-0.072713, .. , 0.586691, .. , 0.950978, .. , 0.950821, .. , 0.587512, .. , -0.001110, .. , -0.587826, .. , -0.952611, .. , -0.951786, .. , -0.587916, .. ,
0.000000, .. , 0.587785, .. , 0.951057, .. , 0.951057, .. , 0.587785, .. , 0.000000, .. , -0.587785, .. , -0.951057, .. , -0.951057, .. , -0.587785, .. ,
*/

/*
* N = 500, eps = 0.0000000000001, b = A(rand)*u(sin)
-0.003859, .. , 0.587774, .. , 0.951056, .. , 0.951054, .. , 0.587782, .. , -0.000011, .. , -0.587786, .. , -0.951073, .. , -0.951064, .. , -0.587787, .. ,
0.000000, .. , 0.587785, .. , 0.951057, .. , 0.951057, .. , 0.587785, .. , 0.000000, .. , -0.587785, .. , -0.951057, .. , -0.951057, .. , -0.587785, .. ,
Iteration count: 272254
*/
