#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cblas.h>

#define N 2048
#define M 10

float* mem_matrix() {
    float* matrix = (float*)malloc(N * N * sizeof(float));
    if (!matrix) exit(1);
    return matrix;
}

void identity_matrix(float* Ed) {
    for (int i = 0; i < N; i++) {
        float* p = Ed + i * N;
        for (int j = 0; j < N; j++) p[j] = (float)(i == j);
    }
}

float circulation_ratio(float* matrix) {
    float max_0 = 0;
    float max_inf = 0;
    for (int i = 0; i < N; i++) {
        float sum_0 = cblas_sasum(N, matrix + N * i, 1);
        float sum_inf = cblas_sasum(N, matrix + i, N);
        if (sum_0 > max_0) max_0 = sum_0;
        if (sum_inf > max_inf) max_inf = sum_inf;
    }
    return 1 / (max_0 * max_inf);
}

double matrix_inverse_time(float* matrix, float* result) {
    clock_t c_start = clock();
    float* Ed = mem_matrix();
    float* R = mem_matrix();
    float* B = mem_matrix();
    float* tmp = mem_matrix();
    identity_matrix(Ed);
    cblas_scopy(N * N, Ed, 1, R, 1);
    cblas_scopy(N * N, Ed, 1, tmp, 1);
    float ratio = circulation_ratio(matrix);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, N, N, N, ratio, Ed, N, matrix, N, 0.0, B, N);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, -1.0, matrix, N, B, N, 1.0, R, N);
    for (int k = 0; k < M; ++k) {
        cblas_scopy(N * N, Ed , 1, k & 1 ? tmp : result, 1);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, k & 1 ? result : tmp, N, R, N, 1.0, k & 1 ? tmp : result, N);
    }
    if (M & 1) cblas_scopy(N * N, result, 1,tmp, 1);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, tmp, N, B, N, 0.0, result, N);
    clock_t c_end = clock();
    free(B);
    free(tmp);
    free(R);
    free(Ed);
    return (double)(c_end - c_start) / CLOCKS_PER_SEC;;
}

int main() {
    float* matrix = mem_matrix();
    float* backMatrix = mem_matrix();
    for (int k = 0; k < N; ++k)
        for (int i = 0; i < N; ++i) {
            matrix[k * N + i] = (float)(rand() % 10);
            backMatrix[k * N + i] = 0;
        }
    double running_time = matrix_inverse_time(matrix, backMatrix);
    printf("Algorithm running time: %f s", running_time);
    free(matrix);
    free(backMatrix);
    return 0;
}