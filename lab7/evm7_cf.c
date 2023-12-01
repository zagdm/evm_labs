#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <xmmintrin.h>

#define N 2048
#define M 10

float* mem_matrix() {
    float* matrix = (float*)malloc(N * N * sizeof(float));
    if (!matrix) exit(1);
    return matrix;
}

void identity_matrix(float* I) {
    for (int i = 0; i < N; i++) {
        float* p = I + i * N;
        for (int j = 0; j < N; j++) p[j] = (float)(i == j);
    }
}

float scalar(float* x, float* y) {
    __m128* xx, * yy;
    __m128 p, s;
    xx = (__m128*) x;
    yy = (__m128*) y;
    s = _mm_setzero_ps();
    for (int i = 0; i < N / 4; ++i) {
        p = _mm_mul_ps(xx[i], yy[i]);
        s = _mm_add_ps(s, p);
    }
    p = _mm_movehl_ps(p, s);
    s = _mm_add_ps(s, p);
    p = _mm_shuffle_ps(s, s, 1);
    s = _mm_add_ss(s, p);
    float sum;
    _mm_store_ss(&sum, s);
    return sum;
}

void mod_transpose_matrix(float* matrix, float* result, float ratio) {
    for (int i = 0; i < N; i++) {
        float* in = matrix + N * i;
        float* out = result + i;
        for (int j = 0; j < N; j++) out[N * j] = in[j] * ratio;
    }
}

void mult_matrix(float* A, float* B, float* res) {
    float* trans = mem_matrix();
    mod_transpose_matrix(B, trans, 1.0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            res[i * N + j] = scalar(&A[i * N], &trans[j * N]);
    free(trans);
}

void summ_matrix(float* A, float* B, float* res, int sign) {
    __m128 p;
    for (int i = 0; i < N; i++) {
        __m128* xx, * yy;
        xx = (__m128*) & A[i * N];
        yy = (__m128*) & B[i * N];
        for (int k = 0; k < N / 4; k++) {
            if (sign) p = _mm_add_ps(xx[k], yy[k]);
            else p = _mm_sub_ps(xx[k], yy[k]);
            _mm_store_ps(&res[i * N + k * 4], p);
        }
    }
}

void copy_matrix(float* inp, float* out) {
    for (int i = 0; i < N; i++) {
        float* a = inp + i * N;
        float* b = out + i * N;
        for (int j = 0; j < N; j++) b[j] = a[j];
    }
}

float circulation_ratio(float* matrix) {
    float max_0 = 0;
    float max_inf = 0;
    for (int i = 0; i < N; i++) {
        float sum_0 = 0;
        float sum_inf = 0;
        float* p_0 = matrix + N * i;
        float* p_inf = matrix + i;
        for (int j = 0; j < N; j++) {
            sum_0 += fabsf(p_0[j]);
            sum_inf += fabsf(p_inf[j * N]);
        }
        if (sum_inf > max_inf) max_inf = sum_inf;
        if (sum_0 > max_0) max_0 = sum_0;
    }
    return 1 / (max_0 * max_inf);
}

double matrix_inverse_time(float* matrix, float* result) {
    clock_t c_start = clock();
    float ratio = circulation_ratio(matrix);
    float* B = mem_matrix();
    mod_transpose_matrix(matrix, B, ratio);
    float* tmp = mem_matrix();
    mult_matrix(B, matrix, tmp);
    float* I = mem_matrix();
    identity_matrix(I);
    float* R = mem_matrix();
    summ_matrix(I, tmp, R, 0);
    copy_matrix(I, tmp);
    for (int k = 0; k < M; ++k) {
        mult_matrix(tmp, R, result);
        summ_matrix(result, I, tmp, 1);
    }
    mult_matrix(tmp, B, result);
    clock_t c_end = clock();
    free(B);
    free(tmp);
    free(I);
    free(R);
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