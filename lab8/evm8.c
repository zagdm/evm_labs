#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const long long Nmin = 1024 / sizeof(int);
const long long Nmax = 32 * 1024 * 1024LL / sizeof(int);
const int K = 10;

void gen_ascending(int* data, const long long N) {
    for (long long i = 0; i < N - 1; i++) {
        data[i] = i + 1;
    }
    data[N] = 0;
}

void gen_descending(int* data, const long long N) {
    data[0] = N - 1;
    for (long long i = 1; i < N; i++) {
        data[i] = i - 1;
    }
}

void gen_random(int* data, const long long N) {
    int temp;
    long long i = N;
    srand(time(NULL));
    while (i > 1) {
        i--;
        int j = rand() % i;
        temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }
}
int test(int* data, const long long N, const int K, void (*gen)(int* data, const long long N)) {
    gen(data, N);

    // warm-up
    int k;
    long long i;

    for (k = 0, i = 0; i < N; i++) {
        k = data[k];
    }
    unsigned long long start = __builtin_ia32_rdtsc();
    for (k = 0, i = 0; i < N * K; i++) {
        k = data[k];
    }
    unsigned long long end = __builtin_ia32_rdtsc();
    printf("N = %lld, Tick = %llu\n", N * 4, (end - start) / N / K);
    return k;
}

int main() {
    int* data = (int*)malloc(Nmax * sizeof(int));
    for (int N = Nmin; N <= Nmax; N *= 2) {
        printf("Ascending: ");
        test(data, N, K, gen_ascending);
        printf("Descending: ");
        test(data, N, K, gen_descending);
        printf("SattoloAlgo: ");
        test(data, N, K, gen_random);
    }
    free(data);
    return 0;
}