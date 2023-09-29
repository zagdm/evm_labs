#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define input_error 1
#define step 50000000
#define init_quant 550000000

double e_sin(double x) {
	return exp(x) * sin(x);
}

double calc_integral(double f(double), int N, double a, double b) {
	double m_riem_sum = 0;
	const double dx = (double)(b - a) / N;
	for (int i = 0; i < N; i++) m_riem_sum += f(a + i * dx) + f(a + (i + 1) * dx);
	return m_riem_sum * dx / 2;
}

int main(int argc, char* argv[]) {
	if (argc < 2) return input_error;
	double integral = calc_integral(e_sin, init_quant + step * atoi(argv[1]), 0, M_PI);
	printf("%20.19f", integral);
	return 0;
}
