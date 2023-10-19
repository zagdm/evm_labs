#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define partition_degree 1000000

double calc_integral(int N, double a, double b) {
	const double dx = (b - a) / N;
	double m_riem_sum = (exp(a)*sin(a) + exp(b)*sin(b)) / 2;
	for (int i = 1; i < N; i++) m_riem_sum += exp(a + i * dx) * sin(a + i * dx);
	return m_riem_sum * dx;
}

int main() {
	double integral = calc_integral(partition_degree, 0, M_PI);
	printf("%f", integral);
	return 0;
}
