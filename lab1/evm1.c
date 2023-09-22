#include <stdio.h>

#define terms_count 3180000000

double calc_pi(long long int n) {
	double pi_4 = 0;
	for (long long int i = 0; i <= n; i++) {
		pi_4 += ((i&1)? -1 : 1) * (((double)1) / (2 * i + 1));
	}
	return 4 * pi_4;
}

int main() {
	double pi = calc_pi(terms_count);
	printf("%10.9f", pi);
	return 0;
}
