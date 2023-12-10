#include <iostream>
#include <cmath>
#include <x86intrin.h>

static const size_t KB = 1024;

int set_size(char* n_ascii) {
	return KB * (int)(pow(2, (double)atoi(n_ascii) / 2) + 0.5);
}

class Loop {
private:
	static const size_t round_count = 100;
	int size;
	int* data;
	void cache_warm() const {
		for (int i = 0, k = 0; i < size + 10; i++) k = data[k];
	}
public:
	Loop(int loop_size) {
		data = new int[loop_size]{ 0 };
		size = loop_size;
	}
	void forward_order() {
		for (int i = 1; i < size; i++) data[i - 1] = i;
		data[size - 1] = 0;
	}
	void revers_order() {
		for (int i = 1; i < size; i++) data[i] = i - 1;
		data[0] = size - 1;
	}
	void random_order() {
		for (int i = 0; i < size; i++) data[i] = i;
		for (int i = size - 1; i != 0; i--) std::swap(data[i], data[rand() % i]);
	}
	double tact_appeal() const {
		cache_warm();
		auto c_begin = __rdtsc();
		volatile int k = 0;
		for (size_t i = 0; i < size * round_count; i++) k = data[k];
		auto c_end = __rdtsc();
		return (double)(c_end - c_begin) / (round_count * size);
	}
	~Loop() {
		delete[] data;
	}
};

int main(int argc, char** argv) {
	if (argc == 1) return 1;
	int data_size = set_size(argv[1]);
	printf("Data size: %lu KB\n", data_size/KB);
	Loop loop(data_size/sizeof(int));
	loop.forward_order();
	double forward_at = loop.tact_appeal();
	loop.revers_order();
	double revers_at = loop.tact_appeal();
	loop.random_order();
	double random_at = loop.tact_appeal();
	printf("Forward: %3.2lf , revers: %3.2lf , random: %3.2lf" , forward_at, revers_at, random_at);
	return 0;
}
