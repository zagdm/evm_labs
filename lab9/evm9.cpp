#include <iostream>
#include <cmath>
#include <x86intrin.h>

static const size_t KB = 1024;
static const size_t L1_size = 32 * KB;
static const size_t L3_size = 4 * KB * KB;

class Loop {
private:
	static const size_t offset = L3_size / sizeof(size_t);
	static const size_t round_count = 10000;
	size_t fragments;
	size_t* data;
	void cache_warm() const {
		for (size_t k = data[0]; k != 0; k = data[k]);
	}
	size_t fragment_size() const {
		return L1_size / (fragments * sizeof(size_t));
	}
	size_t element_count() const {
		return fragments * fragment_size();
	}
public:
	Loop(size_t fragment_count) {
		fragments = fragment_count;
		data = new size_t[fragments * offset]{ 0 };
	}
	void set_order() {
		for (size_t i = 0; i < fragment_size(); i++) {
			for (size_t f = 0; f < fragments - 1; f++) data[i + f * offset] = i + (f + 1) * offset;
			data[i + (fragments - 1) * offset] = i + 1;
		}
		data[fragment_size() - 1 + (fragments - 1) * offset] = 0;
	}
	double tact_appeal() const {
		cache_warm();
		auto c_begin = __rdtsc();
		volatile size_t k = 0;
		for (size_t i = 0; i < element_count() * round_count; i++) k = data[k];
		auto c_end = __rdtsc();
		return (double)(c_end - c_begin) / (element_count() * round_count);
	}
	~Loop() {
		delete[] data;
	}
};

int main(int argc, char** argv) {
	if (argc == 1) return 1;
	size_t fragment_count = (size_t)atoi(argv[1]);
	Loop loop(fragment_count);
	printf("Fragment count: %lu\n", fragment_count);
	loop.set_order();
	printf("Tact count: %3.2lf", loop.tact_appeal());
	return 0;
}
