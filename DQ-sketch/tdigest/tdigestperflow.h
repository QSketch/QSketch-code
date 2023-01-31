#include <bits/stdc++.h>
#include "hash.h"
#include "tdigest.h"

#define inf 2147483647

template<typename ID_TYPE, typename DATA_TYPE>
class TdigestPerflow {
public:
	TdigestPerflow() {}
	TdigestPerflow(uint32_t memory, double compression) {
		d = 3;
		size = memory * 1024 / d / td_required_buf_size(compression);
		assert(size >= 5);
		//std::cout << size << "\n";
		frequency = new uint32_t* [d];
		latency = new td_histogram** [d];
		for (int i = 0; i < d; ++i) {
			frequency[i] = new uint32_t [size];
			latency[i] = new td_histogram* [size];
			for (int j = 0; j < size; ++j) {
				latency[i][j] = td_new(compression);
			}
			memset(frequency[i], 0, size * sizeof(uint32_t));
		} 
	}
	~TdigestPerflow() {
		for (int i = 0; i < d; ++i) {
			for (int j = 0; j < size; ++j) {
				delete latency[i][j];
			}
			delete[] frequency[i];
			delete[] latency[i];
		}
		delete[] frequency;
		delete[] latency;
	}
	void insert(ID_TYPE id, DATA_TYPE time) {
		uint32_t minf = inf;
		bool flag;
		for (int i = 0; i < d; i++)
		{
			uint32_t index = hash(id, i + 10) % size;
			if (minf > frequency[i][index]) minf = frequency[i][index];
		}
		for (int i = 0; i < d; ++i) {
			uint32_t index = hash(id, i + 10) % size;
			if (frequency[i][index] != minf) continue;
			frequency[i][index]++;
			td_add(latency[i][index], time, 1);
		}
	}
	DATA_TYPE query(ID_TYPE id, double w) {

		uint32_t min_frequency = inf, min_i;
		for (int i = 0; i < d; ++i) {
			uint32_t index = hash(id, i + 10) % size;
			if (frequency[i][index] < min_frequency) {
				min_frequency = frequency[i][index];
				min_i = i;
			}
		}
		uint32_t index_i = hash(id, min_i + 10) % size;
		return td_value_at(latency[min_i][index_i], w);

	}
private:
	uint32_t** frequency;
	td_histogram*** latency;
	uint32_t size;
	uint32_t d;
};