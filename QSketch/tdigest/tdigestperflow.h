#include <bits/stdc++.h>
#include "hash.h"
#include "tdigest.h"

#define inf 2147483647

template<typename ID_TYPE, typename DATA_TYPE>
class TdigestPerflow {
public:
	TdigestPerflow() {}
	TdigestPerflow(uint32_t memory, double compression) {
		size = memory * 1024 / td_required_buf_size(compression);
		assert(size >= 5);
		//std::cout << size << "\n";
		latency = new td_histogram* [size];
		for (int i = 0; i < size; i++) 
		{
			latency[i] = td_new(compression);
		}
	}
	~TdigestPerflow() {
		for (int i = 0; i < size; i++)
		{
			delete latency[i];
		}
		delete[] latency;
	}
	void insert(ID_TYPE id, DATA_TYPE time) {
		uint32_t index = hash(id, 1) % size;

		td_add(latency[index], time, 1);
	}
	DATA_TYPE query(ID_TYPE id, double w) {

		uint32_t index = hash(id, 1) % size;
		return td_value_at(latency[index], w);
	}
private:
	td_histogram** latency;
	uint32_t size;
};