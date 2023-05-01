#include <bits/stdc++.h>
#include "hash.h"
#include "tdigest.h"

#define inf 2147483647

template<typename ID_TYPE, typename DATA_TYPE>
class TdigestPerflow {
public:
	TdigestPerflow() {}
	TdigestPerflow(uint32_t memory, double compression) {
		latency = td_new(compression);
	}
	~TdigestPerflow() {
		delete latency;
	}
	void insert(ID_TYPE id, DATA_TYPE time) {
		td_add(latency, time, 1);
	}
	DATA_TYPE query(ID_TYPE id, double w) {

		return td_value_at(latency, w);

	}
private:
	td_histogram* latency;
	uint32_t size;
	uint32_t d;
};