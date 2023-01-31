#ifndef _QUANTILE_SKETCH_H_
#define _QUANTILE_SKETCH_H_

#include <bits/stdc++.h>
#include "hash.h"

#define inf 2147483647

const int bits[] = {2, 4, 8};
const int max_value[] = {2, 14, 254};

template<typename ID_TYPE>
class TowerSketch {
public:
	TowerSketch() {}
	TowerSketch(uint32_t memory) {
		c_num = memory * 1024 / 3 / sizeof(uint32_t);
		counter = new uint32_t* [3];
		for (int i = 0; i < 3; ++i) {
			size[i] = c_num * 32 / bits[i];
			counter[i] = new uint32_t [c_num];
			memset(counter[i], 0, c_num * sizeof(uint32_t));
		}
	}
	~TowerSketch() {
		for (int i = 0; i < 3; ++i) {
			delete[] counter[i];
		}
		delete[] counter;
	}
	uint32_t get(int i, int j, int start, int end) {
		return (counter[i][j] & (((1 << (end - start)) - 1) << start)) >> start;
	}
	void set(int i, int j, int start, int end, uint32_t value) {
		counter[i][j] &= (~(((1 << (end - start)) - 1) << start));
		counter[i][j] |= (value << start);
	}
	void insert(ID_TYPE id) {
		for (int i = 0; i < 3; ++i) {
			uint32_t tag = hash(id, i + 5) % size[i];
			uint32_t index = tag * bits[i] / 32, res = tag * bits[i] - index * 32;
			uint32_t value = get(i, index, res, res + bits[i]);
			if (value == max_value[i]) {
				continue;
			}
			set(i, index, res, res + bits[i], value + 1);
			assert(get(i, index, res, res + bits[i]) == value + 1);
		}
	}
	uint32_t query(ID_TYPE id) {
		int result = inf;
		for (int i = 0; i < 3; ++i) {
			uint32_t tag = hash(id, i + 5) % size[i];
			uint32_t index = tag * bits[i] / 32, res = tag * bits[i] - index * 32;
			uint32_t value = get(i, index, res, res + bits[i]);
			if (value < max_value[i]) 
				result = MIN(result, value);
		}
		return result;
	}
private:
	int size[3];
	int c_num;
	uint32_t** counter;
};

template<typename ID_TYPE, typename DATA_TYPE>
class LatencyRecord {
public:
	LatencyRecord() {}
	LatencyRecord(uint32_t memory, int dd, int mm) {
		d = dd;
		M = mm;
		size = memory * 1024 / d / (M * sizeof(DATA_TYPE) + sizeof(uint32_t));
		frequency = new uint32_t* [d];
		latency = new DATA_TYPE** [d];
		for (int i = 0; i < d; ++i) {
			frequency[i] = new uint32_t [size];
			latency[i] = new DATA_TYPE* [size];
			for (int j = 0; j < size; ++j) {
				latency[i][j] = new DATA_TYPE [M];
				memset(latency[i][j], 0, M * sizeof(DATA_TYPE));
			}
			memset(frequency[i], 0, size * sizeof(uint32_t));
		} 
	}
	~LatencyRecord() {
		for (int i = 0; i < d; ++i) {
			for (int j = 0; j < size; ++j) {
				delete[] latency[i][j];
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
		int pos;
		for (int i = 0; i < d; ++i) {
			uint32_t index = hash(id, i + 10) % size;
			if (frequency[i][index] != minf) continue;
			frequency[i][index]++;
			if (frequency[i][index] <= M) {
				pos = 0;
				for (int j = frequency[i][index] - 2; j >= 0; j--)
				{
					if (latency[i][index][j] < time)
					{
						pos = j + 1;
						break;
					}
					latency[i][index][j + 1] = latency[i][index][j];
				}
				latency[i][index][pos] = time;
			}
			else {
				uint32_t x = ((238947893ll * rand()) + rand()) % frequency[i][index];
				if (x < M)
				{
					latency[i][index][x] = time;
					while(x < M - 1 && latency[i][index][x] > latency[i][index][x + 1])
					{
						std::swap(latency[i][index][x], latency[i][index][x + 1]);
						x++;
					}
					while(x > 0 && latency[i][index][x] < latency[i][index][x - 1])
					{
						std::swap(latency[i][index][x], latency[i][index][x - 1]);
						x--;
					}
				}
			}
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
		return latency[min_i][index_i][int(w * (M - 1))];

	}
private:
	uint32_t** frequency;
	DATA_TYPE*** latency;
	uint32_t size;
	uint32_t M;
	uint32_t d;
};

template<typename ID_TYPE, typename DATA_TYPE>
class QuantileSketch {
public:
	QuantileSketch() {}
	QuantileSketch(uint32_t memory, uint32_t thre, double tower, uint32_t d, uint32_t mm) {
		towersketch = new TowerSketch<ID_TYPE>(tower * memory);
		latencyrecord = new LatencyRecord<ID_TYPE, DATA_TYPE>((1.00 - tower) * memory, d, mm);
		threshold = thre; 
	}
	~QuantileSketch() {
		delete towersketch;
		delete latencyrecord;
	}
	void insert(ID_TYPE id, DATA_TYPE time) {
		if (towersketch->query(id) < threshold) {
			towersketch->insert(id);
			return;
		}
		latencyrecord->insert(id, time);
	}
	DATA_TYPE query(ID_TYPE id, double w) {
		return latencyrecord->query(id, w);
	}

private:
	TowerSketch<ID_TYPE>* towersketch;
	LatencyRecord<ID_TYPE, DATA_TYPE>* latencyrecord;
	uint32_t threshold;
};



#endif