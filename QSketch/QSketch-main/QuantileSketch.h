
#include <bits/stdc++.h>
#include "hash.h"
#include "GroundTruth.h"
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
	LatencyRecord(uint32_t memory, int mm, int dd, double l, std::vector<int> &ln) {
		d = dd;
		lambda = l;
		layers = (int)ln.size() + 1;
		len = new uint32_t [layers];
		for (int i = 0; i < layers - 1; i++) len[i] = ln[i];
		len[layers - 1] = mm;
		uint32_t sum = 0;
		for (int i = 0; i < layers; i++) sum += len[i];
		sum *= sizeof(DATA_TYPE);
		sum += ((4 + layers) * sizeof(uint32_t) + sizeof(ID_TYPE));
		size = memory * 1024 / d / sum;
		std::cerr << size << std::endl;

		vote = new uint32_t* [size];
		frequency = new uint32_t* [size];
		num = new uint32_t* [size];
		ids = new ID_TYPE* [size];
		buffer_size = new uint32_t** [size];
		buffer = new DATA_TYPE*** [size];
		flags = new uint32_t* [size];

		for (int i = 0; i < size; ++i)
		{
			vote[i] = new uint32_t [d];
			memset(vote[i], 0, sizeof(uint32_t) * d);
			frequency[i] = new uint32_t [d];
			memset(frequency[i], 0, sizeof(uint32_t) * d);
			num[i] = new uint32_t [d];
			memset(num[i], 0, sizeof(uint32_t) * d);
			ids[i] = new ID_TYPE [d];
			memset(ids[i], 0, sizeof(ID_TYPE) * d);
			flags[i] = new uint32_t [d];
			memset(flags[i], 0, sizeof(uint32_t) * d);
			buffer_size[i] = new uint32_t* [d];
			buffer[i] = new DATA_TYPE** [d];
			for (int j = 0; j < d; j++)
			{
				buffer_size[i][j] = new uint32_t [layers];
				memset(buffer_size[i][j], 0, sizeof(uint32_t) * layers);
				buffer[i][j] = new DATA_TYPE* [layers];
				for (int k = 0; k < layers; k++)
				{
					buffer[i][j][k] = new DATA_TYPE [len[k]];
					memset(buffer[i][j][k], 0, sizeof(DATA_TYPE) * len[k]);
				}
			}
		} 
	}
	~LatencyRecord() {
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < d; j++)
			{
				delete[] buffer_size[i][j];
				for (int k = 0; k < layers; k++)
				{
					delete[] buffer[i][j][k];
				}
				delete[] buffer[i][j];
			}
			delete[] buffer_size[i];
			delete[] buffer[i];
			delete[] vote[i];
			delete[] frequency[i];
			delete[] ids[i];
			delete[] num[i];
			delete[] flags[i];
		}
		delete[] vote;
		delete[] frequency;
		delete[] ids;
		delete[] buffer_size;
		delete[] buffer;
		delete[] num;
		delete[] len;
		delete[] flags;
	}
	void set(double w) {
		quantile = w;
		if (w > 0.5) p = (2.00 * w - 1.00) / (2.00 * w);
		else p = (1.00 - 2.00 * w) / (2.00 * (1.00 - w));
	}
	void insert(ID_TYPE id, DATA_TYPE time) {
		uint32_t index = hash(id, 10) % size;
		std::vector<DATA_TYPE> v, vv; v.clear(); vv.clear();
		for (int i = 0; i < d; i++)
		{
			if (ids[index][i] == id)
			{
				frequency[index][i]++;
				buffer[index][i][0][buffer_size[index][i][0] ++] = time;
				while(rand() < RAND_MAX * p) num[index][i] ++;
				if (buffer_size[index][i][0] + num[index][i] >= len[0])
				{
					bool flag = 0; v.clear();
					std::sort(buffer[index][i][0], buffer[index][i][0] + buffer_size[index][i][0]);
					if (quantile < 0.5) 
					{
						if ((int)(num[index][i] + buffer_size[index][i][0]) / 2 >= num[index][i]) 
							time = buffer[index][i][0][(num[index][i] + buffer_size[index][i][0]) / 2 - num[index][i]];
						else time = 0;
					}
					else 
					{
						if ((int)(num[index][i] + buffer_size[index][i][0]) / 2 < buffer_size[index][i][0]) 
							time = buffer[index][i][0][(num[index][i] + buffer_size[index][i][0]) / 2];
						else time = -1ull;
					}
					v.push_back(time);
					num[index][i] = 0;
					buffer_size[index][i][0] = 0;
				
					for (int k = 1; k < layers - 1; k++)
					{
						vv.clear();
						for (int r = 0; r < v.size(); r++)
						{
							time = v[r];
							if (time > buffer[index][i][k][0])
							{
								bool flg = 0;
								for (int j = 1; j < len[k]; j++)
								{
									if (buffer[index][i][k][j] >= time)
									{
										flg = 1;
										buffer[index][i][k][j - 1] = time;
										break;
									}
									buffer[index][i][k][j - 1] = buffer[index][i][k][j];
								}
								if (!flg) buffer[index][i][k][len[k] - 1] = time;
							}
							buffer_size[index][i][k] ++;
							if (buffer_size[index][i][k] == 2 * len[k] - 1)
							{
								time = buffer[index][i][k][0];
								memset(buffer[index][i][k], 0, sizeof(DATA_TYPE) * len[k]);
								buffer_size[index][i][k] = 0;
								vv.push_back(time);
							}
						}
						if (vv.empty()) 
						{
							flag = 1;
							break;
						}
						v = vv;
					}
					if (!flag)
					{
						if (buffer[index][i][layers - 1][0] == 0) 
						{
							bool flg = 0;
							for (int j = 1; j < len[layers - 1]; j++)
							{
								if (buffer[index][i][layers - 1][j] > time) 
								{
									buffer[index][i][layers - 1][j - 1] = time; 
									flg = 1;
									break;
								}
								buffer[index][i][layers - 1][j - 1] = buffer[index][i][layers - 1][j];
							}
							if (!flg) buffer[index][i][layers - 1][len[layers - 1] - 1] = time;
						}
						else {
							if (!flags[index][i])
							{
								// replace minimum value
								if (time < buffer[index][i][layers - 1][0]) return;
								bool flg = 0;
								for (int j = 1; j < len[layers - 1]; j++)
								{
									if (buffer[index][i][layers - 1][j] > time) 
									{
										buffer[index][i][layers - 1][j - 1] = time; 
										flg = 1;
										break;
									}
									buffer[index][i][layers - 1][j - 1] = buffer[index][i][layers - 1][j];
								}
								if (!flg) buffer[index][i][layers - 1][len[layers - 1] - 1] = time;
							}
							else 
							{
								// replace maximum value
								if (time > buffer[index][i][layers - 1][len[layers - 1] - 1]) return;
								bool flg = 0;
								for (int j = len[layers - 1] - 2; j >= 0; j--)
								{
									if (buffer[index][i][layers - 1][j] < time) 
									{
										buffer[index][i][layers - 1][j + 1] = time; 
										flg = 1;
										break;
									}
									buffer[index][i][layers - 1][j + 1] = buffer[index][i][layers - 1][j];
								}
								if (!flg) buffer[index][i][layers - 1][0] = time;
							}
							flags[index][i] ^= 1;
						}
					}
				}
				/*if (id == 1722875553ull)
				{
					for (int j = 0; j < layers; j++)
					{
						for (int k = 0; k < len[j]; k++)
						{
							std::cout << gt -> query(id, buffer[index][i][j][k]) << " ";
						}
						std::cout << std::endl;
					}
					std::cout << std::endl;
				}*/
				return;
			}
		}
		for (int i = 0; i < d; i++)
		{
			if (!ids[index][i])
			{
				frequency[index][i]++;
				buffer[index][i][0][buffer_size[index][i][0] ++] = time;
				ids[index][i] = id;
				while(rand() < RAND_MAX * p) num[index][i] ++;
				//assert(buffer_size[index][i][0] + num[index][i] < len[0]);
				return;
			}
		}
		int mn_i, mn = 1e9;
		for (int i = 0; i < d; i++)
		{
			if (frequency[index][i] < mn)
			{
				mn = frequency[index][i];
				mn_i = i;
			}
		}
		vote[index][mn_i] ++;
		if (vote[index][mn_i] >= lambda * frequency[index][mn_i])
		{
			vote[index][mn_i] = frequency[index][mn_i] = ids[index][mn_i] = 0;
			num[index][mn_i] = 0;
			flags[index][mn_i] = 0;
			for (int i = 0; i < layers; i++) 
			{
				buffer_size[index][mn_i][i] = 0;
				memset(buffer[index][mn_i][i], 0, sizeof(DATA_TYPE) * len[i]);
			}
		}
	}
	DATA_TYPE query(ID_TYPE id, double w)
	{
		uint32_t index = hash(id, 10) % size;
		for (int i = 0; i < d; i++)
		{
			if (ids[index][i] == id) return buffer[index][i][layers - 1][int(0.5 * (len[layers - 1] - 1) + 0.1)];
		}
		return buffer[index][0][layers - 1][int(0.5 * (len[layers - 1] - 1) + 0.1)];
	}
private:
	uint32_t** vote;
	uint32_t** frequency;
	uint32_t** num;
	uint32_t*** buffer_size;
	uint32_t** flags;
	uint32_t* len;
	ID_TYPE** ids;
	DATA_TYPE**** buffer;
	uint32_t size;
	uint32_t d;
	uint32_t layers;
	double quantile;
	double p;
	double lambda;
};

template<typename ID_TYPE, typename DATA_TYPE>
class QuantileSketch {
public:
	QuantileSketch() {}
	QuantileSketch(uint32_t memory, uint32_t thre, double tower, uint32_t mm, uint32_t dd, double lambda, std::vector<int> &len) {
		towersketch = new TowerSketch<ID_TYPE>(tower * memory);
		latencyrecord = new LatencyRecord<ID_TYPE, DATA_TYPE>((1.00 - tower) * memory, mm, dd, lambda, len);
		threshold = thre; 
	}
	~QuantileSketch() {
		delete towersketch;
		delete latencyrecord;
	}
	void set(double w) {
		latencyrecord->set(w);
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

