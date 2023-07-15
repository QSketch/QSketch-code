
#include <bits/stdc++.h>
#include "hash.h"
#include "GroundTruth.h"
#define inf 2147483647

template<typename ID_TYPE, typename DATA_TYPE>
class LatencyRecord {
public:
	LatencyRecord() {}
	LatencyRecord(int memory, std::vector<int> &ln) {
		uint32_t size = memory * 1024 / sizeof(DATA_TYPE);
		for (int i = 0; i < (int)ln.size(); i++)
		{
			size -= ln[i];
		}
		assert(size > 10);
		layers = (int)ln.size() + 1;
		len = new uint32_t [layers];
		for (int i = 0; i < layers - 1; i++) len[i] = ln[i];
		len[layers - 1] = size;
		buffer = new DATA_TYPE* [layers];
		for (int i = 0; i < layers; i++)
		{
			buffer[i] = new DATA_TYPE [len[i]];
			memset(buffer[i], 0, sizeof(DATA_TYPE) * len[i]);
		}
		buffer_size = new uint32_t [layers];
		memset(buffer_size, 0, sizeof(uint32_t) * layers);
		num = 0;
		//for (int i = 0; i < layers; i++) std::cout << len[i] << " ";
		//std::cout << std::endl;

		ss = 0.00;
	}
	~LatencyRecord() {
		for (int i = 0; i < layers; i++) delete[] buffer[i];
		delete[] buffer;
		delete[] len;
		delete[] buffer_size;
	}
	void set(double w) {
		quantile = w;
		if (w > 0.5) p = (2.00 * w - 1.00) / (2.00 * w);
		else p = (1.00 - 2.00 * w) / (2.00 * (1.00 - w));
	}
	void insert(DATA_TYPE time){
		while(rand() < RAND_MAX * p) num++;
		buffer[0][buffer_size[0] ++] = time;
		std::vector<DATA_TYPE> v, vv; v.clear();
		if (buffer_size[0] + num >= len[0])
		{
			bool flag = 0;
			//std::cout << buffer_size[0] << " " << num << " " << 0.6 * buffer_size[0] << " " << num - 0.6 * buffer_size[0] << std::endl;
			std::sort(buffer[0], buffer[0] + buffer_size[0]);
			for (int i = -10; i <= 10; i++)
			{
				if (quantile < 0.5) 
				{
					if ((int)(num + buffer_size[0]) / 2 + i >= num) time = buffer[0][(num + buffer_size[0]) / 2 + i - num];
					else time = 0;
				}
				else 
				{
					if ((int)(num + buffer_size[0]) / 2 + i < buffer_size[0]) time = buffer[0][(num + buffer_size[0]) / 2 + i];
					else time = -1ull;
				}
				v.push_back(time);
			}
			num = 0;
			buffer_size[0] = 0;
			for (int i = 1; i < layers - 1; i++)
			{
				vv.clear();
				for (int r = 0; r < v.size(); r++)
				{
					time = v[r];
					if (time > buffer[i][0])
					{
						bool flg = 0;
						for (int j = 1; j < len[i]; j++)
						{
							if (buffer[i][j] >= time)
							{
								flg = 1;
								buffer[i][j - 1] = time;
								break;
							}
							buffer[i][j - 1] = buffer[i][j];
						}
						if (!flg) buffer[i][len[i] - 1] = time;
					}
					buffer_size[i] ++;
					if (buffer_size[i] == 2 * len[i] - 1)
					{
						time = buffer[i][0];
						memset(buffer[i], 0, sizeof(DATA_TYPE) * len[i]);
						buffer_size[i] = 0;
						vv.push_back(time);
					}
				}
				v = vv;
				if (v.empty())
				{
					flag = 1;
					break;
				}
			}
			if (!flag)
			{
				//ss += (gt -> query(1, time) - quantile);
				//std::cout << gt -> query(1, time) << std::endl;
				time = v[0];
				if (buffer[layers - 1][0] == 0) 
				{
					bool flg = 0;
					for (int j = 1; j < len[layers - 1]; j++)
					{
						if (buffer[layers - 1][j] > time) 
						{
							buffer[layers - 1][j - 1] = time; 
							flg = 1;
							break;
						}
						buffer[layers - 1][j - 1] = buffer[layers - 1][j];
					}
					if (!flg) buffer[layers - 1][len[layers - 1] - 1] = time;
				}
				else {
					if (rand() < 0.5 * RAND_MAX)
					{
						// replace minimum value
						if (time < buffer[layers - 1][0]) return;
						flag = 0;
						for (int j = 1; j < len[layers - 1]; j++)
						{
							if (buffer[layers - 1][j] > time) 
							{
								buffer[layers - 1][j - 1] = time; 
								flag = 1;
								break;
							}
							buffer[layers - 1][j - 1] = buffer[layers - 1][j];
						}
						if (!flag) buffer[layers - 1][len[layers - 1] - 1] = time;
					}
					else 
					{
						// replace maximum value
						if (time > buffer[layers - 1][len[layers - 1] - 1]) return;
						flag = 0;
						for (int j = len[layers - 1] - 2; j >= 0; j--)
						{
							if (buffer[layers - 1][j] < time) 
							{
								buffer[layers - 1][j + 1] = time; 
								flag = 1;
								break;
							}
							buffer[layers - 1][j + 1] = buffer[layers - 1][j];
						}
						if (!flag) buffer[layers - 1][0] = time;
					}
				}
				/*for (int i = 0; i < len[layers - 1]; i++)
				{
					printf("%.4lf ", gt -> query(1, buffer[layers - 1][i]));
				}
				printf("\n");*/
			}
		}
		/*for (int i = 0; i < buffer_size[0]; i++) printf("%.4lf ", gt -> query(1, buffer[0][i]));
		printf("\n");
		for (int i = 1; i < layers; i++)
		{
			for (int j = 0; j < len[i]; j++)
			{
				printf("%.4lf ", gt -> query(1, buffer[i][j]));
			}
			printf("\n");
		}
		printf("\n");*/
		/*if (rand() % 100000 < 10)
		{
			for (int i = 0; i < len[layers - 1]; i++)
			{
				printf("%.4lf ", gt -> query(1, buffer[layers - 1][i]));
			}
			printf("\n");
		}*/
	}
	DATA_TYPE query(double w)
	{
		//for (int i = 0; i < len[layers - 1]; i++) std::cout << buffer[layers - 1][i] << " ";
		//std::cout << std::endl;
		//std::cout << ss << std::endl;
		return buffer[layers - 1][int(0.50 * (len[layers - 1] - 1))];
	}
private:
	DATA_TYPE** buffer;
	uint32_t* buffer_size;
	uint32_t* len;
	uint32_t num, layers;
	double quantile;
	double p;
	double ss;
};

template<typename ID_TYPE, typename DATA_TYPE>
class QuantileSketch {
public:
	QuantileSketch() {}
	QuantileSketch(uint32_t memory, std::vector<int> &len) {
		latencyrecord = new LatencyRecord<ID_TYPE, DATA_TYPE>(memory, len);
	}
	~QuantileSketch() {
		delete latencyrecord;
	}
	void set(double w) {
		latencyrecord->set(w);
	}
	void insert(ID_TYPE id, DATA_TYPE time) {
		latencyrecord->insert(time);
	}
	DATA_TYPE query(double w) {
		return latencyrecord->query(w);
	}

private:
	LatencyRecord<ID_TYPE, DATA_TYPE>* latencyrecord;
};

