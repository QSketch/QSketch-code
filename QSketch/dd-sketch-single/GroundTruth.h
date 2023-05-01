#ifndef _GROUND_TRUTH_H_
#define _GROUND_TRUTH_H_

#include <bits/stdc++.h>
#include "hash.h"

#define MAXN 1000000
#define mod 10000019

template<typename ID_TYPE, typename DATA_TYPE>
class GroundTruth {
public:
	GroundTruth()
	{
		nxt.resize(MAXN + 5);
		tid.resize(MAXN + 5);
		head.resize(mod + 5);
		id_latency.resize(MAXN + 5);
		cnt = 0;
		for (int i = 1; i <= MAXN; i++) 
		{
			id_latency[i].clear();
			nxt[i] = 0;
		}
		for (int i = 0; i < mod; i++) head[i] = 0;
	}
	~GroundTruth() {}
	void insert(ID_TYPE id, DATA_TYPE time) {
		int g = id % mod;	
		for (int i = head[g]; i; i = nxt[i])
		{
			if (tid[i] == id)
			{
				id_latency[i].push_back(time);
				return;
			} 
		}
		++cnt;
		tid[cnt] = id; nxt[cnt] = head[g]; head[g] = cnt;
		id_latency[cnt].push_back(time);
	}
	double query(ID_TYPE id, DATA_TYPE time) {
		int g = id % mod;
		for (int i = head[g]; i; i = nxt[i])
		{
			if (tid[i] == id)
			{
				std::sort(id_latency[i].begin(), id_latency[i].end());
				int index = lower_bound(id_latency[i].begin(), id_latency[i].end(), time) - id_latency[i].begin();
				double quantile = 1.0 * index / id_latency[i].size();
				return quantile;
			}
		}
		assert(0);
		return 0;
	}

private:
	int cnt;
	std::vector<std::vector<DATA_TYPE> > id_latency;
	std::vector<int> head;
	std::vector<int> nxt;
	std::vector<ID_TYPE> tid; 
};


#endif