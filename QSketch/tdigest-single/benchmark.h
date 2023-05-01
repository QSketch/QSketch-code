#include <bits/stdc++.h>
#include "tdigestperflow.h"
#include "GroundTruth.h"
#include "hash.h"
#include "Mmap.h"


struct CAIDA_Tuple {
    uint64_t timestamp;
    uint64_t id;
};

#define MAXN 1000000
#define mod 10000019

class Benchmark {
public:
	Benchmark() {}
	Benchmark(std::string path) {
		load_result = Load(path.c_str());
        dataset = (CAIDA_Tuple*)load_result.start;
        length = load_result.length / sizeof(CAIDA_Tuple);
	}
	~Benchmark() {}
    void Init()
    {
        nxt.resize(MAXN + 5);
		tid.resize(MAXN + 5);
		head.resize(mod + 5);
        id_map.resize(MAXN + 5);
        last_time.resize(MAXN + 5);
		cnt = 0;
		for (int i = 1; i <= MAXN; i++) 
		{
			nxt[i] = 0;
            id_map[i] = 0;
		}
		for (int i = 0; i < mod; i++) head[i] = 0;
    }
	std::pair<double, double> Run(double w, uint32_t memory, double compression) {
        Init();
		uint32_t run_length = 20000000;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
		TdigestPerflow<uint64_t, uint64_t> td(memory, compression);

        std::vector<uint64_t> ins; ins.clear();
		for (int i = 0; i < run_length; ++i) {
            int g = dataset[i].id % mod;
            bool flag = 0;
            for (int j = head[g]; j; j = nxt[j])
            {
                if (tid[j] == dataset[i].id)
                {
                    flag = 1;
                    id_map[j] ++;
                    if (dataset[i].timestamp > last_time[j])
                    {
                        gt.insert(1, dataset[i].timestamp - last_time[j]);

                        ins.push_back(dataset[i].timestamp - last_time[j]);
                    }
                    last_time[j] = dataset[i].timestamp;
                    break;
                }
            }
            if (!flag)
            {
                cnt++;
                tid[cnt] = dataset[i].id; nxt[cnt] = head[g];
                head[g] = cnt;
                id_map[cnt] = 1;
                last_time[cnt] = dataset[i].timestamp;
            }
		}

        tt = clock();
        
        for (int i = 0; i < ins.size(); i++)
        {
            td.insert(1, ins[i]);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << "insertion ends\n";

		double error_td = 0;

        uint64_t predict_td = td.query(1, query_quantile);
        double predict_quantile_td = gt.query(1, predict_td);
        error_td += abs(predict_quantile_td - query_quantile);

        return std::make_pair(throughput, error_td);
	}

private:
	std::string filename;
    LoadResult load_result;
    CAIDA_Tuple *dataset;
    uint64_t length;
    int cnt;
	std::vector<int> head, nxt, id_map;
    std::vector<uint64_t> last_time;
	std::vector<uint64_t> tid; 
};

struct synthetic_Tuple {
    uint64_t stamp;
};


class synthetic_Benchmark {
public:
	synthetic_Benchmark() {}
	synthetic_Benchmark(std::string path) {
		freopen(path.c_str(), "r", stdin);
        length = 20000000;
        dataset = new synthetic_Tuple [length];
        for (int i = 0; i < length; i++)
        {
            scanf("%llu", &dataset[i].stamp);
        }
	}
	~synthetic_Benchmark() {}
    std::pair<double, double> Run(double w, uint32_t memory, double compression) {
        uint32_t run_length = 20000000;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
		TdigestPerflow<uint64_t, uint64_t> td(memory, compression);

        std::vector<uint64_t> ins; ins.clear();
		for (int i = 0; i < run_length; ++i) {
            gt.insert(1, dataset[i].stamp);
            ins.push_back(dataset[i].stamp);
		}

        tt = clock();
        
        for (int i = 0; i < ins.size(); i++)
        {
            td.insert(1, ins[i]);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << "insertion ends\n";

		double error_td = 0;

        uint64_t predict_td = td.query(1, query_quantile);
        double predict_quantile_td = gt.query(1, predict_td);
        error_td += abs(predict_quantile_td - query_quantile);

        return std::make_pair(throughput, error_td);
	}

private:
    synthetic_Tuple *dataset;
    uint64_t length;
};