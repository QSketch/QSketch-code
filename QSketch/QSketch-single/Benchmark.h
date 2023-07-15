#include <bits/stdc++.h>
#include "QuantileSketch.h"
#include "GroundTruth.h"
#include "hash.h"
#include "Mmap.h"


struct CAIDA_Tuple {
    uint64_t timestamp;
    uint64_t id;
};

#define MAXN 1000000
#define mod 10000019
#define web_len 383360

class CAIDA_Benchmark {
public:
	CAIDA_Benchmark() {}
	CAIDA_Benchmark(std::string path) {
		load_result = Load(path.c_str());
        dataset = (CAIDA_Tuple*)load_result.start;
        length = load_result.length / sizeof(CAIDA_Tuple);
	}
	~CAIDA_Benchmark() {}
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
	std::pair<double, double> Run(double w, uint32_t memory, std::vector<int> &len) {
        Init();
		uint32_t run_length = 20000000;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
		QuantileSketch<uint64_t, uint64_t> qs(memory, len);
        qs.set(query_quantile);
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

        /*for (int t = 1; t <= 3; t ++)
        {
            for (int i = 0; i < ins.size(); i++)
            {
                std::swap(ins[i], ins[(32894923ll * rand() + rand()) % (i + 1)]);
            }
        }*/

        for (int i = 0; i < ins.size(); i++)
        {
            gt.insert(1, ins[i]);
        }

        gt.build();

        tt = clock();

        for (int i = 0; i < ins.size(); i++)
        {
            qs.insert(1, ins[i]);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << totaltime << std::endl;

		double error_qs = 0;

        uint64_t predict_qs = qs.query(query_quantile);
        double predict_quantile_qs = gt.query(1, predict_qs);

        //std::cout << predict_quantile_qs << std::endl;

        error_qs += abs(predict_quantile_qs - query_quantile);


        return std::make_pair(throughput, error_qs);
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

struct zipf_Tuple {
    uint32_t id;
};


class zipf_Benchmark {
public:
	zipf_Benchmark() {}
	zipf_Benchmark(std::string path) {
		load_result = Load(path.c_str());
        dataset = (zipf_Tuple*)load_result.start;
        length = load_result.length / sizeof(zipf_Tuple);
	}
	~zipf_Benchmark() {}
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
    std::pair<double, double> Run(double w, uint32_t memory, std::vector<int> &len) {        Init();
		uint32_t run_length = 20000000;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
		QuantileSketch<uint64_t, uint64_t> qs(memory, len);
        qs.set(query_quantile);
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
                    ins.push_back(i - last_time[j]);
                    last_time[j] = i;
                    break;
                }
            }
            if (!flag)
            {
                cnt++;
                tid[cnt] = dataset[i].id; nxt[cnt] = head[g];
                head[g] = cnt;
                id_map[cnt] = 1;
                last_time[cnt] = i;
            }
		}

        for (int i = 0; i < ins.size(); i++)
        {
            std::swap(ins[i], ins[(32894923ll * rand() + rand()) % (i + 1)]);
        }

        std::cout << (int)ins.size() << std::endl;
        for (int i = 0; i < ins.size(); i++)
        {
            gt.insert(1, ins[i]);
        }

        //std::cout << ins.size() << std::endl;

        gt.build();

        tt = clock();

        for (int i = 0; i < ins.size(); i++)
        {
            qs.insert(1, ins[i]);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << totaltime << std::endl;

		double error_qs = 0;

        uint64_t predict_qs = qs.query(query_quantile);
        std::cout << predict_qs << std::endl;
        double predict_quantile_qs = gt.query(1, predict_qs);

        //std::cout << predict_quantile_qs << std::endl;

        error_qs += abs(predict_quantile_qs - query_quantile);


        return std::make_pair(throughput, error_qs);
	}

private:
	std::string filename;
    LoadResult load_result;
    zipf_Tuple *dataset;
    uint64_t length;
    int cnt;
	std::vector<int> head, nxt, id_map;
    std::vector<uint64_t> last_time;
	std::vector<uint64_t> tid; 
};

/*struct tnsm_Tuple {
    uint32_t unit_id;
    uint32_t fetch_time;
};

uint32_t stoi(char *start, char *end)
{
    uint32_t res = 0;
    for (char *c = start; c != end; c++) 
    {
        res = res * 10 + ((*c) - '0');
    }
    return res;
}

class tnsm_Benchmark {
public:
	tnsm_Benchmark() {}
	tnsm_Benchmark(std::string path) {
        freopen(path.c_str(), "r", stdin);
        length = web_len;
        char s[1111];
        std::vector<int> pos;
        int ln;
        dataset = new tnsm_Tuple [length];
        scanf("%[^\n]", s); getchar();
        for (int i = 0; i < length; i++)
        {
            scanf("%[^\n]", s); getchar();
            ln = strlen(s);
            pos.clear();
            for (int j = 0; j < ln; j++)
            {
                if (s[j] == '\"') pos.push_back(j);
            }
            dataset[i].unit_id = stoi(s + pos[0] + 1, s + pos[1]);
            dataset[i].fetch_time = stoi(s + pos[8] + 1, s + pos[9]);
        }
	}
	~tnsm_Benchmark()
    {
        delete dataset;
    }
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
	std::pair<double, double> Run(double w, double tower, uint32_t threshold, uint32_t memory, uint32_t d, uint32_t mm) {
        Init();
		uint32_t run_length = web_len;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
		QuantileSketch<uint64_t, uint64_t> qs(memory, threshold, tower, d, mm);
        qs.set(query_quantile);
		for (int i = 0; i < run_length; ++i) {
            int g = dataset[i].unit_id % mod;
            bool flag = 0;
            for (int j = head[g]; j; j = nxt[j])
            {
                if (tid[j] == dataset[i].unit_id)
                {
                    flag = 1;
                    id_map[j] ++;
                    gt.insert(dataset[i].unit_id, dataset[i].fetch_time);

                    tt = clock();
                    qs.insert(dataset[i].unit_id, dataset[i].fetch_time);
                    tottime += (clock() - tt);
                    
                    break;
                }
            }
            if (!flag)
            {
                cnt++;
                tid[cnt] = dataset[i].unit_id; nxt[cnt] = head[g];
                head[g] = cnt;
                id_map[cnt] = 1;

                gt.insert(dataset[i].unit_id, dataset[i].fetch_time);

                tt = clock();
                qs.insert(dataset[i].unit_id, dataset[i].fetch_time);
                tottime += (clock() - tt);
            }
		}


        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << totaltime << std::endl;

		double error_qs = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) {
            if (id_map[i] < 1000)
                continue;
            num++;

            uint64_t predict_qs = qs.query(tid[i], query_quantile);
            double predict_quantile_qs = gt.query(tid[i], predict_qs);
            error_qs += abs(predict_quantile_qs - query_quantile);
            //std::cout << predict_quantile_qs << "\n";
        }

        return std::make_pair(throughput, error_qs / num);
	}

private:
	std::string filename;
    LoadResult load_result;
    tnsm_Tuple *dataset;
    uint64_t length;
    int cnt;
	std::vector<int> head, nxt, id_map;
    std::vector<uint64_t> last_time;
	std::vector<uint64_t> tid; 
};*/

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
    std::pair<double, double> Run(double w, uint32_t memory, std::vector<int> &len) {
		uint32_t run_length = 20000000;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
		QuantileSketch<uint64_t, uint64_t> qs(memory, len);
        qs.set(query_quantile);
        std::vector<uint64_t> ins; ins.clear();
		for (int i = 0; i < run_length; ++i) {
            ins.push_back(dataset[i].stamp);
		}

        for (int i = 0; i < ins.size(); i++)
        {
            gt.insert(1, ins[i]);
        }

        gt.build();

        tt = clock();

        for (int i = 0; i < ins.size(); i++)
        {
            qs.insert(1, ins[i]);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << totaltime << std::endl;

		double error_qs = 0;

        uint64_t predict_qs = qs.query(query_quantile);
        double predict_quantile_qs = gt.query(1, predict_qs);

        //std::cout << predict_quantile_qs << std::endl;

        error_qs += abs(predict_quantile_qs - query_quantile);


        return std::make_pair(throughput, error_qs);
	}

private:
    synthetic_Tuple *dataset;
    uint64_t length;
};