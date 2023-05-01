#include <bits/stdc++.h>
#include "DDSketch.h"
#include "GroundTruth.h"
#include "hash.h"
#include "Mmap.h"


struct CAIDA_Tuple {
    uint64_t timestamp;
    uint64_t id;
};

#define MAXN 20000000
#define mod 10000019
#define n_slice 688
#define MOD 13908490328490391

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
	std::pair<double, double> Run(double w, double tower, uint32_t threshold, uint32_t memory, uint32_t d, uint32_t mm) {
        Init();
		uint32_t run_length = 20000000;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
        DDSketch<uint64_t, uint64_t> dd(memory);
        std::vector<std::pair<uint64_t, uint64_t> > ins; ins.clear();

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
                        gt.insert(dataset[i].id, dataset[i].timestamp - last_time[j]);

                        ins.push_back(std::make_pair(dataset[i].id, dataset[i].timestamp - last_time[j]));

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
            dd.insert(ins[i].first, ins[i].second);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << totaltime << std::endl;

		double error_dd = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) {
            if (id_map[i] < 5000)
                continue;
            num++;
            
            uint64_t predict_dd = dd.query(tid[i], query_quantile);
            double predict_quantile_dd = gt.query(tid[i], predict_dd);
            //printf("%.6lf\n", predict_quantile_dd);
            error_dd += abs(predict_quantile_dd - query_quantile);

        }
        return std::make_pair(throughput, error_dd / num);
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
	std::pair<double, double> Run(double w, double tower, uint32_t threshold, uint32_t memory, uint32_t d, uint32_t mm) {
        Init();
		uint32_t run_length = 20000000;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
        DDSketch<uint64_t, uint64_t> dd(memory);
        std::vector<std::pair<uint64_t, uint64_t> > ins; ins.clear();
        
		for (int i = 0; i < run_length; ++i) {
            int g = dataset[i].id % mod;
            bool flag = 0;
            for (int j = head[g]; j; j = nxt[j])
            {
                if (tid[j] == dataset[i].id)
                {
                    flag = 1;
                    id_map[j] ++;
                    gt.insert(dataset[i].id, i - last_time[j]);

                    ins.push_back(std::make_pair(dataset[i].id, i - last_time[j]));
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

        tt = clock();

        for (int i = 0; i < ins.size(); i++)
        {
            dd.insert(ins[i].first, ins[i].second);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << totaltime << std::endl;

		double error_dd = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) {
            if (id_map[i] < 5000)
                continue;
            num++;

            uint64_t predict_dd = dd.query(tid[i], query_quantile);
            double predict_quantile_dd = gt.query(tid[i], predict_dd);
            error_dd += abs(predict_quantile_dd - query_quantile);
            //std::cout << predict_quantile_dd << "\n";
        }

        return std::make_pair(throughput, error_dd / num);
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

struct Seattle_Tuple {
    uint32_t id;
    uint32_t fetch_time;
};

std::string itos(int n)
{
    std::string res = "";
    while(n)
    {
        res = (char)('0' + (n % 10)) + res;
        n /= 10;
    }
    return res;
}

class Seattle_Benchmark {
public:
	Seattle_Benchmark() {}
	Seattle_Benchmark(std::string path) {
        int cc;
        length = 0;
        dataset = new Seattle_Tuple [n_slice * 99 * 98];
        std::string path1;
        double cur;
        for (int t = 1; t <= n_slice; t++)
        {
            path1 = path + itos(t);
            freopen(path1.c_str(), "r", stdin);
            for (int i = 1; i <= 99; i++)
            {
                for (int j = 1; j <= 99; j++)
                {
                    scanf("%lf", &cur);
                    cc = (int)(cur * 100.00) + 1;
                    if (i != j)
                    {
                        dataset[length].id = i;
                        dataset[length].fetch_time = cc;
                        length++;
                    }
                }
            }
        }
	}
	~Seattle_Benchmark()
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
		uint32_t run_length = length;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
        DDSketch<uint64_t, uint64_t> dd(memory);
        std::vector<std::pair<uint64_t, uint64_t> > ins; ins.clear();
        
		for (int i = 0; i < run_length; ++i) {
            int g = dataset[i].id % mod;
            bool flag = 0;
            for (int j = head[g]; j; j = nxt[j])
            {
                if (tid[j] == dataset[i].id)
                {
                    flag = 1;
                    id_map[j] ++;
                    gt.insert(dataset[i].id, dataset[i].fetch_time);

                    ins.push_back(std::make_pair(dataset[i].id, dataset[i].fetch_time));
                    break;
                }
            }
            if (!flag)
            {
                cnt++;
                tid[cnt] = dataset[i].id; nxt[cnt] = head[g];
                head[g] = cnt;
                id_map[cnt] = 1;
                gt.insert(dataset[i].id, dataset[i].fetch_time);

                ins.push_back(std::make_pair(dataset[i].id, dataset[i].fetch_time));
            }
		}

        tt = clock();

        for (int i = 0; i < ins.size(); i++)
        {
            dd.insert(ins[i].first, ins[i].second);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double(run_length) / totaltime;

		//std::cout << totaltime << std::endl;

		double error_dd = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) {
            if (id_map[i] < 5000)
                continue;
            num++;

            uint64_t predict_dd = dd.query(tid[i], query_quantile);
            double predict_quantile_dd = gt.query(tid[i], predict_dd);
            error_dd += abs(predict_quantile_dd - query_quantile);
            //std::cout << predict_quantile_dd << "\n";
        }

        return std::make_pair(throughput, error_dd / num);
	}

private:
    Seattle_Tuple *dataset;
    uint64_t length;
    int cnt;
	std::vector<int> head, nxt, id_map;
    std::vector<uint64_t> last_time;
	std::vector<uint64_t> tid; 
};

struct five_tuple
{
    char s[13];
    bool operator == (const five_tuple &u) const
    {
        for (int i = 0; i < 13; i++)
        {
            if (s[i] != u.s[i]) return 0;
        }
        return 1;
    }
    bool operator < (const five_tuple &u) const
    {
        for (int i = 0; i < 13; i++)
        {
            if (s[i] != u.s[i]) return s[i] < u.s[i];
        }
        return 0;
    }
    uint64_t get_hash(uint64_t md)
    {
        uint64_t res = 0;
        for (int i = 0; i < 13; i++)
        {
            res = (131ll * res + (int)s[i]) % md;
        }
        return res;
    }
};

struct Tap_Tuple {
    five_tuple id;
    uint64_t timestamp;
};


class Tap_Benchmark {
public:
	Tap_Benchmark() {}
	Tap_Benchmark(std::string path) {
		load_result = Load(path.c_str());
        dataset = (Tap_Tuple*)load_result.start;
        length = load_result.length / sizeof(Tap_Tuple);
	}
	~Tap_Benchmark() {}
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
		uint32_t run_length = 20000000;
        double query_quantile = w;
        double tottime = 0.00, tt;
		assert(run_length <= length);

		GroundTruth<uint64_t, uint64_t> gt;
        DDSketch<uint64_t, uint64_t> dd(memory);
        std::vector<std::pair<uint64_t, uint64_t> > ins; ins.clear();

		for (int i = 0; i < run_length; ++i) {
            uint64_t id = dataset[i].id.get_hash(MOD);
            int g = id % mod;
            bool flag = 0;
            for (int j = head[g]; j; j = nxt[j])
            {
                if (tid[j] == id)
                {
                    flag = 1;
                    id_map[j] ++;
                    if (dataset[i].timestamp > last_time[j])
                    {
                        gt.insert(id, dataset[i].timestamp - last_time[j]);

                        ins.push_back(std::make_pair(id, dataset[i].timestamp - last_time[j]));

                    }
                    last_time[j] = dataset[i].timestamp;
                    break;
                }
            }
            if (!flag)
            {
                cnt++;
                tid[cnt] = id; nxt[cnt] = head[g];
                head[g] = cnt;
                id_map[cnt] = 1;
                last_time[cnt] = dataset[i].timestamp;
            }
		}

        tt = clock();

        for (int i = 0; i < ins.size(); i++)
        {
            dd.insert(ins[i].first, ins[i].second);
        }

        tottime = clock() - tt;

        double totaltime = (double)(tottime) / CLOCKS_PER_SEC;
        double throughput = double((int)ins.size()) / totaltime;

		//std::cout << totaltime << std::endl;

		double error_dd = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) {
            if (id_map[i] < 5000)
                continue;
            num++;
            
            uint64_t predict_dd = dd.query(tid[i], query_quantile);
            double predict_quantile_dd = gt.query(tid[i], predict_dd);
            //printf("%.6lf\n", predict_quantile_dd);
            error_dd += abs(predict_quantile_dd - query_quantile);

        }
        return std::make_pair(throughput, error_dd / num);
	}

private:
	std::string filename;
    LoadResult load_result;
    Tap_Tuple *dataset;
    uint64_t length;
    int cnt;
	std::vector<int> head, nxt, id_map;
    std::vector<uint64_t> last_time;
	std::vector<uint64_t> tid; 
};