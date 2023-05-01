#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#include <bits/stdc++.h>
#include <hash.h>
#include <Mmap.h>
#include "CorrectDetector.h"
#include "kll_sketch.hpp"
#include "Param.h"
#include<time.h>

#define MAXN 20000000
#define mod 10000019
#define n_slice 688
#define MOD 13908490328490391


template<typename ID_TYPE>
class compare_kll
{
public:
    uint64_t max_memory;
    uint32_t k_of_sketch;
    uint64_t max_memory_per_bucket;
    uint64_t bucket_num;
    uint64_t item_inserted;
    
    std::vector<datasketches::kll_sketch<uint64_t>> array_kll;
    
    compare_kll(uint64_t mem, uint32_t KK, uint64_t total_item)
    {
        max_memory=mem*1024;
        k_of_sketch=KK;

        {
            datasketches::kll_sketch<uint64_t> sketch_temp(k_of_sketch);
            for (int i=0;i<total_item;++i)
                sketch_temp.update(i);
            max_memory_per_bucket = sketch_temp.get_serialized_size_bytes();
    
            double temp = double(max_memory) / double( max_memory_per_bucket ); 
            bucket_num = uint64_t(floor(temp));
        }
        
        //compute max_memory_per_bucket and bucket_num again. a bucket does not need to insert all item
        {
            uint64_t new_total_item = total_item/bucket_num;
            datasketches::kll_sketch<uint64_t> sketch_temp1(k_of_sketch);
            for (int i=0;i<new_total_item;++i)
                sketch_temp1.update(i);
            max_memory_per_bucket = sketch_temp1.get_serialized_size_bytes();
    
            double temp1 = double(max_memory) / double( max_memory_per_bucket ); 
            bucket_num = uint64_t(floor(temp1));
        }

        for (int i=0;i<bucket_num;++i)
            array_kll.push_back(datasketches::kll_sketch<uint64_t>(k_of_sketch));

        item_inserted = 0;
    } 
    
    void insert(ID_TYPE id, uint64_t timestamp)
    {
        uint32_t index = hash(id, 1024) % array_kll.size();
        array_kll[index].update(timestamp);
        item_inserted++;
    }
    uint64_t query(ID_TYPE id, double w)
    {
        uint32_t index = hash(id, 1024) % array_kll.size();

        const double fractions[1] {w};
        auto quantiles = array_kll[index].get_quantiles(fractions, 1);

        return quantiles[0];
    }
    uint32_t get_index(ID_TYPE id)
    {
        return hash(id, 1024) % array_kll.size();
    }
    uint32_t actual_len()
    {
        uint64_t res=0;
        for (int i=0;i<bucket_num;++i)
        {
            res += array_kll[i].get_serialized_size_bytes();
        }
        return res;
    }
    void print_status()
    {
        uint32_t res=0;
        for (int i=0;i<bucket_num;++i)
        {
            res += array_kll[i].get_serialized_size_bytes();
        }
        
        printf("-------compare_kll status------------\n");

        printf("----max_memory            : %lu bytes = %lu Kb\n", max_memory, max_memory/1024);
        printf("----k_of_sketch           : %u \n",k_of_sketch);
        printf("----max_memory_per_bucket : %lu bytes = %lu Kb\n", max_memory_per_bucket, max_memory_per_bucket/1024);
        printf("----bucket_num            : %lu \n",bucket_num);
        printf("----item_inserted         : %lu \n",item_inserted);
        printf("----total size            : %u bytes = %u Kb\n", res, res/1024);
        printf("----nominal size - actual_size = %ld bytes = %ld Kb\n", max_memory-res, (max_memory-res)/1024);

        printf("-------compare_kll status end--------\n");
    }
};

struct CAIDA_Tuple {
    uint64_t timestamp;
    uint64_t id;
};

class CAIDA_Benchmark 
{
public:
    CAIDA_Benchmark(std::string PATH) 
    {
        //std::cout<<"dataset = "<<PATH<<std::endl;
        load_result = Load(PATH.c_str());
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
    std::pair<double,double> Run(uint32_t memory, int k_of_sketch, double query_w) 
    {
        Init();
        uint32_t running_length = 20000000;

        compare_kll<uint64_t>* KLL_sketch = new compare_kll<uint64_t>(memory,k_of_sketch,running_length);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        clock_t total=0;
        double totaltime, tt;

        std::vector<std::pair<uint64_t, uint64_t> > ins; ins.clear();


        for (int i = 0; i < running_length; ++i) 
        {
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
                        correct_detector->insert(dataset[i].id, dataset[i].timestamp - last_time[j]);

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
            KLL_sketch -> insert(ins[i].first, ins[i].second);
        }

        total = clock() - tt;
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double error = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) 
        {
            if (id_map[i] < 5000)
                continue;
            num++;
            
            uint64_t predict = KLL_sketch->query(tid[i], query_w);
            double predict_quantile = correct_detector->query(tid[i], predict);
            error += fabs( predict_quantile - query_w );

        }

        return std::make_pair(throughput, error / num);
    }
//private:
    std::string filename;
    LoadResult load_result;
    CAIDA_Tuple* dataset;
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
	std::pair<double, double> Run(uint32_t memory, int k_of_sketch, double query_w) 
    {
        Init();
        uint32_t running_length = 20000000;

        compare_kll<uint64_t>* KLL_sketch = new compare_kll<uint64_t>(memory,k_of_sketch,running_length);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        clock_t total=0;
        double totaltime, tt;

        std::vector<std::pair<uint64_t, uint64_t> > ins; ins.clear();

        for (int i = 0; i < running_length; ++i) 
        {
            int g = dataset[i].id % mod;
            bool flag = 0;
            for (int j = head[g]; j; j = nxt[j])
            {
                if (tid[j] == dataset[i].id)
                {
                    flag = 1;
                    id_map[j] ++;
                    if (i > last_time[j])
                    {
                        correct_detector->insert(dataset[i].id, i - last_time[j]);

                        ins.push_back(std::make_pair(dataset[i].id, i - last_time[j]));
                    }
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
            KLL_sketch -> insert(ins[i].first, ins[i].second);
        }

        total = clock() - tt;
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double error = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) 
        {
            if (id_map[i] < 5000)
                continue;
            num++;
            
            uint64_t predict = KLL_sketch->query(tid[i], query_w);
            double predict_quantile = correct_detector->query(tid[i], predict);
            error += fabs( predict_quantile - query_w );

        }

        return std::make_pair(throughput, error / num);
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
                    cc = (int)(cur * 100.00);
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
	std::pair<double, double> Run(uint32_t memory, int k_of_sketch, double query_w) 
    {
        Init();
        uint32_t running_length = length;

        compare_kll<uint64_t>* KLL_sketch = new compare_kll<uint64_t>(memory,k_of_sketch,running_length);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        clock_t total=0;
        double totaltime, tt;

        std::vector<std::pair<uint64_t, uint64_t> > ins; ins.clear();

        for (int i = 0; i < running_length; ++i) 
        {
            int g = dataset[i].id % mod;
            bool flag = 0;
            for (int j = head[g]; j; j = nxt[j])
            {
                if (tid[j] == dataset[i].id)
                {
                    flag = 1;
                    id_map[j] ++;
                    correct_detector->insert(dataset[i].id, dataset[i].fetch_time);

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
                last_time[cnt] = i;
                correct_detector->insert(dataset[i].id, dataset[i].fetch_time);

                ins.push_back(std::make_pair(dataset[i].id, dataset[i].fetch_time));
            }

        }

        tt = clock();

        for (int i = 0; i < ins.size(); i++)
        {
            KLL_sketch -> insert(ins[i].first, ins[i].second);
        }

        total = clock() - tt;
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double error = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) 
        {
            if (id_map[i] < 5000)
                continue;
            num++;
            
            uint64_t predict = KLL_sketch->query(tid[i], query_w);
            double predict_quantile = correct_detector->query(tid[i], predict);
            error += fabs( predict_quantile - query_w );

        }

        return std::make_pair(throughput, error / num);
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
	std::pair<double,double> Run(uint32_t memory, int k_of_sketch, double query_w) 
    {
        Init();
        uint32_t running_length = 20000000;

        compare_kll<uint64_t>* KLL_sketch = new compare_kll<uint64_t>(memory,k_of_sketch,running_length);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        clock_t total=0;
        double totaltime, tt;

        std::vector<std::pair<uint64_t, uint64_t> > ins; ins.clear();


        for (int i = 0; i < running_length; ++i) 
        {
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
                        correct_detector->insert(id, dataset[i].timestamp - last_time[j]);

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
            KLL_sketch -> insert(ins[i].first, ins[i].second);
        }

        total = clock() - tt;
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double((int)ins.size()) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double error = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) 
        {
            if (id_map[i] < 5000)
                continue;
            num++;
            
            uint64_t predict = KLL_sketch->query(tid[i], query_w);
            double predict_quantile = correct_detector->query(tid[i], predict);
            error += fabs( predict_quantile - query_w );

        }

        return std::make_pair(throughput, error / num);
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

#endif
