#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#include <bits/stdc++.h>
#include <linux/types.h>
#include <hash.h>
#include <Mmap.h>
#include "CorrectDetector.h"
#include "kll_sketch.hpp"
#include "Param.h"
#include<time.h>

__u64 rdtsc()
{
        __u32 lo,hi;


        __asm__ __volatile__
        (
         "rdtsc":"=a"(lo),"=d"(hi)
        );
        return (__u64)hi<<32|lo;
}


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

        std::cerr << max_memory_per_bucket << std::endl;

        bucket_num = 1;
        max_memory_per_bucket = max_memory;

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

class CAIDABenchmark 
{
public:
    CAIDABenchmark(std::string PATH) 
    {
        std::cout<<"dataset = "<<PATH<<std::endl;
        load_result = Load(PATH.c_str());
        dataset = (CAIDA_Tuple*)load_result.start;
        length = load_result.length / sizeof(CAIDA_Tuple);
    }
    ~CAIDABenchmark() {}
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
    std::pair<std::pair<double, double>, double> Run(uint32_t memory, int k_of_sketch, double query_w) 
    {
        Init();
        uint32_t running_length = 20000000;

        compare_kll<uint64_t>* KLL_sketch = new compare_kll<uint64_t>(memory,k_of_sketch,running_length);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        double total = 0.00;
        double totaltime, tt;

        std::vector<uint64_t> ins; ins.clear();

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
                        correct_detector->insert(1, dataset[i].timestamp - last_time[j]);

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
            KLL_sketch -> insert(1, ins[i]);
        }

        total = clock() - tt;

        double totaltime1 = (double)(total) / CLOCKS_PER_SEC;
        double throughput1 = double((int)ins.size()) / totaltime1;

		//std::cout << totaltime << std::endl;

		double error_qs = 0;

        __u64 begins = rdtsc();
        uint64_t predict_qs;
        for (int t = 1; t <= 1000; t++)
        {
            predict_qs = KLL_sketch -> query(1, query_w);
        }
        __u64 ends = rdtsc();
        double throughput2 = 1000000000.00 / (double)(ends - begins);
        double predict_quantile_qs = correct_detector -> query(1, predict_qs);

        //std::cout << predict_quantile_qs << std::endl;

        error_qs += fabs(predict_quantile_qs - query_w);


        return std::make_pair(std::make_pair(throughput1, throughput2), error_qs);
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
    std::pair<std::pair<double, double>, double> Run(uint32_t memory, int k_of_sketch, double query_w) {
        uint32_t running_length = 20000000;

        compare_kll<uint64_t>* KLL_sketch = new compare_kll<uint64_t>(memory,k_of_sketch,running_length);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        double total = 0.00;
        double totaltime, tt;

        std::vector<uint64_t> ins; ins.clear();

        for (int i = 0; i < running_length; ++i) 
        {
            correct_detector->insert(1, dataset[i].stamp);
            ins.push_back(dataset[i].stamp);
        }
        
        tt = clock();
        
        for (int i = 0; i < ins.size(); i++)
        {
            KLL_sketch -> insert(1, ins[i]);
        }

        total = clock() - tt;

        double totaltime1 = (double)(total) / CLOCKS_PER_SEC;
        double throughput1 = double((int)ins.size()) / totaltime1;

		//std::cout << totaltime << std::endl;

		double error_qs = 0;

        __u64 begins = rdtsc();
        uint64_t predict_qs;
        for (int t = 1; t <= 1000; t++)
        {
            predict_qs = KLL_sketch -> query(1, query_w);
        }
        __u64 ends = rdtsc();
        double throughput2 = 1000000000.00 / (double)(ends - begins);
        double predict_quantile_qs = correct_detector -> query(1, predict_qs);

        //std::cout << predict_quantile_qs << std::endl;

        error_qs += fabs(predict_quantile_qs - query_w);


        return std::make_pair(std::make_pair(throughput1, throughput2), error_qs);
	}

private:
    synthetic_Tuple *dataset;
    uint64_t length;
};

#endif
