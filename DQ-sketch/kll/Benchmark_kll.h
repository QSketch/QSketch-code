#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#include <bits/stdc++.h>
#include <hash.h>
#include <Mmap.h>
#include "CorrectDetector.h"
#include "kll_sketch.hpp"
#include "Param.h"
#include<time.h>

#define MAXN 1000000
#define mod 10000019
#define web_len 383360


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
        double totaltime;

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

                        begin=clock();
                        KLL_sketch->insert(dataset[i].id, dataset[i].timestamp - last_time[j]);
                        finish=clock();

                        total += finish-begin;
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
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double error = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) 
        {
            if (id_map[i] < 1000)
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
        double totaltime;

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

                        begin=clock();
                        KLL_sketch->insert(dataset[i].id, i - last_time[j]);
                        finish=clock();

                        total += finish-begin;
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
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double error = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) 
        {
            if (id_map[i] < 100)
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

struct tnsm_Tuple {
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
	std::pair<double, double> Run(uint32_t memory, int k_of_sketch, double query_w) {
        Init();
        uint32_t running_length = web_len;

        compare_kll<uint64_t>* KLL_sketch = new compare_kll<uint64_t>(memory,k_of_sketch,running_length);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        clock_t total=0;
        double totaltime;

        for (int i = 0; i < running_length; ++i) 
        {
            int g = dataset[i].unit_id % mod;
            bool flag = 0;
            for (int j = head[g]; j; j = nxt[j])
            {
                if (tid[j] == dataset[i].unit_id)
                {
                    flag = 1;
                    id_map[j] ++;
                    correct_detector->insert(dataset[i].unit_id, dataset[i].fetch_time);

                    begin=clock();
                    KLL_sketch->insert(dataset[i].unit_id, dataset[i].fetch_time);
                    finish=clock();

                    total += finish-begin;
                    break;
                }
            }
            if (!flag)
            {
                cnt++;
                tid[cnt] = dataset[i].unit_id; nxt[cnt] = head[g];
                head[g] = cnt;
                id_map[cnt] = 1;

                correct_detector->insert(dataset[i].unit_id, dataset[i].fetch_time);

                begin=clock();
                KLL_sketch->insert(dataset[i].unit_id, dataset[i].fetch_time);
                finish=clock();

                total += finish-begin;
            }

        }
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double error = 0;
        int num = 0;
        for (int i = 1; i <= cnt; i++) 
        {
            if (id_map[i] < 1000)
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
    tnsm_Tuple *dataset;
    uint64_t length;
    int cnt;
	std::vector<int> head, nxt, id_map;
    std::vector<uint64_t> last_time;
	std::vector<uint64_t> tid; 
};

#endif
