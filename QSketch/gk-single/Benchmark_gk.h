#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#include <bits/stdc++.h>
#include <hash.h>
#include <Mmap.h>
#include "CorrectDetector.h"
#include "gk.hpp"
#include "Param.h"
#include<time.h>


template<typename ID_TYPE>
class compare_gk
{
public:
    uint64_t max_memory;
    double epsilon;
    uint32_t m_one_over_2e;
    uint64_t max_memory_per_bucket;
    uint64_t bucket_num;
    uint64_t item_inserted;
    
    std::vector<gk<uint64_t>> array_gk;
    
    compare_gk(uint64_t mem, double eps, uint64_t total_item)
    {
        max_memory=mem*1024;
        epsilon=eps;

        m_one_over_2e = uint32_t( round((1.0/(2.0 * epsilon))) );
        max_memory_per_bucket = uint64_t(round(1.0/epsilon * log(eps * double(total_item)) /5.2 )) * 16; 
        std::cerr << max_memory_per_bucket << std::endl;
        
        //16 = size of a tuple, 1/eps * log(eps*N) is the worst case of gk
        //according to experimental experience divided by 5.2 is more reasonable

        bucket_num = 1; 
        max_memory_per_bucket = max_memory; 

        for (int i=0;i<bucket_num;++i)
            array_gk.push_back(gk<uint64_t>(epsilon));

        item_inserted = 0;
    } 

    void insert(ID_TYPE id, uint64_t timestamp)
    {
        uint32_t index = hash(id, 1024) % array_gk.size();
        array_gk[index].insert(timestamp);
        item_inserted++;
    }
    uint64_t query(ID_TYPE id, double w)
    {
        uint32_t index = hash(id, 1024) % array_gk.size();
        return array_gk[index].quantile(w);
    }
    uint32_t get_index(ID_TYPE id)
    {
        return hash(id, 1024) % array_gk.size();
    }
    uint32_t actual_len()
    {
        uint64_t res=0;
        for (int i=0;i<bucket_num;++i)
        {
            //printf("len of bucket %d: %ld\n", i, array_gk[i].m_S.size());
            res += array_gk[i].m_S.size();
        }
        return res;
    }
    void print_status()
    {
        uint32_t res=0;
        for (int i=0;i<bucket_num;++i)
        {
            //printf("len of bucket %d: %ld\n", i, array_gk[i].m_S.size());
            res += array_gk[i].m_S.size();
        }
        
        printf("-------compare_gk status------------\n");

        printf("----max_memory            : %lu bytes = %lu Kb\n", max_memory, max_memory/1024);
        printf("----epsilon               : %f \n",epsilon);
        printf("----m_one_over_2e         : %u \n",m_one_over_2e);
        printf("----max_memory_per_bucket : %lu bytes = %lu Kb\n", max_memory_per_bucket, max_memory_per_bucket/1024);
        printf("----bucket_num            : %lu \n",bucket_num);
        printf("----item_inserted         : %lu \n",item_inserted);
        printf("----total len             : %u\n", res);
        //printf("----ratio                 : %f\n", double(res)/double(bucket_num)/double(m_one_over_2e));
        printf("----total size            : %u bytes = %u Kb\n", res*16, res*16/1024);
        printf("----nominal size - actual_size = %ld bytes = %ld Kb\n", max_memory-res*16, (max_memory-res*16)/1024);

        printf("-------compare_gk status end--------\n");
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
    std::pair<double, double> Run(uint32_t memory, double eps, double query_w) 
    {
        Init();
        uint32_t running_length = 20000000;

        compare_gk<uint64_t>* gk_sketch = new  compare_gk<uint64_t>(memory,eps,20000000);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        clock_t total=0;
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
            gk_sketch -> insert(1, ins[i]);
        }

        total = clock() - tt;
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double log_error = 0;
        uint64_t predict = gk_sketch->query(1, query_w);
        double predict_quantile = correct_detector->query(1, predict);
        log_error += fabs( predict_quantile - query_w );
        //std::cout << "Average Error: " << log_error / num << "\n";

        //gk_sketch->print_status();

        //std::cout << "\n";

        return std::make_pair(throughput, log_error);
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
    std::pair<double, double> Run(uint32_t memory, double eps, double query_w) {
        uint32_t running_length = 20000000;

        compare_gk<uint64_t>* gk_sketch = new  compare_gk<uint64_t>(memory,eps,20000000);

        CorrectDetector<uint64_t, uint64_t>* correct_detector = new CorrectDetector<uint64_t, uint64_t>(); 

        clock_t begin,finish;
        clock_t total=0;
        double totaltime, tt;

        std::vector<uint64_t> ins; ins.clear();
        for (int i = 0; i < running_length; ++i) 
        {
            ins.push_back(dataset[i].stamp);
            correct_detector->insert(1, dataset[i].stamp);
        }

        tt = clock();
        
        for (int i = 0; i < ins.size(); i++)
        {
            gk_sketch -> insert(1, ins[i]);
        }

        total = clock() - tt;
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double log_error = 0;
        uint64_t predict = gk_sketch->query(1, query_w);
        double predict_quantile = correct_detector->query(1, predict);
        log_error += fabs( predict_quantile - query_w );
        //std::cout << "Average Error: " << log_error / num << "\n";

        //gk_sketch->print_status();

        //std::cout << "\n";

        return std::make_pair(throughput, log_error);
	}

private:
    synthetic_Tuple *dataset;
    uint64_t length;
};

#endif
