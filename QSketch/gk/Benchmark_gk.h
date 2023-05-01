#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#include <bits/stdc++.h>
#include <hash.h>
#include <Mmap.h>
#include "CorrectDetector.h"
#include "gk.hpp"
#include "Param.h"
#include<time.h>

#define MAXN 20000000
#define mod 10000019
#define n_slice 688
#define MOD 13908490328490391


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
        //16 = size of a tuple, 1/eps * log(eps*N) is the worst case of gk
        //according to experimental experience divided by 5.2 is more reasonable

        double temp = double(max_memory) / double( max_memory_per_bucket ); 
        bucket_num = uint64_t(floor(temp)); 

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
    std::pair<double, double> Run(uint32_t memory, double eps, double query_w) 
    {
        Init();
        uint32_t running_length = 20000000;

        compare_gk<uint64_t>* gk_sketch = new  compare_gk<uint64_t>(memory,eps,20000000);

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
            gk_sketch -> insert(ins[i].first, ins[i].second);
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
            
            uint64_t predict = gk_sketch->query(tid[i], query_w);
            double predict_quantile = correct_detector->query(tid[i], predict);
            error += fabs( predict_quantile - query_w );

        }
        //std::cout << "Average Error: " << error / num << "\n";

        //gk_sketch->print_status();

        //std::cout << "\n";

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
	std::pair<double, double> Run(uint32_t memory, double eps, double query_w) 
    {
        Init();
        uint32_t running_length = 20000000;

        compare_gk<uint64_t>* gk_sketch = new  compare_gk<uint64_t>(memory,eps,20000000);

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
            gk_sketch -> insert(ins[i].first, ins[i].second);
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
            
            uint64_t predict = gk_sketch->query(tid[i], query_w);
            double predict_quantile = correct_detector->query(tid[i], predict);
            error += fabs( predict_quantile - query_w );

        }

        //gk_sketch->print_status();

        return std::make_pair(throughput, error / num);
    }
//private:
    std::string filename;
    LoadResult load_result;
    zipf_Tuple* dataset;
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
	std::pair<double, double> Run(uint32_t memory, double eps, double query_w) 
    {
        Init();
        uint32_t running_length = length;

        compare_gk<uint64_t>* gk_sketch = new  compare_gk<uint64_t>(memory,eps,length);

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
            gk_sketch -> insert(ins[i].first, ins[i].second);
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
            
            uint64_t predict = gk_sketch->query(tid[i], query_w);
            double predict_quantile = correct_detector->query(tid[i], predict);
            error += fabs( predict_quantile - query_w );

        }

        //gk_sketch->print_status();

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
	std::pair<double, double> Run(uint32_t memory, double eps, double query_w) 
    {
        Init();
        uint32_t running_length = 20000000;

        compare_gk<uint64_t>* gk_sketch = new  compare_gk<uint64_t>(memory,eps,20000000);

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
            gk_sketch -> insert(ins[i].first, ins[i].second);
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
            
            uint64_t predict = gk_sketch->query(tid[i], query_w);
            double predict_quantile = correct_detector->query(tid[i], predict);
            error += fabs( predict_quantile - query_w );

        }
        //std::cout << "Average Error: " << error / num << "\n";

        //gk_sketch->print_status();

        //std::cout << "\n";

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