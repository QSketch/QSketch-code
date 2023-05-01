#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#include <bits/stdc++.h>
#include <hash.h>
#include <Mmap.h>
#include "CorrectDetector.h"
#include "Param.h"
#include <time.h>
#include "testSquadGk.h"
#define web_len 383360
#define n_slice 688

/*
struct Seattle_Tuple {
    uint64_t timestamp;
    uint64_t id;
};

class SeattleBenchmark 
{
public:

    SeattleBenchmark(std::string PATH) 
    {
        std::cout<<"dataset = "<<PATH<<std::endl;

        int datalen=0;
        std::vector<double> num;

        std::ifstream file(PATH.c_str());
        while( ! file.eof() )
        {
            double temp;
            file>>temp;
            num.push_back(temp*100000);
            datalen++;
        }
        file.close();
        std::cout<<datalen<<" "<<num.size()<<"\n";

        //seattle: 688file, each file has 99*99 item

        dataset = new Seattle_Tuple[datalen];
        int index = 0;
        for (int filenum=0;filenum<688;filenum++)
        {
            for (int flow_id=0;flow_id<99;flow_id++)
            {
                for (int i=0;i<99;i++)
                {
                    dataset[index].id=index % (9801) * 41 + 41;
                    dataset[index].timestamp=num[index];
                    index++;
                }   
                  
            }
        }
        std::cout<<"index="<<index<<", datalen="<<datalen<<std::endl;
        //printf("dataset_len=%d, index=%d, large=%d, mid=%d, small=%d\n",datalen,index,large_stream,mid_stream,small_stream);
        
        length = datalen;

    }
    ~SeattleBenchmark() {}

    void Run(uint32_t memory) 
    {
        uint32_t running_length = length;

        IcedSketch<uint64_t>* icde_sketch = new IcedSketch<uint64_t>(memory);

        CorrectDetector<uint64_t>* correct_detector = new CorrectDetector<uint64_t>(); 
        std::map<uint64_t, uint32_t> id_map;
        std::set<uint64_t> id_set;

        clock_t begin,finish;
        clock_t total=0;
        double totaltime;
        uint32_t actual_length=0;

        for (int i = 0; i < running_length; ++i) 
        {
            if (i % (1000000) ==0) printf("i = %d\n",i);
            if (dataset[i].timestamp == 0) continue; //把自己到自己的rtt去掉

            actual_length++;

            if (id_set.find(dataset[i].id) == id_set.end()) 
            {
                id_set.insert(dataset[i].id);
                id_map[dataset[i].id] = 0;
            }
            id_map[dataset[i].id]++;

            #ifdef TIME_BASED
            if (id_map[dataset[i].id] > 1) 
            {
                correct_detector->insert(dataset[i].id, 0, dataset[i].timestamp);

                begin=clock();
                icde_sketch->insert(dataset[i].id, 0, dataset[i].timestamp);
                finish=clock();

                total += finish-begin;
            }
            #else
            if (id_map[dataset[i].id] > 1) 
            {
                //correct_detector->insert(dataset[i].id, last_time[dataset[i].id], i);
                //icde_sketch->insert(dataset[i].id, last_time[dataset[i].id], i);
            }
            #endif
        }
        
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        std::cout <<"time taken: "<<totaltime<<" seconds"<< std::endl;
        double throughput = double(actual_length) / totaltime;
        std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double log_error = 0;
        int num = 0;
        double query_w = 0.95;
        for (auto i : id_map) 
        {
            
            if (i.second < 500)//1000可以改大，例如2000，不要出现1000 和 2000效果差别巨大的情况
                continue;
            num++;

            uint64_t predict = icde_sketch->query(i.first, query_w);
            uint64_t truth = correct_detector->query(i.first, query_w);
            log_error += fabs( log2(predict) - log2(truth) );
            
            //std::cout<<"id="<<i.first<<" "<<"num="<<i.second<<std::endl;
            // std::cout << log2(predict) << " " << log2(truth) << " " << "\n";
            // assert(0);
            // printf("id=%lu, predict=%lu, truth=%lu\n",i.first,predict,truth);
            // if (num==1) break;
        }
        std::cout << "Estimate: " << num << "\n";
        std::cout << "Average Log Error: " << log_error / num << "\n";
        std::cout << "query_w: " << query_w<< "\n";

        icde_sketch->print_status();

        std::cout << "\n";

    }

//private:
    std::string filename;
    LoadResult load_result;
    Seattle_Tuple* dataset;
    uint64_t length;
};


struct Webget_Tuple {
    uint64_t id;
    uint64_t timestamp;
};

class WebgetBenchmark 
{
public:
    WebgetBenchmark(std::string PATH) 
    {
        
        std::cout<<"dataset = "<<PATH<<std::endl;

        int datalen=0;
        std::vector<uint64_t> num1;
        std::vector<uint64_t> num2;

        for (int i=0;i<2;++i)
        {
            std::ifstream file(PATH.c_str());
            while( ! file.eof() )
            {
                uint64_t temp1;
                uint64_t temp2;
                file>>temp1>>temp2;
                if (temp2==0) continue;
    
                num1.push_back(temp1);
                //num.push_back(temp2 * 123465);
                num2.push_back(temp2);
                
                datalen+=1;
            }
            file.close();
        }
        
        std::cout<<datalen<<" "<<num1.size()<<" "<<num2.size()<<"\n";

        //Webget: 9368 flow,some timestamp is 0

        dataset = new Webget_Tuple[datalen];
        for (int i=0;i<datalen;i++)
        {
            dataset[i].id=num1[i];
            dataset[i].timestamp=num2[i];
        }
        
        length = datalen;
    }
    ~WebgetBenchmark() {}

    void Run(uint32_t memory) 
    {
        uint32_t running_length = length;

        IcedSketch<uint64_t>* icde_sketch = new IcedSketch<uint64_t>(memory);

        CorrectDetector<uint64_t>* correct_detector = new CorrectDetector<uint64_t>(); 

        std::map<uint64_t, uint32_t> id_map;
        std::set<uint64_t> id_set;

        clock_t begin,finish;
        clock_t total=0;
        double totaltime;
        uint32_t actual_length=0;

        for (int i = 0; i < running_length; ++i) 
        {
            if (i % (1000000) ==0) printf("i = %d\n",i);
            if (dataset[i].timestamp == 0) continue;

            actual_length++;

            if (id_set.find(dataset[i].id) == id_set.end()) 
            {
                id_set.insert(dataset[i].id);
                id_map[dataset[i].id] = 0;
            }
            id_map[dataset[i].id]++;

            #ifdef TIME_BASED
            if (id_map[dataset[i].id] > 1) 
            {
                correct_detector->insert(dataset[i].id, 0, dataset[i].timestamp);
                
                begin=clock();
                icde_sketch->insert(dataset[i].id, 0, dataset[i].timestamp);
                finish=clock();

                total += finish-begin;
            }
            #else
            if (id_map[dataset[i].id] > 1) 
            {
                //correct_detector->insert(dataset[i].id, last_time[dataset[i].id], i);
                //icde_sketch->insert(dataset[i].id, last_time[dataset[i].id], i);
            }
            #endif
        }

        //std::cout<<"id_map[1722174100]: "<<id_map[1722174100]<<std::endl;
        //std::cout<<"actual_length: "<<actual_length<<std::endl;
        totaltime=(double)(total)/CLOCKS_PER_SEC;
        std::cout <<"time taken: "<<totaltime<<" seconds"<< std::endl;
        double throughput = double(actual_length) / totaltime;
        std::cout <<"throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        double log_error = 0;
        int num = 0;
        double query_w = 0.95;
        for (auto i : id_map) 
        {
            if (i.second < 1000)//1000可以改大，例如2000，不要出现1000 和 2000效果差别巨大的情况
                continue;
            num++;

            uint64_t predict = icde_sketch->query(i.first, query_w);
            uint64_t truth = correct_detector->query(i.first, query_w);
            log_error += fabs( log2(predict) - log2(truth) );

            //std::cout << log2(predict) << " " << log2(truth) << " " << "\n";
            // assert(0);
            //printf("id=%lu, predict=%lu, truth=%lu, sub=%lf \n",i.first,predict,truth,fabs( log2(predict) - log2(truth) ));

            // if (num==20) break;
        }
        std::cout << "Estimate: " << num << "\n";
        std::cout << "Average Log Error: " << log_error / num << "\n";
        std::cout << "query_w: " << query_w<< "\n";

        icde_sketch->print_status();

        std::cout << "\n";

    }
//private:
    std::string filename;
    LoadResult load_result;
    Webget_Tuple* dataset;
    uint64_t length;
};
*/

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
        std::cout<<"lenght="<<length<<"\n";
    }
    ~CAIDABenchmark() {}

    uint64_t getLowBits(uint64_t num, int n) {
        uint64_t mask = (1 << n) - 1;
        return num & mask;
    }

    pair<double, double> Run(uint32_t memory, double eps, double quantile, double theta) 
    {
        uint32_t running_length = 20000000;

        //保证在1-delta的概率 查询误差在eps以内

        reqquantile = quantile;
        //theta = 0.00006; //1000 / 2000_0000
        //theta = 0.0001;
        std::cout << "theta: " << theta << "\n";

        epsilon = eps;       
        long long internal_counters = (int) ((double) (pow(theta, -1) * (double) pow(epsilon,-0.5)));
        long long samples_num = 4* ((int) ((double) (pow(theta, -2) * (double) pow(epsilon,-2)) * pow((double) internal_counters, -1)));
        double qsketch_eps = epsilon / (double) 2.0;
        
        LCU_type* squad = LCU_Init((double) 1/ (double) internal_counters, qsketch_eps, samples_num);
        double size_squad = LCU_Size(squad)/1024;
        std::cout<<"size: "<<size_squad<<" kb"<<std::endl;
        CorrectDetector<uint64_t>* correct_detector = new CorrectDetector<uint64_t>(); 
        std::map<uint64_t, uint32_t> id_map;
        std::map<uint64_t, uint64_t> last_time;
        std::set<uint64_t> id_set;

        clock_t begin,finish;
        clock_t total=0;
        double totaltime;

        std::vector<pair<uint32_t, double> > v; v.clear();

        for (int i = 0; i < running_length; ++i) 
        {
            uint32_t id = dataset[i].id & ((1 << 30) - 1);

            if (i%2000000==0) printf("i = %d\n",i);
            if (id_set.find(id) == id_set.end()) 
            {
                id_set.insert(id);
                id_map[id] = 0;
            }
            id_map[id]++;

            #ifdef TIME_BASED
            if (id_map[id] > 1 && dataset[i].timestamp > last_time[id]) 
            {
                uint64_t latency = dataset[i].timestamp - last_time[id];
                latency >>= 20;

                correct_detector->insert(id, 0, latency);
                //if (latency == 0) std::cout<<i<<"\n";

                v.push_back(std::make_pair(id, (double)latency));

            }
            last_time[id] = dataset[i].timestamp;
            #else
            if (id_map[id] > 1) 
            {
                correct_detector->insert(id, last_time[id], i);
                //icde_sketch->insert(id, last_time[id], i);
            }
            last_time[id] = i;
            #endif
        }

        begin = clock();
        for (int i = 0; i < v.size(); i++)
        {
            LCU_UpdateLatency(squad, v[i].first, v[i].second);
        }
        total = clock() - begin;


        totaltime=(double)(total)/CLOCKS_PER_SEC;
        //std::cout <<"insert time taken: "<<totaltime<<" seconds"<< std::endl;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"insert throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;
        
        double error = 0;
        int num = 0;
        double query_w = quantile;        
        total = 0;
        uint64_t total_id = 0;
        uint32_t zero_p = 0;
        for (auto i : id_map) 
        {
            total_id++;
            if (i.second < 1000)//1000可以改大，例如2000，不要出现1000 和 2000效果差别巨大的情况
                continue;
            num++;


            begin = clock();
            double est = LCU_QuantileQuery(squad, (i.first & ((1 << 30) - 1)), query_w); //是latency
            finish = clock();
            total += finish-begin;

            double predict = correct_detector->query_quantile((i.first & ((1 << 30) - 1)), (uint64_t)est);
            error += abs( predict - query_w );

            //std::cout << predict << std::endl;

            //if (predict == 0) zero_p++;
            /*
            std::cout << log2(predict) << " " << log2(truth) << " " << "\n";
            // assert(0);
            
            printf("id=%llu, predict=%llu, truth=%llu, abs=%f, fabs=%f \n",i.first,predict,truth,
                                                                        abs( log2(predict) - log2(truth) ),
                                                                        fabs( log2(predict) - log2(truth) ));
                                                                        */
            //if (num==20) break;
            
            
        }
        /*totaltime=(double)(total)/CLOCKS_PER_SEC;
        std::cout <<"query time taken: "<<totaltime<<" seconds"<< std::endl;
        throughput = double(num) / totaltime;
        std::cout <<"query throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        std::cout << "Total id: " << total_id << "\n";
        std::cout << "Estimate: " << num << "\n";
        std::cout << "Average Log Error: " << log_error / num << "\n";
        std::cout << "query_w: " << query_w<< "\n";
        std::cout << "predice = 0: " << zero_p << "\n";*/

        return std::make_pair(throughput, error / num);

    }
//private:
    std::string filename;
    LoadResult load_result;
    CAIDA_Tuple* dataset;
    uint64_t length;
};

struct zipf_Tuple {
    uint32_t id;
};

class zipf_Benchmark 
{
public:
    zipf_Benchmark(std::string PATH) 
    {
        std::cout<<"dataset = "<<PATH<<std::endl;
        load_result = Load(PATH.c_str());
        dataset = (zipf_Tuple*)load_result.start;
        length = load_result.length / sizeof(zipf_Tuple);
        std::cout<<"lenght="<<length<<"\n";
    }
    ~zipf_Benchmark() {}

    uint64_t getLowBits(uint64_t num, int n) {
        uint64_t mask = (1 << n) - 1;
        return num & mask;
    }

    pair<double, double> Run(uint32_t memory, double eps, double quantile, double theta) 
    {
        uint32_t running_length = 20000000;

        //保证在1-delta的概率 查询误差在eps以内

        reqquantile = quantile;
        //theta = 0.00006; //1000 / 2000_0000
        //theta = 0.0001;
        std::cout << "theta: " << theta << "\n";

        epsilon = eps;       
        long long internal_counters = (int) ((double) (pow(theta, -1) * (double) pow(epsilon,-0.5)));
        long long samples_num = 4* ((int) ((double) (pow(theta, -2) * (double) pow(epsilon,-2)) * pow((double) internal_counters, -1)));
        double qsketch_eps = epsilon / (double) 2.0;
        
        LCU_type* squad = LCU_Init((double) 1/ (double) internal_counters, qsketch_eps, samples_num);
        double size_squad = LCU_Size(squad)/1024;
        std::cout<<"size: "<<size_squad<<" kb"<<std::endl;
        CorrectDetector<uint64_t>* correct_detector = new CorrectDetector<uint64_t>(); 
        std::map<uint64_t, uint32_t> id_map;
        std::map<uint64_t, uint64_t> last_time;
        std::set<uint64_t> id_set;

        clock_t begin,finish;
        clock_t total=0;
        double totaltime;

        std::vector<pair<uint32_t, double> > v; v.clear();

        for (int i = 0; i < running_length; ++i) 
        {
            uint32_t id = dataset[i].id;

            if (i%2000000==0) printf("i = %d\n",i);
            if (id_set.find(id) == id_set.end()) 
            {
                id_set.insert(id);
                id_map[id] = 0;
            }
            id_map[id]++;

            if (id_map[id] > 1 && i > last_time[id]) 
            {
                uint32_t latency = i - last_time[id];

                correct_detector->insert(id, 0, latency);
                //if (latency == 0) std::cout<<i<<"\n";

                v.push_back(std::make_pair(id, (double)latency));

            }
            last_time[id] = i;

        }

        std::cout << v.size() << std::endl;

        begin = clock();
        for (int i = 0; i < v.size(); i++)
        {
            LCU_UpdateLatency(squad, v[i].first, v[i].second);
        }
        total = clock() - begin;


        totaltime=(double)(total)/CLOCKS_PER_SEC;
        //std::cout <<"insert time taken: "<<totaltime<<" seconds"<< std::endl;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"insert throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;
        
        double error = 0;
        int num = 0;
        double query_w = quantile;        
        total = 0;
        uint64_t total_id = 0;
        uint32_t zero_p = 0;
        for (auto i : id_map) 
        {
            total_id++;
            if (i.second < 1000)//1000可以改大，例如2000，不要出现1000 和 2000效果差别巨大的情况
                continue;
            num++;


            begin = clock();
            double est = LCU_QuantileQuery(squad, i.first, query_w); //是latency
            finish = clock();
            total += finish-begin;

            double predict = correct_detector->query_quantile(i.first, (uint32_t)est);
            error += abs( predict - query_w );

            //if (predict == 0) zero_p++;
            /*
            std::cout << log2(predict) << " " << log2(truth) << " " << "\n";
            // assert(0);
            
            printf("id=%llu, predict=%llu, truth=%llu, abs=%f, fabs=%f \n",i.first,predict,truth,
                                                                        abs( log2(predict) - log2(truth) ),
                                                                        fabs( log2(predict) - log2(truth) ));
            if (num==20) break;
            */
            
        }
        /*totaltime=(double)(total)/CLOCKS_PER_SEC;
        std::cout <<"query time taken: "<<totaltime<<" seconds"<< std::endl;
        throughput = double(num) / totaltime;
        std::cout <<"query throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        std::cout << "Total id: " << total_id << "\n";
        std::cout << "Estimate: " << num << "\n";
        std::cout << "Average Log Error: " << log_error / num << "\n";
        std::cout << "query_w: " << query_w<< "\n";
        std::cout << "predice = 0: " << zero_p << "\n";*/

        return std::make_pair(throughput, error / num);

    }
//private:
    std::string filename;
    LoadResult load_result;
    zipf_Tuple* dataset;
    uint64_t length;
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
        /*for (int i = 0; i < length; i++)
        {
            std::cout << dataset[i].id << " " << dataset[i].fetch_time << std::endl;
        }*/
	}
	~Seattle_Benchmark()
    {
        delete dataset;
    }
	pair<double, double> Run(uint32_t memory, double eps, double quantile, double theta) 
    {
        uint32_t running_length = length;

        //保证在1-delta的概率 查询误差在eps以内

        reqquantile = quantile;
        //theta = 0.00006; //1000 / 2000_0000
        //theta = 0.0001;
        theta = 1000.00 / (double)length;
        std::cout << "theta: " << theta << "\n";

        epsilon = eps;       
        long long internal_counters = (int) ((double) (pow(theta, -1) * (double) pow(epsilon,-0.5)));
        long long samples_num = 4* ((int) ((double) (pow(theta, -2) * (double) pow(epsilon,-2)) * pow((double) internal_counters, -1)));
        double qsketch_eps = epsilon / (double) 2.0;
        
        LCU_type* squad = LCU_Init((double) 1/ (double) internal_counters, qsketch_eps, samples_num);
        double size_squad = LCU_Size(squad)/1024;
        std::cout<<"size: "<<size_squad<<" kb"<<std::endl;
        CorrectDetector<uint64_t>* correct_detector = new CorrectDetector<uint64_t>(); 
        std::map<uint64_t, uint32_t> id_map;
        std::map<uint64_t, uint64_t> last_time;
        std::set<uint64_t> id_set;

        clock_t begin,finish;
        clock_t total=0;
        double totaltime;

        std::vector<pair<uint32_t, double> > v; v.clear();

        for (int i = 0; i < running_length; ++i) 
        {
            uint32_t id = dataset[i].id;

            if (i%2000000==0) printf("i = %d\n",i);
            if (id_set.find(id) == id_set.end()) 
            {
                id_set.insert(id);
                id_map[id] = 0;
            }
            id_map[id]++;

            //std::cout << id << " " << dataset[i].fetch_time << std::endl;

            correct_detector->insert(id, 0, dataset[i].fetch_time);
            //if (latency == 0) std::cout<<i<<"\n";

            v.push_back(std::make_pair(id, (double)dataset[i].fetch_time));

        }

        //std::cout << v.size() << std::endl;

        begin = clock();
        for (int i = 0; i < v.size(); i++)
        {
            LCU_UpdateLatency(squad, v[i].first, v[i].second);
        }
        total = clock() - begin;


        totaltime=(double)(total)/CLOCKS_PER_SEC;
        //std::cout <<"insert time taken: "<<totaltime<<" seconds"<< std::endl;
        double throughput = double(running_length) / totaltime;
        //std::cout <<"insert throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;
        
        double error = 0;
        int num = 0;
        double query_w = quantile;        
        total = 0;
        uint64_t total_id = 0;
        uint32_t zero_p = 0;
        for (auto i : id_map) 
        {
            total_id++;
            if (i.second < 1000)//1000可以改大，例如2000，不要出现1000 和 2000效果差别巨大的情况
                continue;
            num++;


            begin = clock();
            double est = LCU_QuantileQuery(squad, i.first, query_w); //是latency
            finish = clock();
            total += finish-begin;

            double predict = correct_detector->query_quantile(i.first, (uint32_t)est);
            error += abs( predict - query_w );

            //if (predict == 0) zero_p++;
            /*
            std::cout << log2(predict) << " " << log2(truth) << " " << "\n";
            // assert(0);
            
            printf("id=%llu, predict=%llu, truth=%llu, abs=%f, fabs=%f \n",i.first,predict,truth,
                                                                        abs( log2(predict) - log2(truth) ),
                                                                        fabs( log2(predict) - log2(truth) ));
            if (num==20) break;
            */
            
        }
        /*totaltime=(double)(total)/CLOCKS_PER_SEC;
        std::cout <<"query time taken: "<<totaltime<<" seconds"<< std::endl;
        throughput = double(num) / totaltime;
        std::cout <<"query throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        std::cout << "Total id: " << total_id << "\n";
        std::cout << "Estimate: " << num << "\n";
        std::cout << "Average Log Error: " << log_error / num << "\n";
        std::cout << "query_w: " << query_w<< "\n";
        std::cout << "predice = 0: " << zero_p << "\n";*/

        return std::make_pair(throughput, error / num);

    }

private:
    std::string filename;
    LoadResult load_result;
    Seattle_Tuple* dataset;
    uint64_t length;
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
	pair<double, double> Run(uint32_t memory, double eps, double quantile, double theta) 
    {
        uint32_t running_length = 20000000;

        //保证在1-delta的概率 查询误差在eps以内

        reqquantile = quantile;
        //theta = 0.00006; //1000 / 2000_0000
        //theta = 0.0001;
        theta = 1000.00 / (double)running_length;
        std::cout << "theta: " << theta << "\n";

        epsilon = eps;       
        long long internal_counters = (int) ((double) (pow(theta, -1) * (double) pow(epsilon,-0.5)));
        long long samples_num = 4* ((int) ((double) (pow(theta, -2) * (double) pow(epsilon,-2)) * pow((double) internal_counters, -1)));
        double qsketch_eps = epsilon / (double) 2.0;
        
        LCU_type* squad = LCU_Init((double) 1/ (double) internal_counters, qsketch_eps, samples_num);
        double size_squad = LCU_Size(squad)/1024;
        std::cout<<"size: "<<size_squad<<" kb"<<std::endl;
        CorrectDetector<uint64_t>* correct_detector = new CorrectDetector<uint64_t>(); 
        std::map<uint64_t, uint32_t> id_map;
        std::map<uint64_t, uint64_t> last_time;
        std::set<uint64_t> id_set;

        clock_t begin,finish;
        clock_t total=0;
        double totaltime;

        std::vector<pair<uint32_t, double> > v; v.clear();

        for (int i = 0; i < running_length; ++i) 
        {
            uint32_t id = (int)(dataset[i].id.get_hash(MOD) & ((1ll << 30) - 1));

            if (i%2000000==0) printf("i = %d\n",i);
            if (id_set.find(id) == id_set.end()) 
            {
                id_set.insert(id);
                id_map[id] = 0;
            }
            id_map[id]++;

            if (id_map[id] > 1)
            {
                uint64_t latency = dataset[i].timestamp - last_time[id];
                correct_detector->insert(id, 0, latency);
                //if (latency == 0) std::cout<<i<<"\n";

                v.push_back(std::make_pair(id, (double)latency));
            }
            
            last_time[id] = dataset[i].timestamp;
        }

        //std::cout << v.size() << std::endl;

        begin = clock();
        for (int i = 0; i < v.size(); i++)
        {
            LCU_UpdateLatency(squad, v[i].first, v[i].second);
        }
        total = clock() - begin;


        totaltime=(double)(total)/CLOCKS_PER_SEC;
        std::cout << running_length << " " << v.size() << std::endl;
        //std::cout <<"insert time taken: "<<totaltime<<" seconds"<< std::endl;
        double throughput = double((int)v.size()) / totaltime;
        //std::cout <<"insert throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;
        
        double error = 0;
        int num = 0;
        double query_w = quantile;        
        total = 0;
        uint64_t total_id = 0;
        uint32_t zero_p = 0;
        for (auto i : id_map) 
        {
            total_id++;
            if (i.second < 1000)//1000可以改大，例如2000，不要出现1000 和 2000效果差别巨大的情况
                continue;
            num++;


            begin = clock();
            double est = LCU_QuantileQuery(squad, i.first, query_w); //是latency
            finish = clock();
            total += finish-begin;

            double predict = correct_detector->query_quantile(i.first, (uint32_t)est);
            error += abs( predict - query_w );

            //if (predict == 0) zero_p++;
            /*
            std::cout << log2(predict) << " " << log2(truth) << " " << "\n";
            // assert(0);
            
            printf("id=%llu, predict=%llu, truth=%llu, abs=%f, fabs=%f \n",i.first,predict,truth,
                                                                        abs( log2(predict) - log2(truth) ),
                                                                        fabs( log2(predict) - log2(truth) ));
            if (num==20) break;
            */
            
        }
        /*totaltime=(double)(total)/CLOCKS_PER_SEC;
        std::cout <<"query time taken: "<<totaltime<<" seconds"<< std::endl;
        throughput = double(num) / totaltime;
        std::cout <<"query throughput: "<<std::fixed<<std::setprecision(4)<< throughput <<" item/s"<< std::endl;

        std::cout << "Total id: " << total_id << "\n";
        std::cout << "Estimate: " << num << "\n";
        std::cout << "Average Log Error: " << log_error / num << "\n";
        std::cout << "query_w: " << query_w<< "\n";
        std::cout << "predice = 0: " << zero_p << "\n";*/

        return std::make_pair(throughput, error / num);

    }

private:
    std::string filename;
    LoadResult load_result;
    Tap_Tuple* dataset;
    uint64_t length;
};

#endif
