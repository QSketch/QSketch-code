#include <bits/stdc++.h>
#include "Benchmark.h"
using namespace std;

int mm;
double f[11][11];
std::vector<int> v;

int main(int argc, const char** argv) {
	freopen("data.out", "w", stdout);
	string path1 = "../CAIDA2016/1.dat";
	string path2 = "../zipf_2022/zipf_1.0.dat";
	string path3 = "../tnsm-19-dataset/dataset/webget-all.csv";
    string path4 = "../1.dat";
    
    std::pair<double, double> res[333], cur;
    int cnt = 0;

    v.push_back(230);
    v.push_back(2);
    v.push_back(2);
    CAIDA_Benchmark benchmark(path1);
    for (double w = 0.10; w < 0.95; w += 0.10)
    {
        for (int memory = 10; memory <= 10; memory += 10)
        {
			++cnt;
            cur = benchmark.Run(w, memory, v);
            res[cnt] = cur;
        }
    }

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].first << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].second << "	";
    }

    std::cout << std::endl;

    synthetic_Benchmark benchmark4(path4);

    v.clear();
    v.push_back(920);
    v.push_back(2);
    v.push_back(2);
    v.push_back(2);
    v.push_back(2);
    v.push_back(2);
    v.push_back(2);

    cnt = 0;

    for (double w = 0.10; w < 0.95; w += 0.10)
    {
        for (int memory = 10; memory <= 10; memory += 10)
        {
            ++cnt;
            cur = benchmark4.Run(w, memory, v);
            res[cnt] = cur;
        }
    }

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].first << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].second << "	";
    }

    std::cout << std::endl;

	return 0;
}