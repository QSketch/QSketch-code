#include <bits/stdc++.h>
#include "Benchmark_gk.h"
#include <cmath>
#include <random>
#include <ctime>
#include <chrono>
using namespace std;

int main(int argc, const char** argv) {
	freopen("data.out", "w", stdout);
	string path1 = "../CAIDA2016/1.dat";
	string path2 = "../CAIDA2016/formatted02.dat";
	string path3 = "../CAIDA2016/formatted03.dat";
    string path4 = "../1.dat";

	std::pair<std::pair<double, double>, double> res[33], cur;
    int cnt = 0;

    CAIDABenchmark benchmark(path1);
    for (double w = 0.1; w < 0.95; w += 0.10)
    {
		++cnt;
        cur = benchmark.Run(5000, 0.004, w);
        res[cnt] = cur;
    }

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].first.first << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].first.second << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].second << "	";
    }

    std::cout << std::endl;
    
    cnt = 0;

    synthetic_Benchmark benchmark4(path4);
    for (double w = 0.1; w < 0.95; w += 0.10)
    {
		++cnt;
        cur = benchmark4.Run(5000, 0.004, w);
        res[cnt] = cur;
    }

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].first.first << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].first.second << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i].second << "	";
    }

    std::cout << std::endl;
    
    cnt = 0;

    return 0;
}
