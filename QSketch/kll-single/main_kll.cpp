#include <bits/stdc++.h>
#include "Benchmark_kll.h"
#include <cmath>
#include <random>
#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

int main(int argc, char **argv)
{
    freopen("data.out", "w", stdout);

    std::string path_caida1 = "../CAIDA2016/1.dat";
    std::string path_caida2 = "../CAIDA2016/formatted02.dat";
    std::string path_caida3 = "../CAIDA2016/formatted03.dat";
    std::string path4 = "../1.dat";

    std::pair<double, double> res[33], cur;
    int cnt = 0;

    CAIDABenchmark benchmark(path_caida1);
    for (double w = 0.10; w < 0.95; w += 0.10)
    {
        cur = benchmark.Run(5000, 500, w);
        res[++cnt] = cur;
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

    cnt = 0;

    synthetic_Benchmark benchmark4(path4);
    for (double w = 0.1; w < 0.95; w += 0.10)
    {
		++cnt;
        cur = benchmark4.Run(5000, 500, w);
        res[cnt] = cur;
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
    
    cnt = 0;

    return 0;
    
}
