#include <bits/stdc++.h>
#include "Benchmark.h"
using namespace std;

int threshold[33], mm[33], d[33], y;

double tower[33], x, z;

int main(int argc, const char** argv) {
	freopen("data.out", "w", stdout);
	string path1 = "../CAIDA2016/1.dat";
	string path2 = "../zipf_2022/zipf_2.0.dat";
	string path3 = "../tnsm-19-dataset/dataset/webget-all.csv";
    string path4 = "../1.dat";

	std::pair<std::pair<double, double>, double> res[33], cur;
    int cnt = 0;

    CAIDA_Benchmark benchmark(path1);
    for (double w = 0.10; w < 0.95; w += 0.10)
    {
		++cnt;
        cur = benchmark.Run(w, 10);
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
    for (double w = 0.10; w < 0.95; w += 0.10)
    {
		++cnt;
        cur = benchmark4.Run(w, 10);
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