#include <bits/stdc++.h>
#include "Benchmark.h"
using namespace std;

int threshold[33], mm[33], d[33], y;

double tower[33], x, z;

int main(int argc, const char** argv) {
	//freopen("data.out", "w", stdout);
	string path1 = "../CAIDA2016/1.dat";
	string path2 = "../zipf_2022/zipf_2.0.dat";
	string path3 = "../tnsm-19-dataset/dataset/webget-all.csv";

	double res[33], tput[33];
    int cnt1 = 0, cnt2 = 0;

    CAIDA_Benchmark benchmark(path1);
    for (int memory = 5000; memory <= 25000; memory += 5000)
    {
        benchmark.Run(0.1, 40, memory, 2, sqrt(memory) / 2, res, tput, cnt1, cnt2);
    }


    for (int i = 1; i <= cnt2; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << tput[i] << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt1; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i] << "	";
    }

    std::cout << std::endl;

    cnt1 = cnt2 = 0;

    zipf_Benchmark benchmark2(path2);
    for (int memory = 5000; memory <= 25000; memory += 5000)
    {
        benchmark2.Run(0.1, 40, memory, 2, sqrt(memory) * 5, res, tput, cnt1, cnt2);
    }


    for (int i = 1; i <= cnt2; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << tput[i] << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt1; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i] << "	";
    }

    std::cout << std::endl;

    cnt1 = cnt2 = 0;

    tnsm_Benchmark benchmark3(path3);
    for (int memory = 5000; memory <= 25000; memory += 5000)
    {
        benchmark3.Run(0.1, 40, memory, 2, 1000, res, tput, cnt1, cnt2);
    }


    for (int i = 1; i <= cnt2; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << tput[i] << "	";
    }

    std::cout << std::endl;

    for (int i = 1; i <= cnt1; i++) 
    {
        std::cout << std::fixed<<std::setprecision(4) << res[i] << "	";
    }

    std::cout << std::endl;

	return 0;
}