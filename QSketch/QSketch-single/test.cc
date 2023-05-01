#include <bits/stdc++.h>
#include "Benchmark.h"
using namespace std;

double w, tower, res, mn_f, mn;
int memory, d, threshold, m, mn_i, dd, M;

int main(int argc, const char** argv) {

    freopen("test.out", "w", stdout);
	string path = "../CAIDA2016/1.dat";

    CAIDA_Benchmark benchmark(path); // w, tower, threshold, memory, d, m
    
    /*printf("d:\n");
    for (memory = 5000; memory <= 15000; memory += 5000)
    {
        m = sqrt(memory) / 2; mn = 1e9;
        for (d = 1; d <= 5; d++)
        {
            res = benchmark.Run(w, 0.1, 40, memory, d, m).second;
            std::cout << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_i = d;
            }
        }
        std::cout << mn_i << std::endl;
    }

    printf("tower:\n");
    for (memory = 5000; memory <= 15000; memory += 5000)
    {
        m = sqrt(memory) / 2; mn = 1e9;
        for (tower = 0.05; tower < 0.31; tower += 0.05)
        {
            res = benchmark.Run(w, tower, 40, memory, 3, m).second;
            std::cout << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_f = tower;
            }
        }
        std::cout << mn_f << std::endl;
    }

    printf("threshold:\n");
    for (memory = 5000; memory <= 15000; memory += 5000)
    {
        m = sqrt(memory) / 2; mn = 1e9;
        
        for (threshold = 10; threshold <= 50; threshold += 10)
        {
            res = benchmark.Run(w, 0.1, threshold, memory, 3, m).second;
            std::cout << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_i = threshold;
            }
        }
        std::cout << mn_i << std::endl;
    }*/

    /*printf("M:\n");
    for (w = 0.1; w < 0.95; w += 0.1)
    {
        for (memory = 5000; memory <= 15000; memory += 5000)
        {
            mn = 1e9;
            for (M = 10; M <= 60; M ++)
            {
                res = benchmark.Run(w, 0.1, 40, memory, 3, M, 11).second;
                std::cout << res << "	";
                if (res < mn)
                {
                    mn = res;
                    mn_i = M;
                }
            }
            std::cout << w << " " << memory << " " << mn_i << std::endl;
        }
    }
    */

	return 0;
}