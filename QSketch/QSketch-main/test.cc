#include <bits/stdc++.h>
#include "Benchmark.h"
using namespace std;

double w, tower, res, mn_f, mn, lambda;
int memory, d, threshold, m, r, s, mn_i, dd, M, layer;

int main(int argc, const char** argv) {

    freopen("test.out", "w", stdout);
	string path = "../CAIDA2016/1.dat";

    CAIDA_Benchmark benchmark(path); // w, tower, threshold, memory, d, r, m, s, lambda

    w = 0.5; tower = 0.1; threshold = 40; d = 7; r = 15; s = 11; lambda = 4.00;

    //layer = 3; m = 2;
    
    std::vector<int> v;
    
    /*printf("d:\n");
    for (memory = 500; memory <= 1500; memory += 500)
    {
        mn = 1e9;
        for (int dd = 3; dd <= 13; dd += 2)
        {
            v.clear();
            v.push_back(r); v.push_back(2);v.push_back(2);v.push_back(2);
            res = benchmark.Run(w, tower, threshold, memory, s, dd, lambda, v).second;
            std::cout << std::fixed<<std::setprecision(4) << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_i = dd;
            }
        }
        std::cout << mn_i << std::endl;
    }*/

    /*printf("tower:\n");
    for (memory = 500; memory <= 1500; memory += 500)
    {
        mn = 1e9;
        for (double ttower = 0.05; ttower < 0.31; ttower += 0.05)
        {
            v.clear();
            v.push_back(r); v.push_back(2);v.push_back(2);v.push_back(2);
            res = benchmark.Run(w, ttower, threshold, memory, s, d, lambda, v).second;
            std::cout << std::fixed<<std::setprecision(4) << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_f = ttower;
            }
        }
        std::cout << mn_f << std::endl;
    }*/

    /*printf("threshold:\n");
    for (memory = 500; memory <= 1500; memory += 500)
    {
        mn = 1e9;
        for (int tthreshold = 10; tthreshold <= 50; tthreshold += 10)
        {
            v.clear();
            v.push_back(r); v.push_back(2);v.push_back(2);v.push_back(2);
            res = benchmark.Run(w, tower, tthreshold, memory, s, d, lambda, v).second;
            std::cout << std::fixed<<std::setprecision(4) << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_i = tthreshold;
            }
        }
        std::cout << mn_i << std::endl;
    }*/

    /*printf("lambda:\n");
    for (memory = 500; memory <= 1500; memory += 500)
    {
        mn = 1e9;
        for (double llambda = 1.0; llambda < 6.5; llambda += 1.0)
        {
            v.clear();
            v.push_back(r); v.push_back(2);v.push_back(2);v.push_back(2);
            res = benchmark.Run(w, tower, threshold, memory, s, d, llambda, v).second;
            std::cout << std::fixed<<std::setprecision(4) << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_f = llambda;
            }
        }
        std::cout << mn_f << std::endl;
    }*/

    printf("r:\n");
    for (memory = 500; memory <= 1500; memory += 500)
    {
        mn = 1e9;
        for (int rr = 10; rr <= 40; rr += 6)
        {
            v.clear();
            v.push_back(rr); v.push_back(2);v.push_back(2);v.push_back(2);
            res = benchmark.Run(w, tower, threshold, memory, s, d, lambda, v).second;
            std::cout << std::fixed<<std::setprecision(4) << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_i = rr;
            }
        }
        std::cout << mn_i << std::endl;
    }


    printf("s:\n");
    for (memory = 500; memory <= 1500; memory += 500)
    {
        mn = 1e9;
        for (int ss = 4; ss <= 14; ss += 2)
        {
            v.clear();
            v.push_back(r); v.push_back(2);v.push_back(2);v.push_back(2);
            res = benchmark.Run(w, tower, threshold, memory, ss, d, lambda, v).second;
            std::cout << std::fixed<<std::setprecision(4) << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_i = ss;
            }
        }
        std::cout << mn_i << std::endl;
    }

   /* printf("layers:\n");
    for (memory = 500; memory <= 500; memory += 2000)
    {
        mn = 1e9;
        for (int layer = 1; layer <= 5; layer ++)
        {
            v.clear();
            v.push_back(r);
            for (int i = 0; i < layer; i++) v.push_back(2);
            res = benchmark.Run(w, tower, threshold, memory, s, d, lambda, v).second;
            std::cout << std::fixed<<std::setprecision(4) << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_i = layer;
            }
        }
        std::cout << mn_i << std::endl;
    }

    printf("m:\n");
    for (memory = 500; memory <= 500; memory += 2000)
    {
        mn = 1e9;
        for (int m = 2; m <= 8; m++)
        {
            v.clear();
            v.push_back(r); v.push_back(m); v.push_back(m); v.push_back(m);
            res = benchmark.Run(w, tower, threshold, memory, s, d, lambda, v).second;
            std::cout << std::fixed<<std::setprecision(4) << res << "	";
            if (res < mn)
            {
                mn = res;
                mn_i = m;
            }
        }
        std::cout << mn_i << std::endl;
    }*/

	return 0;
}