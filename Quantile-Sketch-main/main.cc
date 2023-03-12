#include <bits/stdc++.h>
#include "Benchmark.h"
using namespace std;

int mm;
double f[11][11];

int main(int argc, const char** argv) {
	//freopen("data.out", "w", stdout);
	string path1 = "../CAIDA2016/1.dat";
	string path2 = "../zipf_2022/zipf_2.0.dat";
	string path3 = "../tnsm-19-dataset/dataset/webget-all.csv";

    /*freopen("data.out", "r", stdin);
    freopen("data2.out", "w", stdout);
    for (int i = 1; i <= 6; i++)
    {
        for (int j = 1; j <= 9; j++)
        {
            for (int k = 1; k <= 5; k++)
            scanf("%lf", &f[j][k]);
        }
        for (int j = 1; j <= 9; j++)
        {
            for (int k = 1; k <= 3; k++)
            {
                std::cout << std::fixed<<std::setprecision(4) << f[j][k] << "	";
            }
        }
        std::cout << std::endl;
        for (int j = 1; j <= 9; j++)
        {
            for (int k = 4; k <= 5; k++)
            {
                std::cout << std::fixed<<std::setprecision(4) << f[j][k] << "	";
            }
        }
        std::cout << std::endl;
    }

    return 0;*/

	std::pair<double, double> res[33], cur;
    int cnt = 0;

    CAIDA_Benchmark benchmark(path1);
    for (double w = 0.10; w < 0.95; w += 0.10)
    {
        for (int memory = 5000; memory <= 25000; memory += 5000)
        {
			++cnt; mm = sqrt(memory) / 2;
            cur = benchmark.Run(w, 0.1, 40, memory, 3, mm + 2 * sqrt(mm));
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

    cnt = 0;
    zipf_Benchmark benchmark2(path2);

    for (double w = 0.10; w < 0.95; w += 0.10)
    {
        for (int memory = 5000; memory <= 25000; memory += 5000)
        {
			++cnt; mm = sqrt(memory) / 2;
            cur = benchmark2.Run(w, 0.1, 40, memory, 3, mm);
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

    cnt = 0;

    tnsm_Benchmark benchmark3(path3);

    for (double w = 0.10; w < 0.95; w += 0.10)
    {
        for (int memory = 5000; memory <= 25000; memory += 5000)
        {
            ++cnt;
            cur = benchmark3.Run(w, 0.1, 40, memory, 3, 1000);
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