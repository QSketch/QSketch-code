#include <bits/stdc++.h>
#include "Benchmark_squad.h"
using namespace std;

int main(int argc, char** argv) {
    freopen("data.out", "w", stdout);
    string path1 = "../../CAIDA2016/1.dat";
	string path2 = "../../zipf_2022/zipf_2.0.dat";
	string path3 = "../../Seattle/SeattleData_";
    string path4 = "../../tap.dat";

    //SeattleBenchmark benchmark(path_seattle);
    //WebgetBenchmark benchmark(path_webget);
    //CAIDABenchmark benchmark(path_caida);
    //benchmark.Run(atoi(argv[1]),stod(argv[2])); //memory单位是kb

    std::pair<double, double> res[333], cur;
    int cnt = 0;

    Tap_Benchmark benchmark(path4);
    double w = 0.9;
	++cnt;
    cur = benchmark.Run(1500, 0.9, w, 1000.00 / 17790002.00);
    res[cnt] = cur;

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