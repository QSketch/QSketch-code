#include <bits/stdc++.h>
#include "QuantileSketch.h"
#include <sys/time.h>
using namespace std;

#define len 19666763

uint64_t id[22222222], val[22222222], a;
double f;
vector<pair<uint64_t, uint64_t> > ans;

int main(int argc, const char** argv) {

    freopen("flow.in", "r", stdin);

    vector<int> v;
    v.push_back(15); v.push_back(2); v.push_back(2);
    QuantileSketch<uint64_t, uint64_t> ql(100, 10, 0.05, 10, 5, 4.0, v);
    QuantileSketch<uint64_t, uint64_t> qr(100, 10, 0.05, 10, 5, 4.0, v);
    ql.set(0.01);
    qr.set(0.99);

    for (int i = 0; i < len; i++)
    {
        scanf("%llu%lf", &a, &f);
        id[i] = a;
        val[i] = (uint64_t)(f * 100.00);
    }

    struct timeval t_start, t_end;
    gettimeofday( &t_start, NULL );

    for (int i = 0; i < len; i++)
    {
        ql.insert(id[i], val[i]);
        qr.insert(id[i], val[i]);
    }
    
    

    for (int i = 0; i < len; i++)
    {
        if (val[i] >= ql.query(id[i], 0.01) && val[i] <= qr.query(id[i], 0.99))
        {
            ans.push_back(make_pair(id[i], val[i]));
        }
    }

    gettimeofday( &t_end, NULL );

    printf("%lld\n", 1000000ll * (t_end.tv_sec - t_start.tv_sec) + t_end.tv_usec - t_start.tv_usec);

	return 0;
}