#include <bits/stdc++.h>
#include "QuantileSketch.h"
#include <sys/time.h>
using namespace std;

#define len 19666763

uint64_t id[22222222], val[22222222], a;
int ln, idx;
double f;
vector<pair<uint64_t, uint64_t> > ans;
map<uint64_t, vector<uint64_t> > mp;
vector<uint64_t> v;

int main(int argc, const char** argv) {

    freopen("flow.in", "r", stdin);

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
        mp[id[i]].push_back(val[i]);
    }

    for (map<uint64_t, vector<uint64_t> >::iterator it = mp.begin(); it != mp.end(); it++)
    {
        sort(mp[it -> first].begin(), mp[it -> first].end());
    }

    for (int i = 0; i < len; i++)
    {
        idx = lower_bound(mp[id[i]].begin(), mp[id[i]].end(), val[i]) - mp[id[i]].begin();
        if (idx >= (int)(0.01 * mp[id[i]].size()) && idx <= (int)(0.99 * mp[id[i]].size()))
        {
            ans.push_back(make_pair(id[i], val[i]));
        }
    }
    
    gettimeofday( &t_end, NULL );

    printf("%lld\n", 1000000ll * (t_end.tv_sec - t_start.tv_sec) + t_end.tv_usec - t_start.tv_usec);

	return 0;
}