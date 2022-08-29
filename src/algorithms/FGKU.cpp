#include "algo.hpp"

#include <queue>
#include <algorithm>
#include <string_view>

using std::queue;
using std::string_view;
using std::min;

LCFwM_result the_algorithm(string_view s1, string_view s2, int k, float) {
    const int n1 = s1.length();
    const int n2 = s2.length();

    int l = 0, r1 = 0, r2 = 0;

    auto diagonal = [&](int i, int j) {
        queue<int> Q;
        int p = 0, s = 0;

        int limit = min(n1 - i, n2 - j);
        while (p < limit) {
            if (s1[i + p] != s2[j + p]) {
                Q.push(p);
                if ((int)Q.size() > k) {
                    s = Q.front() + 1;
                    Q.pop();
                }
            }
            p++;

            if (p - s > l) {
                l = p - s;
                r1 = i + s;
                r2 = j + s;
            }
        }
    };

    for (int j = 0; j < n2; j++)
        diagonal(0, j);
    for (int i = 1; i < n1; i++)
        diagonal(i, 0);

    return {l, r1, r2};
}
