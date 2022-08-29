#include "algo.hpp"
#include "suffix_tree.hpp"

#include <string>
#include <functional>
#include <algorithm>
#include <vector>

using std::string;
using std::function;
using std::min;
using std::max;
using std::vector;
using std::sort;
using std::string_view;
using std::lower_bound;

LCFwM_result the_algorithm(string_view s1, string_view s2, int k, float) {
    const int n1 = s1.length();
    const int n2 = s2.length();

    // suffix tree data structures
    string s12 = string(s1) + "#" + string(s2) + "$";
    suffix_tree tree12(s12);
    st_lca lca12(tree12);

    // regular LCF
    int lcf_value = 0;

    function<int(suffix_tree::node_t*)> lcf = [&](suffix_tree::node_t* node) -> int {
        if (node->child) {
            int mask = 0x0;
            for (auto child = node->child; child; child = child->sibling)
                mask |= lcf(child);
            if (mask == 0x3)
                lcf_value = max(lcf_value, node->depth);
            return mask;
        }
        else
            return node->start < n1 ? 0x1 : 0x2;
    };

    lcf(tree12.root);
    if (lcf_value == 0)
        return { min(k, min(n1, n2)), 0, 0 };

    // algorithm propper
    struct nei_entry {
        int start;
    };

    int found_i1, found_i2;

    auto check_len = [&](int len) -> bool {
        vector<int> dels(k);
        for (int j = 0; j < k; j++)
            dels[j] = j;

        // comparator function
        auto nei_entry_comp = [&](const nei_entry& ne1, const nei_entry& ne2) {
            // lexicographical order
            for (int del_it = 0, str_it = 0; str_it < len; ) {
                int next_stop = del_it < k ? dels[del_it++] : len;
                if (str_it == next_stop) {
                    str_it++;
                    continue;
                }
                int lce = lca12.lca(ne1.start + str_it, ne2.start + str_it)->depth;
                if (lce < next_stop - str_it)
                    return s12[ne1.start + str_it + lce] < s12[ne2.start + str_it + lce];
                str_it = next_stop + 1;
            }

            return false;
        };

        // neighbourhood generation procedures
        auto next_comb = [&]() -> bool {
            int i = k - 1, last = len;
            while (i >= 0 && dels[i] == last - 1)
                last = dels[i--];
            if (i < 0)
                return false;
            dels[i]++;
            while (++i < k)
                dels[i] = dels[i - 1] + 1;
            return true;
        };

        // prepare neighbourhood
        vector<nei_entry> nei1;
        for (int i = 0; i < n1 - len + 1; i++)
            nei1.push_back({i});

        do {
            // sort the neighbourhood
            sort(nei1.begin(), nei1.end(), nei_entry_comp);

            // generate second neighbourhood and look for match
            for (int i = 0; i < n2 - len + 1; i++) {
                nei_entry entry = {n1 + 1 + i};
                auto it = lower_bound(nei1.begin(), nei1.end(), entry, nei_entry_comp);
                if (it != nei1.end() && !nei_entry_comp(entry, *it)){
                    found_i1 = it->start;
                    found_i2 = i;
                    return true;
                }
            }
        } while(next_comb());

        return false;
    };

    int len_at_least = min(max(lcf_value, k), min(n1, n2));
    int len_at_most = min((k + 1) * lcf_value + k, min(n1, n2));
    while (len_at_least < len_at_most) {
        int mid_len = (len_at_least + len_at_most + 1) / 2;
        if (check_len(mid_len))
            len_at_least = mid_len;
        else
            len_at_most = mid_len - 1;
    }

    check_len(len_at_least);
    return { len_at_least, found_i1, found_i2 };
}

