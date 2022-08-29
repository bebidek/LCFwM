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
using std::string_view;


LCFwM_result the_algorithm(string_view s1, string_view s2, int k, float) {
    int n1 = s1.length();
    int n2 = s2.length();

    // suffix tree data structures
    string s12 = string(s1) + "#" + string(s2) + "$";
    string s12r = string(s1.rbegin(), s1.rend()) + "#" + string(s2.rbegin(), s2.rend()) + "$";
    suffix_tree tree12(s12), tree12r(s12r);
    st_lca lca12(tree12), lca12r(tree12r);

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

    // algorithm proper
    auto forward_lce = [&](int i1, int i2) -> int {
        return lca12.lca(i1, n1 + 1 + i2)->depth;
    };
    auto backward_lce = [&](int j1, int j2) -> int {
        return lca12r.lca(n1 - j1 - 1, s12r.length() - j2 - 2)->depth;
    };

    auto next_mismatch_forward = [&](int i1, int i2, int last_pos) {
        if (i1+last_pos >= n1 || i2+last_pos >= n2)
            return last_pos;
        if (i1+last_pos+1 >= n1 || i2+last_pos+1 >= n2)
            return last_pos + 1;
        return last_pos + forward_lce(i1+last_pos+1, i2+last_pos+1) + 1;
    };
    auto next_mismatch_backward = [&](int i1, int i2, int last_pos) {
        if (i1+last_pos < 0 || i2+last_pos < 0)
            return last_pos;
        if (i1+last_pos-1 < 0 || i2+last_pos-1 < 0)
            return last_pos - 1;
        return last_pos - backward_lce(i1+last_pos-1, i2+last_pos-1) - 1;
    };


    int best_len = 0, best_i1 = 0, best_i2 = 0;

    auto diagonal = [&](int i, int j, int h) {
        vector<int>forward_mm(k+1),backward_mm(k+1);
        int limit = min(n1 - i, n2 - j);

        for (int p = 0; p < limit; p += h) {
            int last_f = p - 1, last_b = p;
            for (int it = 0; it <= k; it++) {
                last_f = next_mismatch_forward(i, j, last_f);
                forward_mm[it] = last_f;
                last_b = next_mismatch_backward(i, j, last_b);
                backward_mm[k - it] = last_b;
            }

            for (int it = 0; it <= k; it++) {
                int result = forward_mm[it] - backward_mm[it] - 1;
                if (result > best_len) {
                    best_len = result;
                    best_i1 = i + backward_mm[it] + 1;
                    best_i2 = j + backward_mm[it] + 1;
                }
            }
        }
    };


    int h = min(max(n1, n2), (k+1)*lcf_value+k);
    while(true) {
        for (int j = 0; j < n2; j++)
            diagonal(0, j, h);
        for (int i = 1; i < n1; i++)
            diagonal(i, 0, h);
        if (best_len >= h)
            break;
        h /= 2;
    }

    return {best_len, best_i1, best_i2};
}

