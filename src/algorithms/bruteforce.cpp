#include "algo.hpp"

LCFwM_result the_algorithm(std::string_view s1, std::string_view s2, int k, float) {
    int best_len = 0, best_i = 0, best_j = 0;

    for (int i = 0; i < (int)s1.length(); i++)
        for (int j = 0; j < (int)s2.length(); j++) {
            int len = 0, mm_left = k;

            while (i + len < (int)s1.length() && j + len < (int)s2.length()) {
                if (s1[i + len] != s2[j + len]) {
                    if (mm_left)
                        mm_left--;
                    else
                        break;
                }
                len++;
            }

            if (len > best_len) {
                best_len = len;
                best_i = i;
                best_j = j;
            }
        }

    return {best_len, best_i, best_j};
}
