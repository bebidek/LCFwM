#ifndef ALGO_HPP
#define ALGO_HPP

#include <string_view>
#include <cstdlib>

struct LCFwM_result {
    int len, s1_it, s2_it;
};

LCFwM_result the_algorithm(std::string_view s1, std::string_view s2, int k, float eps);

#endif
