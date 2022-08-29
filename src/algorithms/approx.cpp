#include "algo.hpp"
#include "suffix_tree.hpp"
#include "fft.hpp"

#include <optional>
#include <utility>
#include <string>
#include <stack>
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

using std::optional;
using std::pair;
using std::string;
using std::stack;
using std::min;
using std::max;
using std::swap;
using std::random_device;
using std::mt19937_64;
using std::uniform_int_distribution;
using std::vector;
using std::function;

const int Q = 2013265921; //15*2^27+1
const int R = 440564289;

random_device rd;
mt19937_64 mt(rd());

template<typename T>
struct reservoir_sampler {
    long long stored_samples = 0;
    T the_sample;

    void add(long long num, function<T(long long)> sample_getter) {
        stored_samples += num;
        long long next_sample = uniform_int_distribution<>(0, stored_samples - 1)(mt);
        if (next_sample < num)
            the_sample = sample_getter(next_sample);
    }
};

LCFwM_result the_algorithm(std::string_view s1, std::string_view s2, int k, float eps) {
    int n1 = s1.length();
    int n2 = s2.length();
    int n = min(n1, n2);
    int N = max(n1, n2);


    // FFT preparation
    FFT fft(N + n - 1);
    vector<int> V1(fft.size), V2(fft.size);
    for (int i = 0; i < n1; i++)
        V1[i] = s1[i];
    for (int i = 0; i < n2; i++)
        V2[i] = s2[i];
    fft.fft(V1);
    fft.fft(V2);


    // Suffix tree preparation
    string s12 = string(s1) + "#" + string(s2) + "$";
    suffix_tree tree12(s12);
    st_lca lca12(tree12);


    // Hamming distance computation
    auto distance_leq = [&](int i1, int i2, int len, int thr) -> bool {
        while (len > 0) {
            if (s1[i1] != s2[i2]) {
                i1++;
                i2++;
                len--;
                if (--thr < 0)
                    return false;
            }
            else {
                int lce = lca12.lca(i1, n1 + 1 + i2)->depth;
                i1 += lce;
                i2 += lce;
                len -= lce;
            }
        }

        return true;
    };


    // projection generation
    auto randomize_projection = [&](vector<int>&Pool, int len, int trials) {
        int size = 0;
        while (trials--) {
            int pos = uniform_int_distribution<>(0, len - 1)(mt);
            if (pos < size)
                continue;
            swap(Pool[size], Pool[pos]);
            size++;
        }
        return size;
    };


    // substring fingerprint structure
    struct Fingerprint {
        int pos;
        int value;
        
        bool operator<(const Fingerprint& other) const {
            return value < other.value;
        }
    };


    // algorithm for decision version
    auto decision = [&](int l) -> optional<pair<int,int>> {
        // special trivial cases
        if (l <= (1.0+eps)*k)
            return {{0, 0}};
        
        // compute constants
        float p2 = 1.0 - (1.0 + eps)*k/l;
        int L = ceil(pow(n, 1.0 / (1.0 + eps)) * 0.0625);
        int m = log(1.0 / n) / log(p2);
        int fft_thr = ceil(2.0 * log2(n));
        int relaxed_k = ceil((1.0+eps)*k);

        // prepate to projection generation
        vector<int> Pool;
        for (int i = 0; i < l; i++)
            Pool.push_back(i);

        // prepare structures for collisions
        reservoir_sampler<pair<int, int>> sampler;
        long long subset_checks_left = 4 * (long long)n * L;

        vector<long long>R_pwrs(m);
        R_pwrs[0] = 1;
        for (int i = 1; i < m; i++)
            R_pwrs[i] = (R_pwrs[i-1] * R) % Q;

        // find and process collisions
        vector<int>U(fft.size);
        vector<int>&X = U, Y(fft.size);
        vector<Fingerprint> FX(n1 - l + 1), FY(n2 - l + 1);

        while (L--) {
            int proj_size = randomize_projection(Pool, l, m);
            vector<int>&proj = Pool;

            if (proj_size < fft_thr) {
                // generate fingerprints manually
                for (int i = 0; i <= n1 - l; i++) {
                    FX[i] = {i, 0};
                    for (int j = 0; j < proj_size; j++)
                        FX[i].value = (FX[i].value + R_pwrs[j] * (s1[i + proj[j]])) % Q;
                }
                sort(FX.begin(), FX.end());
                for (int i = 0; i <= n2 - l; i++) {
                    FY[i] = {i, 0};
                    for (int j = 0; j < proj_size; j++)
                        FY[i].value = (FY[i].value + R_pwrs[j] * (s2[i + proj[j]])) % Q;
                }
                sort(FY.begin(), FY.end());
            }
            else {
                // get fingerprints
                fill(U.begin(), U.end(), 0);
                for (int i = 0; i < proj_size; i++)
                    U[l - 1 - proj[i]] = R_pwrs[i];
                fft.fft(U);
                for (int i = 0; i < fft.size; i++) {
                    Y[i] = ((long long)V2[i] * U[i]) % Q;
                    X[i] = ((long long)V1[i] * U[i]) % Q;
                }
                fft.ifft(X);
                fft.ifft(Y);

                // extract fingerprints
                for (int i = 0; i <= n1 - l; i++)
                    FX[i] = {i, X[i + l - 1]};
                sort(FX.begin(), FX.end());
                for (int i = 0; i <= n2 - l; i++)
                    FY[i] = {i, Y[i + l - 1]};
                sort(FY.begin(), FY.end());
            }

            // find collisions
            for (int i1 = 0, i2 = 0; i1 < (int)FX.size() && i2 < (int)FY.size(); )
                if (FX[i1].value < FY[i2].value)
                    i1++;
                else if (FX[i1].value > FY[i2].value)
                    i2++;
                else {
                    int x_colls = 0, y_colls = 0;
                    while (i1 + x_colls < (int)FX.size() && FX[i1 + x_colls].value == FX[i1].value)
                        x_colls++;
                    while (i2 + y_colls < (int)FY.size() && FY[i2 + y_colls].value == FY[i2].value)
                        y_colls++;

                    sampler.add((long long)x_colls * y_colls, [&](long long i){
                        return pair<long long,long long>{FX[i % x_colls].pos, FY[i / x_colls].pos};
                    });

                    if (subset_checks_left)
                        for (int xi = 0; xi < x_colls; xi++)
                            for (int yi = 0; yi < y_colls; yi++)
                                if (subset_checks_left) {
                                    subset_checks_left--;
                                    if (distance_leq(FX[i1 + xi].pos, FY[i2 + yi].pos, l, relaxed_k))
                                        return {{FX[i1 + xi].pos, FY[i2 + yi].pos}};
                                }

                    i1 += x_colls;
                    i2 += y_colls;
                }
        }

        // check the one collision
        if (distance_leq(sampler.the_sample.first, sampler.the_sample.second, l, relaxed_k))
            return {{sampler.the_sample.first, sampler.the_sample.second}};

        return {};
    };


    // 20 questions game
    int questions = 2 * ceil(log(n));
    stack<pair<int,int>> S;
    optional<pair<int,int>> response;
    int best_val = 0, best_i1 = 0, best_i2 = 0;

    while (questions--) {
        if (S.empty())
            S.push({0, n});
        auto [l,r] = S.top();
        int mid = (l + r + 1) / 2;

        response = decision(mid);
        if (response) {
            if (mid > best_val) {
                best_val = mid;
                best_i1 = response->first;
                best_i2 = response->second;
            }
            if (decision(r + 1))
                S.pop();
            else
                S.push({mid, r});
        }
        else {
            if (decision(l))
                S.push({l, mid-1});
            else
                S.pop();
        }
    }

    return {best_val, best_i1, best_i2};
}
