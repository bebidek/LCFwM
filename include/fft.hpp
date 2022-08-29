#ifndef FFT_HPP
#define FFT_HPP

#include <vector>
#include <algorithm>

using std::vector;
using std::swap;


struct FFT {
    static const int MOD = 2013265921; //15*2^27+1
    static const int ROOT = 440564289; //GEN^15

    static int fastpow(int b, int e) {
        long long result = 1, pwr = b;
        while (e) {
            if (e & 1)
                result = (result * pwr) % MOD;
            pwr = (pwr * pwr) % MOD;
            e >>= 1;
        }
        return (int)result;
    }


    int size, inv_size;
    vector<int> permutation;
    vector<vector<int>> omegas;


    FFT(int decl_size) {
        // power of 2
        int real_size = 1, real_exp = 0;
        while (real_size < decl_size) {
            real_size *= 2;
            real_exp++;
        }
        size = real_size;
        inv_size = fastpow(size, MOD - 2);

        // prepare permutation
        permutation.resize(size);
        for (int i = 0; i < size; i++) {
            int inv = 0, ii = i;
            for (int j = 0; j < real_exp; j++) {
                inv = (inv << 1) | (ii & 1);
                ii >>= 1;
            }
            permutation[i] = inv;
        }

        // calcualte omegas
        long long base_omega = fastpow(ROOT, 1 << (27 - real_exp));
        omegas.push_back(vector<int>(size + 1));
        omegas[0].front() = omegas[0].back() = 1;
        for (int i = 1; i < size; i++)
            omegas[0][i] = (omegas[0][i-1] * base_omega) % MOD;
        while (omegas.back().size() > 3) {
            vector<int>omega(omegas.back().size() / 2 + 1);
            for (int i = 0; i < (int)omega.size(); i++)
                omega[i] = omegas.back()[2*i];
            omegas.push_back(omega);
        }
        reverse(omegas.begin(), omegas.end());
    }


    void fft(vector<int> &V, bool inv = false) {
        // rearange input vector
        for (int i = 0; i < size; i++)
            if (permutation[i] < i)
                swap(V[i], V[permutation[i]]);

        // butterfly operation
        auto butterfly = [&](int i, int j, int omega) {
            long long omega_j = ((long long)omega * V[j]) % MOD;
            long long new_i = V[i] + omega_j;
            long long new_j = V[i] - omega_j;
            if (new_i >= MOD)
                new_i -= MOD;
            if (new_j < 0)
                new_j += MOD;
            V[i] = new_i;
            V[j] = new_j;
        };

        // iterations
        if (inv) {
            for (int level = 0, delta = 1; delta < size; level++, delta *= 2)
                for (int i = 0; i < size; )
                    if (i & delta) {
                        butterfly(i - delta, i, omegas[level][2*delta - (i & (delta - 1))]);
                        i++;
                    }
                    else
                        i += delta;
        }
        else {
            for (int level = 0, delta = 1; delta < size; level++, delta *= 2)
                for (int i = 0; i < size; )
                    if (i & delta) {
                        butterfly(i - delta, i, omegas[level][i & (delta - 1)]);
                        i++;
                    }
                    else
                        i += delta;
        }
    }

    void mult(vector<int> &A, vector<int>&B) {
        for (int i = 0; i < size; i++)
            A[i] = ((long long)A[i] * B[i]) % MOD;
    }

    void ifft(vector<int> &V) {
        fft(V, true);
        for (int i = 0; i < size; i++)
            V[i] = ((long long)V[i] * inv_size) % MOD;
    }
};

#endif
