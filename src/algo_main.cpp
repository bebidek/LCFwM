#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>
#include <cassert>

#include "algo.hpp"

using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::getline;

int main(int argc, char** argv) {
    // parse command-line arguments
    // format: LCSwM pathA pathB k epsilon
    assert(argc == 5);
    string path_s1 = argv[1];
    string path_s2 = argv[2];
    int param_k = std::stoi(argv[3]);
    float param_eps = std::stof(argv[4]);

    // load strings from files
    string S1, S2;
    getline(ifstream(path_s1), S1);
    getline(ifstream(path_s2), S2);

    // start time measurement
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    auto time_before = std::chrono::high_resolution_clock::now();

    // run the algorithm
    auto result = the_algorithm(S1, S2, param_k, param_eps);

    // finish time measurement, print results
    auto time_after = std::chrono::high_resolution_clock::now();
    cout << std::chrono::duration_cast<std::chrono::microseconds>(time_after - time_before).count() << endl;
    cout << result.len << " " << result.s1_it << " " << result.s2_it << endl;

    // verify number of mismatches
    int mm = 0;
    for (ssize_t i = 0; i < result.len; i++)
        if (result.s1_it + i >= (int)S1.length() || result.s2_it + i >= (int)S2.length()) {
            cout << "FAIL" << endl;
            return 1;
        }
        else if (S1[result.s1_it + i] != S2[result.s2_it + i])
            mm++;
    cout << mm << endl;

    return 0;
}
