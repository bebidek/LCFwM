#include <random>
#include <iostream>
#include <cassert>
#include <fstream>

std::random_device rd;

int main(int argc, char** argv) {
    // parse command-line arguments
    // format: gen len Sigma_size file1 file2
    assert(argc == 5);
    size_t len = std::stoll(argv[1]);
    int sigma = std::stoi(argv[2]);
    std::string path1 = argv[3];
    std::string path2 = argv[4];
    assert('a' + sigma - 1 <= 'z');
    std::uniform_int_distribution<char> distribution_c('a', 'a' + sigma - 1);
    std::uniform_int_distribution<int> distribution_p(0, len-1);

    std::string s1;
    for (size_t i = 0; i < len; i++)
        s1.push_back(distribution_c(rd));

    std::string s2(s1);
    for (size_t i = 0; i < size_t(0.1*len); i++)
        s2[distribution_p(rd)] = distribution_c(rd);

    std::ofstream of1(path1);
    assert(of1.good());
    of1 << s1 << std::endl;

    std::ofstream of2(path2);
    assert(of2.good());
    of2 << s2 << std::endl;

    return 0;
}
