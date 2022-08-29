#include <random>
#include <iostream>
#include <cassert>

std::random_device rd;

int main(int argc, char** argv) {
    // parse command-line arguments
    // format: gen len Sigma_size
    assert(argc == 3);
    size_t len = std::stoll(argv[1]);
    int sigma = std::stoi(argv[2]);
    assert('a' + sigma - 1 <= 'z');
    std::uniform_int_distribution<char> distribution('a', 'a' + sigma - 1);

    for (size_t i = 0; i < len; i++)
        std::cout << distribution(rd);
    std::cout << std::endl;

    return 0;
}
