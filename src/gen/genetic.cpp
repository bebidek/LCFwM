#include <random>
#include <iostream>
#include <cassert>
#include <fstream>

std::random_device rd;

char translate(char x) {
    switch (x) {
        case 'a':
        case 'A': return 'a';
        case 'c':
        case 'C': return 'b';
        case 'g':
        case 'G': return 'c';
        case 't':
        case 'T': return 'd';
        default: return 0;
    }
}

int main(int argc, char** argv) {
    // parse command-line arguments
    // format: gen len data_file
    assert(argc == 3);
    size_t len = std::stoll(argv[1]);
    std::string path = argv[2];
    
    std::ifstream file(path);
    assert(file.good());
    std::string sequence;
    while (!file.eof()) {
        std::string line;
        file >> line;
        if (line.empty() || line[0]=='>')
            continue;
        for (char c : line) {
            char ct = translate(c);
            if (ct)
                sequence.push_back(translate(c));
        }
    }

    int seq_len = sequence.length();
    int start = rd() % (seq_len - len + 1);
    std::cout << sequence.substr(start, len) << std::endl;

    return 0;
}

