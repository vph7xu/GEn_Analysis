#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <regex>
#include <algorithm>
#include <limits>

// SortFOM.C — ROOT macro to read lines from a text file, sort by FOM_1 descending,
// and write results to an output file.
// Usage in ROOT:
//   root> .L SortFOM.C+
//   root> SortFOM("input.txt", "sorted_output.txt");

void SortFOM(const char* inputFile, const char* outputFile) {
    std::ifstream infile(inputFile);
    if (!infile) {
        std::cerr << "Error: cannot open input file '" << inputFile << "'" << std::endl;
        return;
    }

    std::ofstream outfile(outputFile);
    if (!outfile) {
        std::cerr << "Error: cannot open output file '" << outputFile << "'" << std::endl;
        return;
    }

    std::vector<std::pair<double, std::string>> entries;
    std::string line;
    std::regex re(R"(FOM_1\s*:\s*([\-+]?nan|[\-+]?\d*\.?\d+(?:[eE][\-+]?\d+)?))");
    std::smatch m;

    while (std::getline(infile, line)) {
        if (std::regex_search(line, m, re)) {
            std::string sval = m[1];
            double dval;
            if (sval == "nan" || sval == "+nan" || sval == "-nan") {
                dval = -std::numeric_limits<double>::infinity();
            } else {
                try {
                    dval = std::stod(sval);
                } catch (...) {
                    dval = -std::numeric_limits<double>::infinity();
                }
            }
            entries.emplace_back(dval, line);
        }
    }

    std::sort(entries.begin(), entries.end(),
              [](auto &a, auto &b) { return a.first > b.first; });

    for (auto &p : entries) {
        outfile << p.second << std::endl;
    }

    outfile.close();
    std::cout << "Saved sorted list to '" << outputFile << "'" << std::endl;
}

