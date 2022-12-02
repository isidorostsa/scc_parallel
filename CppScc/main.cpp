#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <array>
#include <string>
#include <unordered_set>
#include <chrono>

#include "colorSCC.hpp"

int main(int argc, char** argv) {
    std::string filename(argc > 1 ? argv[1] : "../matrices/uk-2002/uk-2002.mtx");
/*
    if(argc<2){
        std::cout << "Assumed " << filename <<  " as input" << std::endl;
    }

    std::cout << "Reading file '" << filename << "'\n";

    if(argc>2){
        times = std::stoi(argv[2]);
    }

    if(argc>3){
        DEBUG = std::stoi(argv[3]) == 1;
    }

    if(argc > 4){
        TOO_BIG = std::atoi(argv[4]) == 1;
    }
*/

    size_t times = 1;
    bool DEBUG = false;
    bool TOO_BIG = true;
    // load file into csr and time it

    auto start_load = std::chrono::high_resolution_clock::now();
    Sparse_matrix csr = loadFileToCSC(filename);
    Sparse_matrix csc;
    csr_tocsc(csc, csr);
    auto end_load = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_seconds = end_load-start_load;
    std::cout << "Loading time: " << elapsed_seconds.count() << "s\n";

/*
    Coo_matrix coo = loadFile(filename);
    Sparse_matrix csc;
    Sparse_matrix csr;

    std::cout << "Loaded matrix" << std::endl;
    if(TOO_BIG){
        std::cout << "Too big matrix, will go coo -> csr, csr -> csc instead of coo -> csr, csc" << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        coo_tocsr(coo, csr);
        coo.Ai = std::vector<size_t>();
        coo.Aj = std::vector<size_t>();
        csr_tocsc(csr, csc);
        auto end = std::chrono::high_resolution_clock::now();

        std::cout << "Conversion took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;
    } else {
        DEB("Starting coo -> csr, csc");
        auto now = std::chrono::high_resolution_clock::now();
        # pragma omp parallel sections
        {
            # pragma omp section
            {
                coo_tocsr(coo, csr);
            }
            # pragma omp section
            {
                coo_tocsc(coo, csc);
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - now;

        std::cout << "coo -> csr, csc took " << elapsed.count() << "s" << std::endl;
    }

*/


    std::cout << "Running " << times << " times\n";

    std::vector<size_t> SCC_id;
    auto start = std::chrono::high_resolution_clock::now();
    for(size_t i = 0; i < times; i++) {
        SCC_id = colorSCC_no_conversion(&csc, csr, DEBUG);
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Time w/o conversion: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()/times << "ms" << std::endl;
    std::cout << "SCC count: " << *std::max_element(SCC_id.begin(), SCC_id.end()) << std::endl;

    std::cout << "REAL SCC count: " << std::unordered_set(SCC_id.begin(), SCC_id.end()).size() << std::endl;

    return 0;
}