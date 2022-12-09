#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <array>
#include <string>
#include <unordered_set>
#include <chrono>
#include <filesystem>

#include "sparse_util.hpp"
#include "colorSCC.hpp"

#define DEFAULT_PREFIX "../../matrices/"

void testFile(std::string filename, size_t times, bool DEBUG, bool TOO_BIG, size_t NUM_THREADS) {
    if(!std::filesystem::exists(filename)) {
        std::cout << "File not found: " << filename << std::endl;
        return;
    }

    if(times == 0) {
        std::cout << "Invalid number of times to run" << std::endl;
        return;
    }

    DEB("Loading file into CSC")
    auto start_load_csc = std::chrono::high_resolution_clock::now();
    Sparse_matrix csc = loadFileToCSC(filename);
    auto end_load_csc = std::chrono::high_resolution_clock::now();
    DEB("Loaded file into CSC, toook " << std::chrono::duration_cast<std::chrono::milliseconds>(end_load_csc - start_load_csc).count() << "ms")

    Sparse_matrix csr = Sparse_matrix();
    if (TOO_BIG) {
        DEB("File is too big, skipping CSR")
    } else {
        DEB("Making CSR Version")
        auto start_load_csr = std::chrono::high_resolution_clock::now();
        csc_tocsr(csc, csr);
        auto end_load_csr = std::chrono::high_resolution_clock::now();
        DEB("Making CSR, toook " << std::chrono::duration_cast<std::chrono::milliseconds>(end_load_csr - start_load_csr).count() << "ms")
    }

    DEB("Running " << times << " times")

    std::vector<size_t> SCC_id;
    std::vector<int64_t> times_in_us;

    std::string dataset_name = filename.substr(filename.find_last_of("/") + 1);
    dataset_name = dataset_name.substr(0, dataset_name.find_last_of("."));

    for(size_t i = 0; i < times; i++) {
        DEB("Starting run " << i)
        // csr may a null object, but that's fine, the next variable communicates that
        auto start = std::chrono::high_resolution_clock::now();
        if(!TOO_BIG) {
            SCC_id = colorSCC_no_conversion(csc, csr, true, DEBUG, NUM_THREADS);
        } else {
            SCC_id = colorSCC_no_conversion(csc, csr, false, DEBUG, NUM_THREADS);
        }
        auto end = std::chrono::high_resolution_clock::now();

        auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        times_in_us.push_back(time);
        std::cout << "DATASET: " << dataset_name << "\tTIME: " << time << "us" << std::endl;
        DEB("Finished run " << i)
    }
    DEB("Finished all runs")

    /*
    // print for spreadsheets
    for(size_t i = 0; i < times_in_us.size(); i++) {
        if (i == times - 1) {
            std::cout << times_in_us[i];
        } else {
            std::cout << times_in_us[i] << "\t";
        }
    }
    std::cout << std::endl;
    */

    const size_t real_scc_count = std::unordered_set(SCC_id.begin(), SCC_id.end()).size();

    // compare the real_scc_count to the known number of SCCs

    /*
    celegansneural: 57
    foldoc: 71
    language: 2456
    eu-2005: 90768
    wiki-topcats: 1
    sx-stackoverflow: 953658
    wb-edu: 4269022
    indochina-2004: 1749052
    uk-2002: 3887634
    arabic-2005: 4000414
    uk-2005: 5811041
    twitter7: 8044728
    */

    if(dataset_name == "celegansneural") {
        if (real_scc_count != 57) {
            std::cout << "ERROR: celegansneural found " << real_scc_count << " instead of 57" << std::endl;
        }
    } else if (dataset_name == "foldoc") {
        if (real_scc_count != 71) {
            std::cout << "ERROR: foldoc found " << real_scc_count << " instead of 71" << std::endl;
        }
    } else if (dataset_name == "language") {
        if (real_scc_count != 2456) {
            std::cout << "ERROR: language found " << real_scc_count << " instead of 2456" << std::endl;
        }
    } else if (dataset_name == "eu-2005") {
        if (real_scc_count != 90768) {
            std::cout << "ERROR: eu-2005 found " << real_scc_count << " instead of 90768" << std::endl;
        }
    } else if (dataset_name == "wiki-topcats") {
        if (real_scc_count != 1) {
            std::cout << "ERROR: wiki-topcats found " << real_scc_count << " instead of 1" << std::endl;
        }
    } else if (dataset_name == "sx-stackoverflow") {
        if (real_scc_count != 953658) {
            std::cout << "ERROR: sx-stackoverflow found " << real_scc_count << " instead of 953658" << std::endl;
        }
    } else if (dataset_name == "wb-edu") {
        if (real_scc_count != 4269022) {
            std::cout << "ERROR: wb-edu found " << real_scc_count << " instead of 4269022" << std::endl;
        }
    } else if (dataset_name == "indochina-2004") {
        if (real_scc_count != 1749052) {
            std::cout << "ERROR: indochina-2004 found " << real_scc_count << " instead of 1749052" << std::endl;
        }
    } else if (dataset_name == "uk-2002") {
        if (real_scc_count != 3887634) {
            std::cout << "ERROR: uk-2002 found " << real_scc_count << " instead of 3887634" << std::endl;
        }
    } else if (dataset_name == "arabic-2005") {
        if (real_scc_count != 4000414) {
            std::cout << "ERROR: arabic-2005 found " << real_scc_count << " instead of 4000414" << std::endl;
        }
    } else if (dataset_name == "uk-2005") {
        if (real_scc_count != 5811041) {
            std::cout << "ERROR: uk-2005 found " << real_scc_count << " instead of 5811041" << std::endl;
        }
    } else if (dataset_name == "twitter7") {
        if (real_scc_count != 8044728) {
            std::cout << "ERROR: twitter7 found " << real_scc_count << " instead of 8044728" << std::endl;
        }
    } else {
        std::cout << "ERROR: Unknown dataset " << dataset_name << "found " << real_scc_count << std::endl;
    }
}

int main(int argc, char** argv) {
    size_t times = 1;
    bool DEBUG = false;
    bool TOO_BIG = false;
    size_t NUM_THREADS = 1;

    if(argc == 1) {
        std::cout << "Usage: " << argv[0] << " [relativeFilePath: path to \".mtx\" file or folder] [timesToRun: int] [DEBUG: 1 or 0] [TOO_BIG: 1 or 0] [NUM_THREADS: 1...24]\n" << std::endl;

        std::cout << "Description:\n" << std::endl;
        std::cout << "    relativeFilePath:     The path to the file to run, relative to the matrices folder" << std::endl;
        std::cout << "    timesToRun:           The number of times to run the algorithm, the resulting time will be the average" << std::endl;
        std::cout << "    DEBUG:                If true, will print out debug information" << std::endl;
        std::cout << "    TOO_BIG:              If true, will skip the CSR conversion to save space and run the algorithm only with the CSC Matrix, which will be a little slower" << std::endl;
        std::cout << "    NUM_THREADS:          The number of threads to use in the pthreads implementation" << std::endl;
        std::cout << std::endl;

        std::cout << "Running with relativeFilePath not ending in \".mtx\" means this is the matrix folder the programm will try to run the algorithm on all files of the form:" << std::endl;
        std::cout << "      relativeFilePath/celegansneural/celegansneural.mtx,\n      relativeFilePath/foldoc/foldoc.mtx\n      ... etc ...\n      for all known datasets from the SPARSE MATIX COLLECTION." << std::endl;
        std::cout << std::endl;

        std::cout << "Example: " << argv[0] << DEFAULT_PREFIX << "language/language.mtx 10 1 1\n" << std::endl;
        std::cout << "      This will run the algorithm 10 times on the language matrix, printing out debug information and skipping the CSR conversion\n" << std::endl;

        std::cout << "Example: " << argv[0] << " " << DEFAULT_PREFIX << " 10 1 1\n" << std::endl;
        std::cout << "      This will run the algorithm 10 times on all files in the matrices folder in the way described above, printing out debug information and skipping the CSR conversion" << std::endl;
        std::cout << std::endl;

        std::cout << "Running test with default dataset: " << std::endl;
        std::cout << std::endl;
        testFile(DEFAULT_PREFIX + std::string("language/language.mtx"), times, DEBUG, TOO_BIG, NUM_THREADS);
        return 0;
    }

    if(argc > 2) {
        times = std::stoi(argv[2]);
    }

    if(argc > 3) {
        DEBUG = std::stoi(argv[3]) == 1;
    }

    if(argc > 4) {
        TOO_BIG = std::stoi(argv[4]) == 1;
    }

    if(argc > 5) {
        NUM_THREADS = std::stoi(argv[5]);
    }

    std::vector<std::string> filesToRun;
 
    if(argc > 1) {
        std::string inputFilename = argv[1];
        // check if it ends in .mtx
        if(inputFilename.substr(inputFilename.size() - 4) == ".mtx") {
            filesToRun = {inputFilename};
        } else {
            if(inputFilename.substr(inputFilename.size() - 1) == "/") {
                inputFilename = inputFilename.substr(0, inputFilename.size() - 1);
            }

            // files: celegansneural, foldoc, language, eu-2005, wiki-topcats, sx-stackoverflow, wb-edu, indochina-2004, uk-2002, arabic-2005, uk-2005, twitter7
            filesToRun = {
                inputFilename + "/celegansneural/celegansneural.mtx",
                inputFilename + "/foldoc/foldoc.mtx",
                inputFilename + "/language/language.mtx",
                inputFilename + "/eu-2005/eu-2005.mtx",
                inputFilename + "/wiki-topcats/wiki-topcats.mtx",
                inputFilename + "/sx-stackoverflow/sx-stackoverflow.mtx",
                inputFilename + "/wb-edu/wb-edu.mtx",
                inputFilename + "/indochina-2004/indochina-2004.mtx",
                inputFilename + "/uk-2002/uk-2002.mtx",
                inputFilename + "/arabic-2005/arabic-2005.mtx",
                inputFilename + "/uk-2005/uk-2005.mtx",
                inputFilename + "/twitter7/twitter7.mtx"
            };
        }
    }

    std::cout << "Starting with options: " << std::endl;
    std::cout << "Times: " << times << std::endl;
    std::cout << "Debug: " << DEBUG << std::endl;
    std::cout << "Too big: " << TOO_BIG << std::endl;

    std::cout << std::endl;

    for(auto& filename : filesToRun) {
        testFile(filename, times, DEBUG, TOO_BIG, NUM_THREADS);
    }

    std::cout << std::endl;
    std::cout << "---END---" << std::endl;

    return 0;
}