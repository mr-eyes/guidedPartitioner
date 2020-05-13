#include "colored_kDataFrame.hpp"
#include <string>
#include <iostream>
#include <vector>
#include "kseqReader.hpp"
#include <sys/stat.h>
#include <map>
#include "readClassifier.hpp"

using namespace std;

inline bool file_exists(const std::string &name) {
    struct stat buffer{};
    return (stat(name.c_str(), &buffer) == 0);
}


inline string time_diff(std::chrono::high_resolution_clock::time_point &t1) {
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    long hr = milli / 3600000;
    milli = milli - 3600000 * hr;
    long min = milli / 60000;
    milli = milli - 60000 * min;
    long sec = milli / 1000;
    milli = milli - 1000 * sec;
    string timeDiff;
    timeDiff.append(to_string(min));
    timeDiff.append(":");
    timeDiff.append(to_string(sec));
    timeDiff.append(":");
    timeDiff.append(to_string(milli));

    return timeDiff;
}

inline bool isFusion(string flag) {
    return flag == "clear_fusion" || flag == "ambig_fusion" || flag == "multi_fusion";
}

inline bool isSameRef(string flag) {
    return flag == "unique" || flag == "ambiguous";
}


int main(int argc, char **argv) {


    if (argc != 4) {
        cerr << "run ./peReadsPart <R1> <R2> <index_prefix>" << endl;
        exit(1);
    }

    string PE_1_reads_file = argv[1];
    string PE_2_reads_file = argv[2];
    string index_prefix = argv[3];

    if (!file_exists(PE_1_reads_file)) {
        throw std::runtime_error("Could not open R1 file");
    }

    if (!file_exists(PE_2_reads_file)) {
        throw std::runtime_error("Could not open R2 file");
    }

    if (!file_exists(index_prefix + ".extra")) {
        throw std::runtime_error("Could not open kProcessor index file");
    }

    int batchSize = 10000;
    int kSize = 25;
    int no_of_sequences = 67954363;

    // kProcessor Index Loading
    std::cerr << "Loading kProcessor index ..." << std::endl;
    colored_kDataFrame *ckf = colored_kDataFrame::load(index_prefix);
    kDataFrame *kf = ckf->getkDataFrame();
    assert(kSize == (int) kf->getkSize());
    std::cerr << "kProcessor index loaded successfully ..." << std::endl;


    // Initializations
    int no_chunks = ceil((double) no_of_sequences / (double) batchSize);
    int current_chunk = 0;


    auto *PEReader = new kseqReader(PE_1_reads_file, PE_2_reads_file, batchSize);


    ReadsClassifier classifier(ckf);

    while (!PEReader->end()) {
        std::chrono::high_resolution_clock::time_point _currentTime = std::chrono::high_resolution_clock::now();
        cerr << "Processing chunk " << ++current_chunk << "/" << no_chunks << "... ";

        for (auto const &PE : *PEReader->next_chunk()) {
            classifier.stats.n_fragments++;

            tuple<string, string, string> fragementFlags = classifier.classifyFragment(PE.R1_seq, PE.R2_seq);

            string flag_R1 = get<0>(fragementFlags);
            string flag_R2 = get<1>(fragementFlags);
            string flag_fragment = get<2>(fragementFlags);


//            cout << PE.R1_name << " : " << flag_R1 << endl;
//            cout << PE.R2_name << " : " << flag_R2 << endl;



        }

        cerr << "Done in " << time_diff(_currentTime) << endl;
    }

    classifier.stats.print();

    return 0;

}

