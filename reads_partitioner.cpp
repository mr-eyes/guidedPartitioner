#include "colored_kDataFrame.hpp"
#include <string>
#include <iostream>
#include <vector>
#include "kseqReader.hpp"
#include <sys/stat.h>


using namespace std;

class statistics {
public:
    uint64_t n_unmatched;
    uint64_t n_same;
    uint64_t n_amb_same;
    uint64_t n_clear_fusion;
    uint64_t n_ambig_fusion;
    uint64_t n_mutli_fusion;
    uint64_t n_paired_fusion;
    uint64_t n_fragments;

    statistics() {
        n_fragments = 0;
        n_unmatched = 0;
        n_same = 0;
        n_amb_same = 0;
        n_clear_fusion = 0;
        n_ambig_fusion = 0;
        n_mutli_fusion = 0;
        n_paired_fusion = 0;
    }

    void print() const {
        cout << "No of input fragments: " << n_fragments << endl;
        cout << "unmatched:" << n_unmatched << endl;
        cout << "Unique:" << n_same << endl;
        cout << "Ambiguous:" << n_amb_same << endl;
        cout << "Single read clear fusion:" << n_clear_fusion << endl;
        cout << "Single read ambiguous fusion:" << n_ambig_fusion << endl;
        cout << "Single read multi fusion:" << n_mutli_fusion << endl;
        cout << "paired read fusion:" << n_paired_fusion << endl;

    }
};


inline bool file_exists (const std::string& name) {
    struct stat buffer{};
    return (stat (name.c_str(), &buffer) == 0);
}

static string
readFusion(string read, kDataFrame *DB, vector<vector<int> > *families, vector<int> *shared_kmers, vector<int> *gaps,
           statistics *stats);

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

int main(int argc, char **argv) {


    if (argc != 4) {
        cerr << "run ./readsPartitioner <R1> <R2> <index_prefix>" << endl;
        exit(1);
    }

    string PE_1_reads_file = argv[1];
    string PE_2_reads_file = argv[2];
    string index_prefix = argv[3];

    if (!file_exists(PE_1_reads_file)){
        throw std::runtime_error("Could not open R1 file");
    }

    if (!file_exists(PE_2_reads_file)){
        throw std::runtime_error("Could not open R2 file");
    }

    if (!file_exists(index_prefix + ".extra")){
        throw std::runtime_error("Could not open kProcessor index file");
    }

    int batchSize = 1000;
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


    auto *PEReader = new kseqReader(PE_1_reads_file, PE_2_reads_file, 1000);

    while (!PEReader->end()) {
        std::chrono::high_resolution_clock::time_point _currentTime = std::chrono::high_resolution_clock::now();
        cerr << "Processing chunk " << ++current_chunk << "/" << no_chunks << "... ";

        for (auto const &PE : *PEReader->next_chunk()) {

            // Analysis goes here

        }

        cerr << "Done in " << time_diff(_currentTime) << endl;
    }


    statistics stats;
    vector<vector<int> > families0,families1,families;
    vector<int> shared_kmers0,gaps0,shared_kmers1,gaps1,shared_kmers,gaps;

    return 0;

}

