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


inline bool file_exists(const std::string &name) {
    struct stat buffer{};
    return (stat(name.c_str(), &buffer) == 0);
}

static string readFusion(const string &read, colored_kDataFrame *DB, vector<vector<uint32_t> > &families,
                         vector<uint32_t> &shared_kmers, vector<uint32_t> &gaps, statistics &stats);

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

    if (!file_exists(PE_1_reads_file)) {
        throw std::runtime_error("Could not open R1 file");
    }

    if (!file_exists(PE_2_reads_file)) {
        throw std::runtime_error("Could not open R2 file");
    }

    if (!file_exists(index_prefix + ".extra")) {
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


    auto *PEReader = new kseqReader(PE_1_reads_file, PE_2_reads_file, batchSize);

    statistics stats;
    vector<vector<uint32_t> > families0, families1, families;
    vector<uint32_t> shared_kmers0, gaps0, shared_kmers1, gaps1, shared_kmers, gaps;


    while (!PEReader->end()) {
        std::chrono::high_resolution_clock::time_point _currentTime = std::chrono::high_resolution_clock::now();
        cerr << "Processing chunk " << ++current_chunk << "/" << no_chunks << "... ";

        for (auto const &PE : *PEReader->next_chunk()) {

            families0.clear();
            gaps0.clear();
            shared_kmers0.clear();
            families1.clear();
            gaps1.clear();
            shared_kmers1.clear();

            string flag_R1 = readFusion(PE.R1_seq, ckf, families0, shared_kmers0, gaps0, stats);
            string flag_R2 = readFusion(PE.R2_seq, ckf, families1, shared_kmers1, gaps1, stats);

//            cout << PE.R1_name << " : " << flag_R1 << endl;
//            cout << PE.R2_name << " : " << flag_R2 << endl;


        }

        cerr << "Done in " << time_diff(_currentTime) << endl;
    }


    return 0;

}

inline void print(const string &x) {
    cout << x << endl;
}

inline void dumpVec(const vector<uint32_t> &vec) {
    for (auto const &element : vec) {
        cout << element << ", ";
    }
    cout << endl;
}

static string readFusion(const string &read, colored_kDataFrame *DB, vector<vector<uint32_t> > &families,
                         vector<uint32_t> &shared_kmers, vector<uint32_t> &gaps, statistics &stats) {
    shared_kmers.clear();
    families.clear();
    gaps.clear();
    string flag;
    vector<uint32_t> lf_ids;
    vector<uint32_t> rt_ids;
    uint8_t kSize = DB->getkDataFrame()->ksize();
    if (read.size() < kSize) {
        stats.n_unmatched++;
        return "unmatched";
    }
    //# find a matching k-mer at the beginning of the read

    lf_ids = DB->getKmerSourceFromColor(DB->getkDataFrame()->getCount(read.substr(0, kSize)));

    int idx = 1;

    while (idx < read.size() - kSize + 1 && lf_ids.empty()) {
        lf_ids = DB->getKmerSourceFromColor(DB->getkDataFrame()->getCount(read.substr(idx, kSize)));
        idx++;
    }
    if (lf_ids.empty()) {
//        print("no single match");
        stats.n_unmatched++;
        flag = "unmatched";
    } else if (idx == read.size() - kSize + 1) {
//        print("same, only last kmer matched");
        families.push_back(lf_ids);

        if (lf_ids.size() == 1) {
            stats.n_same += 1;
            flag = "unique";
        } else {
            stats.n_amb_same += 1;
            flag = "ambiguous";
        }
    } else {
        // # len(lf_ids) > 0 & idx < len(hashvals)
        //# find a matching k-mer at the end of the read
        vector<uint32_t> rt_ids = DB->getKmerSourceFromColor(
                DB->getkDataFrame()->getCount(read.substr(read.size() - kSize, kSize)));

        int idy = read.size() - 1;
        while (idy - kSize >= idx - 1 && rt_ids.empty()) {
            rt_ids = DB->getKmerSourceFromColor(DB->getkDataFrame()->getCount(read.substr(idy - kSize, kSize)));
            idy--;
        }

        if (rt_ids.empty()) {
//            print("same, only one non-last kmer matched");
            families.push_back(lf_ids);
            if (lf_ids.size() == 1) {
                stats.n_same += 1;
                flag = "unique";
            } else {
                stats.n_amb_same += 1;
                flag = "ambiguous";
            }
        } else {
            vector<uint32_t> intersect_ids;
            intersect_ids.clear();
            auto it = set_intersection(lf_ids.begin(), lf_ids.end(), rt_ids.begin(), rt_ids.end(),
                                       back_inserter(intersect_ids));
            if (!intersect_ids.empty()) {
                families.push_back(intersect_ids);
                if (intersect_ids.size() == 1) {
                    stats.n_same += 1;
                    flag = "unique";
                } else {
                    stats.n_amb_same += 1;
                    flag = "ambiguous";
                }
            } else {
                //# fusion to be resolved
                uint64_t shared_kmer = 1;
                uint64_t gap_size = 0;
                bool gap = false;
                while (idx <= (idy + 1 - kSize)) {
                    vector<uint32_t> temp_ids = DB->getKmerSourceFromColor(
                            DB->getkDataFrame()->getCount(read.substr(idx, kSize)));

                    if (read.substr(idx, kSize).size() != kSize) {
                        cout << read.substr(idx, kSize) << " " << idx << " < " << read.size() << endl;
                    }
                    if (!temp_ids.empty()) {
                        intersect_ids.clear();
                        it = set_intersection(lf_ids.begin(), lf_ids.end(), temp_ids.begin(), temp_ids.end(),
                                              back_inserter(intersect_ids));

                        if (!intersect_ids.empty()) {
                            lf_ids = intersect_ids;
                            shared_kmer += 1;
                            gap_size = 0;
                        } else {
                            // # len(intersect_ids) == 0

                            families.push_back(lf_ids);
                            shared_kmers.push_back(shared_kmer);
                            lf_ids = temp_ids;
                            shared_kmer = 1;
                            gaps.push_back(gap_size);
                            gap_size = 0;
                        }
                    } else {
                        gap_size += 1;
                    }
                    idx += 1;
                }
                families.push_back(lf_ids);

                shared_kmers.push_back(shared_kmer);

                if (families.size() <= 1) {
                    cerr << "Error" << endl;
                    return "Error";
                }
                if (families.size() == 2) {
                    if ((families)[0].size() == 1 && (families)[1].size() == 1) {
                        stats.n_clear_fusion += 1;
                        flag = "clear_fusion";
                    } else {
                        stats.n_ambig_fusion += 1;
                        flag = "ambig_fusion";
                    }
                } else {
                    //# len(families) > 2
                    stats.n_mutli_fusion += 1;
                    flag = "multi_fusion";
                }

            }
        }
    }

    //#if len(families) == 0:
    //#    families = "-"

    //#if len(shared_kmers) == 0:
    //#    shared_kmers = "-"

    return flag;

}


