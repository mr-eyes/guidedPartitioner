#include <cstdint>
#include "colored_kDataFrame.hpp"

class Statistics {
public:
    uint64_t n_unmatched;
    uint64_t n_same;
    uint64_t n_amb_same;
    uint64_t n_clear_fusion;
    uint64_t n_ambig_fusion;
    uint64_t n_mutli_fusion;
    uint64_t n_paired_fusion;
    uint64_t n_paired_unique;
    uint64_t n_fragments;
    uint64_t n_paired_ambig_unique;
    uint64_t n_paired_unmatched;

    Statistics() {
        n_fragments = 0;
        n_unmatched = 0;
        n_same = 0;
        n_amb_same = 0;
        n_clear_fusion = 0;
        n_ambig_fusion = 0;
        n_mutli_fusion = 0;
        n_paired_fusion = 0;
        n_paired_unique = 0;
        n_paired_ambig_unique = 0;
        n_paired_unmatched = 0;
    }

    void print() const {
        cout << "Summary report" << endl;
        cout << "===============" << endl << endl;
        cout << "Total fragments: " << n_fragments << endl << endl;

        cout << "Reads stats:" << endl;

        cout << "\tunmatched:          " << n_unmatched << endl;
        cout << "\tunique:             " << n_same << endl;
        cout << "\tambiguous unique:   " << n_amb_same << endl;
        cout << "\tclear chimeric:     " << n_clear_fusion << endl;
        cout << "\tambiguous chimeric: " << n_ambig_fusion << endl;
        cout << "\tmulti chimeric:     " << n_mutli_fusion << endl;

        cout << endl;

        cout << "Fragments stats:" << endl;

        cout << "\tunmatched:          " << n_paired_unmatched << endl;
        cout << "\tchimeric:           " << n_paired_fusion << endl;
        cout << "\tunique:             " << n_paired_unique << endl;
        cout << "\tambiguous unique:   " << n_paired_ambig_unique << endl;
    }
};


class ReadsClassifier {

private:

    colored_kDataFrame *DB;
    int kSize;
    unordered_map<int, string> namesMap;

    vector<vector<uint32_t> > families0, families1, families;
    vector<uint32_t> shared_kmers0, gaps0, shared_kmers1, gaps1, shared_kmers, gaps;

public:
    Statistics stats;

    explicit ReadsClassifier(colored_kDataFrame *ckf) {
        this->DB = ckf;
        this->kSize = ckf->getkDataFrame()->ksize();
        namesMap = ckf->names_map();
    }

    static bool isFusion(string &flag) {
        return flag == "clear_fusion" || flag == "ambig_fusion" || flag == "multi_fusion";
    }

    static bool isSameRef(string &flag) {
        return flag == "unique" || flag == "ambiguous";
    }

    tuple<string, string, string> classifyFragment(const string &R1_seq, const string &R2_seq) {

        families0.clear();
        gaps0.clear();
        shared_kmers0.clear();
        families1.clear();
        gaps1.clear();
        shared_kmers1.clear();

        string flag_R1 = readFusion(R1_seq, families0, shared_kmers0, gaps0);
        string flag_R2 = readFusion(R2_seq, families1, shared_kmers1, gaps1);
        string fragmentFlag{};

        /*
         * unique, ambiguous, unmatched
         * multi_fusion, ambig_fusion, clear_fusion
         * */


        if (isFusion(flag_R1) || isFusion(flag_R2)) {
            // One or both are chimeric
            fragmentFlag = "chimeric";
            stats.n_paired_fusion++;

        }else if (isSameRef(flag_R1) && isSameRef(flag_R2)) {
            // Both are ambiguous or unique but may be no intersection
            vector<int> intersect_ids;
            intersect_ids.clear();
            auto it = set_intersection(families0[0].begin(), families0[0].end(),
                                       families1[0].begin(), families1[0].end(), back_inserter(intersect_ids));
            if (intersect_ids.empty()) {
                // If no intersection, then they are chimeric
                stats.n_paired_fusion++;
                fragmentFlag = "chimeric";
            } else {
                // Check if they are unique or ambiguous
                if (flag_R1 == "unique" && flag_R2 == "unique") {
                    fragmentFlag = "unique";
                    stats.n_paired_unique++;
                }else{
                    fragmentFlag = "ambiguous";
                    stats.n_paired_ambig_unique++;
                }
            }
        }else{
            fragmentFlag = "unmatched";
            stats.n_paired_unmatched++;
        }

        return make_tuple(flag_R1, flag_R2, fragmentFlag);
    }

    string readFusion(const string &read, vector<vector<uint32_t> > &families,
                      vector<uint32_t> &shared_kmers, vector<uint32_t> &gaps) {

        shared_kmers.clear();
        families.clear();
        gaps.clear();
        string flag;
        vector<uint32_t> lf_ids;
        vector<uint32_t> rt_ids;

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
        return flag;
    }

};