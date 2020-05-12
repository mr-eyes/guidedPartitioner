#include "colored_kDataFrame.hpp"
#include "algorithms.hpp"
#include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv) {

    if (argc != 3) {
        cerr << "run ./genesIndexing <fasta_file> <names_mape>" << endl;
        exit(1);
    }

    string fasta_file = argv[1];
    string names_file = argv[2];

    cout << "Indexing " << fasta_file << " with k=25" << endl;

    int kSize = 25;

    int chunk_size = 10000;

    auto *kf = new kDataFramePHMAP(kSize, 1);

    colored_kDataFrame *ckf = kProcessor::index(kf, {{"mode", 1}}, fasta_file, chunk_size, names_file);

    ckf->save(fasta_file);

}