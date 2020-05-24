import sys
import gzip


def auto_open(file_name, mode, encoding='utf8'):
    if file_name[-3:] == ".gz":
        return gzip.open(file_name, mode + 't', encoding=encoding)
    else:
        return open(file_name, mode)


fasta_file = str()
if len(sys.argv) < 2:
    sys.exit("run: python extract_interlaced.py <fasta>")
else:
    fasta_file = sys.argv[1]


file1_name = fasta_file.replace(".fa", '') + "_1.fa"
file2_name = fasta_file.replace(".fa", '') + "_2.fa"

with auto_open(fasta_file, 'r') as FASTA_READER, open(file1_name, 'w') as FA1, open(file2_name, 'w') as FA2:
    for line in FASTA_READER:
        seq1_header = line
        seq1_seq = next(FASTA_READER)
        seq2_header = next(FASTA_READER)
        seq2_seq = next(FASTA_READER)
        FA1.write(seq1_header + seq1_seq)
        FA2.write(seq2_header + seq2_seq)
