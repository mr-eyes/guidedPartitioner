import sys
import gzip
import os

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

basename = os.path.basename(fasta_file)
file1_name = basename.replace(".fa", '') + "_1.fastq"
file2_name = basename.replace(".fa", '') + "_2.fastq"

with auto_open(fasta_file, 'r') as FASTA_READER, open(file1_name, 'w') as FA1, open(file2_name, 'w') as FA2:
    for line in FASTA_READER:
        fakeQuality = "+\n" + "?"*151 + "\n"
        seq1_header = line.replace(">","@").replace("_1","")
        seq1_seq = next(FASTA_READER)
        seq2_header = next(FASTA_READER).replace(">","@").replace("_2","")
        seq2_seq = next(FASTA_READER)
        FA1.write(seq1_header + seq1_seq + fakeQuality)
        FA2.write(seq2_header + seq2_seq + fakeQuality)
