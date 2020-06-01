import sys
import subprocess
import gzip
from tqdm import tqdm


def get_lines_count(file_path):
    if file_path[-3:] == ".gz":
        ps = subprocess.Popen(
            f"gzip -cd {file_path} | wc -l", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return int(ps.communicate()[0])

    else:
        return int(subprocess.getoutput('wc -l ' + file_path).split()[0])


if len(sys.argv) == 1:
    exit("run: python fasta_to_fastq.py <seq1> <seq2> <seq3> ....")


for fasta_file in sys.argv[1:]:
    no_seqs = get_lines_count(fasta_file) // 2
    fastq_file = fasta_file[:-1] + "q"
    read_len = int()
    with open(fasta_file) as fastaReader:
        next(fastaReader)
        read_len = len(next(fastaReader).strip())

    fakeQuality = "+\n" + "?"*read_len + "\n"

    with open(fasta_file) as fastaReader, open(fastq_file, 'w') as fastqWriter:
        for line in tqdm(fastaReader, total=no_seqs):
            seq1_header = line.strip().replace(">", "@")[:-2] + "\n"
            seq1_seq = next(fastaReader)
            fastqWriter.write(seq1_header + seq1_seq + fakeQuality)
