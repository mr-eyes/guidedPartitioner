import sys
from Bio import SeqIO

if len(sys.argv) == 1:
    exit("run: python fasta_to_fastq.py <seq1> <seq2> <seq3> ....")

for fasta_file in sys.argv[1:]:
    records = SeqIO.parse(fasta_file, "fasta")
    count = SeqIO.write(records, fasta_file.replace(".fa","fq"), "fastq")
    print("Converted %i records" % count)