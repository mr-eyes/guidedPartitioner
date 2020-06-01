import sys
import os
import subprocess
from tqdm import tqdm

def get_lines_count(file_path):
    if file_path[-3:] == ".gz":
        ps = subprocess.Popen(
            f"gzip -cd {file_path} | wc -l", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        return int(ps.communicate()[0])

    else:
        return int(subprocess.getoutput('wc -l ' + file_path).split()[0])

fasta_file = sys.argv[1]
no_seqs = get_lines_count(fasta_file) // 2
fasta_output = "serialzed_" + os.path.basename(fasta_file)
serial_id = 1
with open(fasta_file, 'r') as FASTA_IN, open(fasta_output, 'w') as FASTA_OUT:
    for line in tqdm(FASTA_IN, total=no_seqs):
        header = line.replace(">",f">{serial_id} ")
        seq = next(FASTA_IN)
        FASTA_OUT.write(header + seq)
        serial_id += 1