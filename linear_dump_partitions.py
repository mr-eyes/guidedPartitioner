import sqlite3
import sys
from collections import defaultdict
import os
from tqdm import tqdm
import multiprocessing


if len(sys.argv) < 4:
    raise ValueError("run: python dump_partitions.py <db_file> <idx_prefix> <no_cores>")

sqlite_db_path = sys.argv[1]
idx_prefix = sys.argv[2]
no_cores = int(sys.argv[3])


genes = dict()
total_genes_count = int()

with open(idx_prefix + ".namesMap") as namesMapReader:
    total_genes_count = int(next(namesMapReader))
    for line in namesMapReader:
        gene_id, gene_name = tuple(line.strip().split(' '))
        genes[int(gene_id)] = gene_name
    assert total_genes_count == len(genes)



output_dir = os.path.basename(sqlite_db_path).replace(".db", '')
os.makedirs(output_dir)

all_params = list()

for gene_id, gene_name in genes.items():
    file_name = os.path.join(output_dir, f"{gene_id}.fa")
    all_params.append((file_name, gene_id, gene_name))

conn = sqlite3.connect(sqlite_db_path)


for param in tqdm(all_params):
    file_path, gene_id, gene_name = param    
    read_sql = f"select * from reads where gene_id = {gene_id}"
    read_curs = conn.execute(read_sql)
    rows = read_curs.fetchall()
    if len(rows):
        with open(file_path, 'w') as fastaWriter:
            for row in rows:
                fastaWriter.write(f">{row[2]}_1\t{gene_name}\n{row[3]}\n")
                fastaWriter.write(f">{row[2]}_2\t{gene_name}\n{row[4]}\n")

conn.close()