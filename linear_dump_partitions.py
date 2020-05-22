import sqlite3
import sys
from collections import defaultdict
import os
from tqdm import tqdm
import multiprocessing
from itertools import islice


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

all_params = dict()

for gene_id, gene_name in genes.items():
    file_name = os.path.join(output_dir, f"{gene_id}.fa")
    all_params[gene_id] = (file_name, gene_name)

conn = sqlite3.connect(sqlite_db_path)

chunks = list()

def chunks(data, SIZE=500):
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k:data[k] for k in islice(it, SIZE)}

all_chunks = list(chunks(all_params, 800))

for chunk in tqdm(all_chunks):
    all_gene_ids = tuple(chunk.keys())
    geneID_to_geneName = dict()
    all_file_paths = set()
    written_files = set()

    files_handlers = dict()
    for gene_id, params in chunk.items():
        file_path, gene_name = params
        all_file_paths.add(file_path)
        geneID_to_geneName[gene_id] = gene_name
        files_handlers[gene_id] = open(file_path, 'w')

    read_sql = "select * from reads where gene_id in ({seq})".format(seq=','.join(['?'] * len(all_gene_ids)))
    curs = conn.execute(read_sql, all_gene_ids)
    rows = curs.fetchall()

    for row in rows:
        written_files.add(chunk[row[1]][0])
        files_handlers[row[1]].write(f">{row[2]}_1\n{row[3]}\n")
        files_handlers[row[1]].write(f">{row[2]}_2\n{row[4]}\n")

    for gene_id, fileHandler in files_handlers.items():
        fileHandler.close()

    empty_files = all_file_paths - written_files

    for _file in empty_files:
        os.remove(_file)
        

conn.close()