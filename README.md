# guidedPartitioner


## 1. Downloading and preparing the mouse transcriptome

```shell script

mkdir -p data

# Download
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz
gunzip gencode.vM25.transcripts.fa.gz
mv gencode.vM25.transcripts.fa data/


```

## 2. Indexing

- RAM: 7.5 GB
- Chunk size: 1k
- Time: 00:03:28
- Type: PHMAP

```shell script

# Creating namesFile
grep ">" data/gencode.vM25.transcripts.fa | cut -c2- |  awk -F'|' '{print $0"\t"$2}' > data/gencode.vM25.transcripts.fa.names

# Indexing

FASTA=data/gencode.vM25.transcripts.fa
NAMES=data/gencode.vM25.transcripts.fa.names

/usr/bin/time -v python genes_indexing.py ${FASTA} ${NAMES}

```

## 3. Partitioning

```text
Total fragments: 67954363

Reads stats:
	unique:             112025611
	ambiguous:          13475355
	chimeric:           2771915
	unmatched:          7635845

Fragments stats:
	unique:             57951403
	ambiguous:          4308078
	chimeric:           2738505
	unmatched:          2956377
```

- command : `./build/peReadsPart SRR11015356_1.fasta SRR11015356_2.fasta idx_gencode.v34.transcripts.fa`
- TIme: 52 mins
- RAM: 4.55 GB

---

## 4. Dumping gene partitions

```bash
DB=genes_partitions.db
IDX_PREFIX=idx_gencode.vM25.transcripts.fa
NO_CORES=1

/usr/bin/time -v python linear_dump_partitions.py ${DB} ${IDX_PREFIX} ${NO_CORES}
```

- **Results:**

```text
Elapsed (wall clock) time (h:mm:ss or m:ss): 4:46:40
Memory: 0.95 GB
```

---

## Assembling the whole data with **PLASS**

```bash
# Convert fasta file to fastq (Required for plass)
python fasta_to_fastq.py SRR11015356_1.fasta SRR11015356_2.fasta

R1=SRR11015356_1.fastq
R2=SRR11015356_2.fastq
THREADS=32
OUTPUT_DIR=tmp_SRR11015356_plass

/usr/bin/time -v plass nuclassemble -v 3 ${R1} ${R2} assembled_SRR11015356.fa ${OUTPUT_DIR} --threads ${THREADS}

```

- NO Cores: 32
- Memory: 162 GB
- Time: 6:27:08

#### **seqkit stats**

```tsv
file                      format  type  num_seqs      sum_len  min_len  avg_len  max_len
assembled_SRR11015356.fa  FASTA   DNA    202,959  391,239,783    1,000  1,927.7   11,144
```

---

