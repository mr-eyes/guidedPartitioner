# guidedPartitioner


### 1. Downloading and preparing the mouse transcriptome

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
