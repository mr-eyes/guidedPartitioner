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

### 3.1 Without DB Insertion

```
Total fragments: 67954363

Reads stats:
	unmatched:          89806872
	unique:             33714348
	ambiguous unique:   6545131
	clear chimeric:     1852951
	ambiguous chimeric: 1699724
	multi chimeric:     2289700

Fragments stats:
	unmatched:          51951662
	chimeric:           6024413
	unique:             8175096
	ambiguous unique:   1803192
```

- command : `./build/peReadsPart SRR11015356_1.fasta SRR11015356_2.fasta idx_gencode.v34.transcripts.fa`
- TIme: 1:26:30
- RAM: 4.5 GB
