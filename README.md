# guidedPartitioner

## 1. Indexing

<p>
    Indexing the human transcriptome by genes.
</p>

### 1.1 Downloading and preparing the human transcriptome

```shell script

mkdir -p data

# Download
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz
gunzip gencode.v34.transcripts.fa.gz
mv gencode.v34.transcripts.fa data/


```

## 2. Indexing

- RAM: 8.7 GB
- Chunk size: 1k
- Time: 00:04:50
- Type: PHMAP

```shell script

# Creating namesFile
grep ">" data/gencode.v34.transcripts.fa | cut -c2- |  awk -F'|' '{print $0"\t"$2}' > data/gencode.v34.transcripts.fa.names

# Indexing

FASTA=data/gencode.v34.transcripts.fa
NAMES=data/gencode.v34.transcripts.fa.names

/usr/bin/time -v python genes_indexing.py ${FASTA} ${NAMES}

```