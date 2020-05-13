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

## 3. Partitioning

### 3.1 Pilot run (Draft)

```
No of input fragments: 67954363
unmatched:89806872
Unique:33712972
Ambiguous:6548029
Single read clear fusion:1852642
Single read ambiguous fusion:1700609
Single read multi fusion:2287602
paired read fusion:988318
```
- command : `./build/peReadsPart SRR11015356_1.fasta SRR11015356_2.fasta idx_gencode.v34.transcripts.fa`
- TIme: 1:27:57
- RAM: 4.5 GB
