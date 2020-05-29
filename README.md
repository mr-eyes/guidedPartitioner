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

## 5. Assembling the whole data with **PLASS**

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

## 6. Whole data assembly assessment (RNAQuast) [min-contig-length=1000]

```bash

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz
gunzip *gz


THREADS=24
GTF=gencode.vM25.chr_patch_hapl_scaff.annotation.gtf
REF_GENOME=GRCm38.p6.genome.fa
TRANSCRIPTS=assembled_SRR11015356.fa

/usr/bin/time -v rnaQUAST.py --disable_infer_transcripts --disable_infer_genes --transcripts ${TRANSCRIPTS} --reference ${REF_GENOME} --gtf ${GTF} -t ${THREADS} -o rnaQuast_plass_SRR11015356

```

```txt
METRICS/TRANSCRIPTS                                    assembled_SRR11015356

 == DATABASE METRICS ==
Genes                                                  55938
Avg. number of exons per isoform                       5.909

 == BASIC TRANSCRIPTS METRICS ==
Transcripts                                            202959
Transcripts > 500 bp                                   202959
Transcripts > 1000 bp                                  202959

 == ALIGNMENT METRICS ==
Aligned                                                202750
Uniquely aligned                                       155842
Multiply aligned                                       1006
Unaligned                                              209

 == ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS ==
Avg. aligned fraction                                  0.941
Avg. alignment length                                  1870.677
Avg. mismatches per transcript                         11.155

 == ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS ==
Misassemblies                                          41522

 == ASSEMBLY COMPLETENESS (SENSITIVITY) ==
Database coverage                                      0.263
Duplication ratio                                      2.999
50%-assembled genes                                    11403
95%-assembled genes                                    4879
50%-covered genes                                      11815
95%-covered genes                                      7985
50%-assembled isoforms                                 22298
95%-assembled isoforms                                 5884
50%-covered isoforms                                   24694
95%-covered isoforms                                   10397
Mean isoform coverage                                  0.742
Mean isoform assembly                                  0.661

 == GeneMarkS-T METRICS ==
Predicted genes                                        104996

 == ASSEMBLY SPECIFICITY ==
50%-matched                                            40920
95%-matched                                            13965
Unannotated                                            0.63
```

```txt
Memory: 100GB
Cores: 24
Time: 3:38:25
```

## 7. Whole data assembly assessment (RNAQuast) [min-contig-length=500]

```txt
 == BASIC TRANSCRIPTS METRICS ==
Transcripts                                            476494
Transcripts > 500 bp                                   476494
Transcripts > 1000 bp                                  212853

 == ALIGNMENT METRICS ==
Aligned                                                475347
Uniquely aligned                                       375177
Multiply aligned                                       2951
Unaligned                                              1147

 == ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS ==
Avg. aligned fraction                                  0.947
Avg. alignment length                                  1182.064
Avg. mismatches per transcript                         8.429

 == ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS ==
Misassemblies                                          86475

 == ASSEMBLY COMPLETENESS (SENSITIVITY) ==
Database coverage                                      0.314
Duplication ratio                                      3.937
50%-assembled genes                                    12632
95%-assembled genes                                    5628
50%-covered genes                                      13299
95%-covered genes                                      9283
50%-assembled isoforms                                 29651
95%-assembled isoforms                                 7177
50%-covered isoforms                                   34346
95%-covered isoforms                                   13765
Mean isoform coverage                                  0.714
Mean isoform assembly                                  0.618
```

```txt
Memory: 102GB
Cores: 24
Time: 5:40:09
```

## 8. Genes partitions assembly

```bash
THREADS=4
TMP_DIR=tmp_plass
GENES_PARTITIONS=genes_partitions

OUTPUT=assembled_transcripts
mkdir ${OUTPUT}

TOTAL_FILES_NUMBER=$(ls -1q ${GENES_PARTITIONS}/* | wc -l)
COUNTER=${TOTAL_FILES_NUMBER}

echo -e "Processing ${COUNTER} Partitions...\n"

PLASS_LOG=plass_assembly.log
MERGED_TRANSCRIPTS=all_transcripts.fa

touch ${MERGED_TRANSCRIPTS}
touch ${PLASS_LOG}

for FASTA in ${GENES_PARTITIONS}/*
do

    # Print remaining partitions
    COUNTER=`expr $COUNTER - 1`
    echo "Processing (${FASTA}): ${COUNTER} / ${TOTAL_FILES_NUMBER}"

    seqs_no=$(grep ">" ${FASTA} | wc -l)

    # Skipping all paritions with number of reads < 5
    if (( seqs_no < 10)); then
        continue
    fi


    # Extract the interlaced Fasta and convert to Fastq
    GENE_ID=$(basename ${FASTA} .fa)
    python extract_interlaced.py ${FASTA}
    R1=${GENE_ID}_1.fastq
    R2=${GENE_ID}_2.fastq

    ${PLASS} nuclassemble -v 0 ${R1} ${R2} assembled_${GENE_ID}.fa ${TMP_DIR} --threads ${THREADS} >>${PLASS_LOG} 2>&1

    # Move the output assembled transcripts into ${OUTPUT}
    cat assembled_${GENE_ID}.fa >> ${MERGED_TRANSCRIPTS}
    mv assembled_${GENE_ID}.fa ${OUTPUT}

    # Remove temporary files
    rm -rf ${TMP_DIR}
    rm -rf ${R1} ${R2}
    rm -rf ${FASTA}

done
```

```txt
time: 1d 1h 38m 58s
memory: 466 M
```

### 8.1 Post-assembly (CDHIT-97%)

```bash
THREADS=32
WORD_SIZE=11
SIM=0.97
REF_FASTA=all_transcripts.fa
MAX_RAM_MB=150000
cd-hit-est -i ${REF_FASTA} -n ${WORD_SIZE} -c ${SIM} -o clusters_${SIM}_${OUT_PREFIX} -d 0 -T ${THREADS} -M ${MAX_RAM_MB}
```

```txt
Memory: 12.18G
Time: 1h
```

### 8.2 MiContigLength=1000 assembly assessment

```bash
THREADS=32
GTF=gencode.vM25.chr_patch_hapl_scaff.annotation.gtf
REF_GENOME=GRCm38.p6.genome.fa
TRANSCRIPTS=assembled_transcripts.fa
/usr/bin/time -v rnaQUAST.py --disable_infer_transcripts --disable_infer_genes --transcripts ${TRANSCRIPTS} --reference ${REF_GENOME} --gtf ${GTF} -t ${THREADS} -o rnaQuast_plass_SRR11015356
```

```txt
time: 3:48:12
memory: 13GB
```

> Something wrong in the total number of transcripts

```txt
METRICS/TRANSCRIPTS                                    assembled_transcripts

 == DATABASE METRICS ==
Genes                                                  55938
Avg. number of exons per isoform                       5.909

 == BASIC TRANSCRIPTS METRICS ==
Transcripts                                            4227
Transcripts > 500 bp                                   4227
Transcripts > 1000 bp                                  4227

 == ALIGNMENT METRICS ==
Aligned                                                4227
Uniquely aligned                                       183684
Multiply aligned                                       861
Unaligned                                              0

 == ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS ==
Avg. aligned fraction                                  0.827
Avg. alignment length                                  1603.279
Avg. mismatches per transcript                         10.745

 == ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS ==
Misassemblies                                          2546

 == ASSEMBLY COMPLETENESS (SENSITIVITY) ==
Database coverage                                      0.255
Duplication ratio                                      4.296
50%-assembled genes                                    10602
95%-assembled genes                                    4484
50%-covered genes                                      11292
95%-covered genes                                      6994
50%-assembled isoforms                                 20947
95%-assembled isoforms                                 5830
50%-covered isoforms                                   23334
95%-covered isoforms                                   9726
Mean isoform coverage                                  0.779
Mean isoform assembly                                  0.694

 == GeneMarkS-T METRICS ==
Predicted genes                                        3780

 == ASSEMBLY SPECIFICITY ==
50%-matched                                            2777
95%-matched                                            8
Unannotated                                            0.74

```

### 8.2 MiContigLength=500 assembly assessment

- Assembly

```txt
memory: 350 MB
time: 1d 5h 22m
```

- Assembly assessment

```txt
memory: 12.6 GB
time: 4:52:43
```

> The same postAssembly clustering and rnaQuast is done, just putting the results

```txt

METRICS/TRANSCRIPTS                                    assembled_transcripts

 == DATABASE METRICS ==
Genes                                                  55938
Avg. number of exons per isoform                       5.909

 == BASIC TRANSCRIPTS METRICS ==
Transcripts                                            14179
Transcripts > 500 bp                                   14179
Transcripts > 1000 bp                                  5554

 == ALIGNMENT METRICS ==
Aligned                                                14178
Uniquely aligned                                       387243
Multiply aligned                                       2029
Unaligned                                              1

 == ALIGNMENT METRICS FOR NON-MISASSEMBLED TRANSCRIPTS ==
Avg. aligned fraction                                  0.808
Avg. alignment length                                  944.198
Avg. mismatches per transcript                         7.516

 == ALIGNMENT METRICS FOR MISASSEMBLED (CHIMERIC) TRANSCRIPTS ==
Misassemblies                                          7182

 == ASSEMBLY COMPLETENESS (SENSITIVITY) ==
Database coverage                                      0.287
Duplication ratio                                      4.707
50%-assembled genes                                    11592
95%-assembled genes                                    4377
50%-covered genes                                      12944
95%-covered genes                                      7601
50%-assembled isoforms                                 26974
95%-assembled isoforms                                 5744
50%-covered isoforms                                   32861
95%-covered isoforms                                   11399
Mean isoform coverage                                  0.712
Mean isoform assembly                                  0.603

 == GeneMarkS-T METRICS ==
Predicted genes                                        12768

 == ASSEMBLY SPECIFICITY ==
50%-matched                                            9365
95%-matched                                            24
Unannotated                                            0.716
```
