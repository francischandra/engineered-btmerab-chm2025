#!/usr/bin/env bash

# DATASET QC AND TRIMMING

GENOME_DIR="genome"
RAW_DIR="raw"
QC_DIR="qc"
FILT_DIR="filt"

# fastqc to generate report

mkdir -p ${QC_DIR}/raw ${QC_DIR}/trim # this is to make the relevant directories for QC

fastqc $RAW_DIR/*.fastq.gz -o $QC_DIR/raw/. -t 24 # this code is to get the sequencing QC for all of the raw reads

# consolidate qc with multiqc

multiqc $QC_DIR/raw -o $QC_DIR/raw

# trim reads with fastp
# fastp can be used with paired end reads. Supply forward w/ -i and reverse w/ -I; output forward -o and rev -O

for file in $RAW_DIR/*R1_001.fastq.gz; do
        echo $file
        base=$(basename $file _R1_001.fastq.gz) #extract the sample name while removing the '_R1_001.fastq.gz'
        echo $base
        fastp -i $RAW_DIR/${base}_R1_001.fastq.gz -I $RAW_DIR/${base}_R2_001.fastq.gz \
              -o $FILT_DIR/${base}_R1.trimmed.fastq.gz -O $FILT_DIR/${base}_R2.trimmed.fastq.gz \
              --detect_adapter_for_pe -g -x -l 25 \ #trim poly g, poly x, reads with less than 25 bp
              --thread 16 \
              -j $FILT_DIR/${base}.fastp.json -h $FILT_DIR/${base}.fastp.html
done

# qc of trimmed reads

fastqc $FILT_DIR/*.trimmed.fastq.gz -o $QC_DIR/trim/. -t 24

# consolidate trim qc

multiqc $QC_DIR/trim -o $QC_DIR/trim

# GENOME INDEXING
# genome fasta and gtf file was downloaded from ENCODE database primary assembly (GRCm39)

STAR --runMode genomeGenerate --genomeDir /path/to/output/genome/dir --genomeFastaFiles /path/to/genome/fasta --sjdbGTFfile path/to/annotation/gtf --sjdbOverhang 149

# STAR ALIGNMENT

# Set up directories
FILT_DIR="filt"
ALIGN_DIR="align"

# Adjust the paths and parameters according to your specific STAR index and settings
STAR_INDEX="/path/to/indexed/genome"

for file in $FILT_DIR/*_R1.trimmed.fastq.gz; do
    base=$(basename $file _R1.trimmed.fastq.gz)
    STAR --genomeDir $STAR_INDEX \
         --runThreadN 24 \
         --readFilesIn $FILT_DIR/${base}_R1.trimmed.fastq.gz $FILT_DIR/${base}_R2.trimmed.fastq.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix $ALIGN_DIR/${base}_ \
         --quantTranscriptomeSAMoutput BanSingleEnd \
         --outSAMtype BAM Unsorted \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --outFilterMultimapNmax 20 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --quantMode TranscriptomeSAM \
         --outSAMattributes NH HI AS NM MD
done

# READ QUANTIFICATION
# generate transcriptome from genome fasta and gtf
gffread -w GRCm39.primary_assembly.transcripts.fa -g GRCm39.primary_assembly.genome.fa gencode.vM35.primary_assembly.annotation.gtf

# quantifying transcripts with Salmon

# Set up directories
ALIGN_DIR="align"
QUANT_DIR="quant"

# Adjust the paths and parameters according to your specific genome_dir and settings
GENOME_INDEX="/path/to/indexed/genome"

for file in $ALIGN_DIR/*_Aligned.toTranscriptome.out.bam; do
        base=$(basename $file _Aligned.toTranscriptome.out.bam)
        salmon quant -t $GENOME_INDEX/GRCm39.primary_assembly.transcripts.fa \
        --libType A \
        -a $ALIGN_DIR/${base}_Aligned.toTranscriptome.out.bam \
        -o $QUANT_DIR/${base}.salmon_quant \
        -p 20 \
        --gcBias --seqBias

done
