#!/bin/bash
# mamba activate /home/mfernand/miniforge3/envs/genomics_ONT
# Parse the input arguments
while [ "$#" -gt 0 ]; do
  case "$1" in
    --dir_reads=*)
        dir_reads="${1#*=}"
        ;;
    --output=*)
        d="${1#*=}"
        ;;
    --illumina=*)
        illumina_dir="${1#*=}"
        ;;
  esac
  shift
done

mkdir -p $d
subdirs="01.RawData 02.Porechop/Nanoplot 03.Chopper 04.Flye 05.Medaka 06.QC"
for subdir in $subdirs ; do
    mkdir -p $d/$subdir
done

cp $dir_reads/*.fastq.gz $d/01.RawData
echo "Reads in $dir_reads copied to $d/01.RawData"

for reads in $d/01.RawData/*.fastq.gz ; do
  sample=$(basename $reads | cut -d'.' -f1)
  # Remove adapters
  # Make sure adapters are in ~/miniforge3/envs/genomics_ont/lib/python3.12/site-packages/porechop/adapters.py
  base=$(basename "$file" .fastq.gz)

  # Run NanoPlot
  NanoPlot -t 72 --fastq $d/01.RawData/*.fastq.gz -o $d/01.RawData/Nanoplot

  # Remove adapters and low quality reads
  porechop -i $reads -o $d/02.Porechop/${sample}_porechopped.fastq.gz -t 40 -v 2
  chopper -i $d/02.Porechop/${sample}_porechopped.fastq.gz -q 7 -l 100 -t 40 | gzip > $d/03.Chopper/${sample}_chopped.fastq.gz

  # Check again QC stats
  NanoPlot -t 72 --fastq $d/03.Chopper/*.fastq.gz -o $d/03.Chopper/Nanoplot
  fastcat --histograms $d/03.Chopper/histograms -f $d/03.Chopper/Fastcat/fastcat_summary.txt \
    -r $d/03.Chopper/Fastcat/read_fastcat_summary.txt $d/03.Chopper/. | gzip > $d/03.Chopper/Fastcat/concat_reads.fastq.gz
  seqkit stats -Ta $d/03.Chopper/*.fastq.gz > $d/03.Chopper/lenght_stats.tsv

  if [ -z "$illumina_dir" ]; then
    # Assemble the reads
    flye --nano-raw $d/03.Chopper/${sample}_chopped.fastq.gz --out-dir $d/04.Flye/${sample} --asm-coverage 125 --threads 72 -g 15k

    # Medaka polishing
    medaka_consensus -i $d/03.Chopper/${sample}_chopped.fastq.gz -d $d/04.Flye/${sample}/assembly.fasta -o $d/05.Medaka/${sample} -t 72

    # Map the reads to the consensus
    sample=$(basename $bar | cut -d'_' -f1)
    reads=$(ls ${d}03.Chopper/$sp/${sample}_*.fastq.gz)
    ref=$d/05.Medaka/concensus.fasta

    # Map the reads to the consensus
    mkdir -p $d/06.QC/Samtools/${sample}
    minimap2 -ax map-ont -N 1 $d/05.Medaka/$sp/${sample}/consensus.fasta \
      $d/03.Chopper/$sp/${sample}_chopped.fastq.gz -t 72 | samtools view -Sbh -@ 72 | samtools sort -@ 72 -o $d/06.QC/Samtools/${sample}/reads_mapped_sorted.bam
    samtools index $d/06.QC/Samtools/${sample}/reads_mapped_sorted.bam
    samtools stats $d/06.QC/Samtools/${sample}/reads_mapped_sorted.bam > $d/06.QC/Samtools/${sample}/reads_mapped_sorted.bam.stats
    grep ^SN $d/06.QC/Samtools/${sample}/reads_mapped_sorted.bam.stats | cut -f 2- > $d/06.QC/Samtools/${sample}.stats

    # Quast
    mkdir -p $d/06.QC/Quast
    quast.py -o $d/06.QC/Quast/${sample} -m 0 -t 75 --circos --fragmented \
      --bam $d/06.QC/Samtools/${sample}/reads_mapped_sorted.bam $ref
  fi
done

if [ -n "$illumina_dir" ]; then
  mkdir -p $d/03.Chopper/merged_reads $d/07.Pilon/bbmap
  # Merge the reads
  zcat $d/03.Chopper/*.fastq.gz > $d/03.Chopper/merged_reads/chopped_reads.fastq.gz

  # Assemble the reads
  flye --nano-raw $d/03.Chopper/merged_reads/chopped_reads.fastq.gz --out-dir $d/04.Flye --asm-coverage 125 --threads 72 -g 15k

  # Polishing
  medaka_consensus -i $d/03.Chopper/merged_reads/chopped_reads.fastq.gz -d $d/04.Flye/assembly.fasta -o $d/05.Medaka -t 72

  # Map Illumina reads to the consensus
  short_reads=$(ls $illumina_dir/*R1*.fastq.gz | sed 's/R1/R#/')
  bbmap.sh \
      ref=$d/05.Medaka/consensus.fasta \
      in=$short_reads \
      out=$d/07.Pilon/bbmap/short_reads.bam \
      minid=0.76 \
      local="t" \
      nodisk=t trimreaddescriptions=t secondarycov=f slow=t unpigz=t \
      covstats=$d/07.Pilon/bbmap/contig_covstats.txt \
      2> $d/07.Pilon/bbmap/bbmap.log

  samtools sort $d/07.Pilon/bbmap/short_reads.bam -o $d/07.Pilon/bbmap/short_reads_sorted.bam
  samtools $d/index 07.Pilon/bbmap/short_reads_sorted.bam

  # Polish with Pilon with Illumina reads
  pilon --genome $d/05.Medaka/consensus.fasta --frags $d/07.Pilon/bbmap/short_reads_sorted.bam --outdir $d/07.Pilon -Xmx16G
fi
