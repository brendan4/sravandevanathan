#!/bin/bash
#SBATCH -J 121317_LIB8-98401312
#SBATCH -o /scratch/06059/ozadamh/dba/run_july_4/slurm_out/121317_LIB8-98401312.out%j
#SBATCH -e /scratch/06059/ozadamh/dba/run_july_4/slurm_err/121317_LIB8-98401312.err%j
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 16              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 12:30:00        # Run time (hh:mm:ss)
#SBATCH -A Ceniklab-CPRIT

input1=/work/06059/ozadamh/lonestar/projects/DBA/fastq/merged/DBA_121317/121317_LIB8-98401312/121317_LIB8-98401312_R1_001.fastq.gz
input2=/work/06059/ozadamh/lonestar/projects/DBA/fastq/merged/DBA_121317/121317_LIB8-98401312/121317_LIB8-98401312_R2_001.fastq.gz




output1=$(basename $input1| awk -F "_" '{print($1"_"$2)}')
output2=$(basename $input2| awk -F "_" '{print($1"_"$2)}')


OUTFOLDER=${SCRATCH}/dba/run_july_4/pipeline

HISAT2INDEX="/work/06059/ozadamh/lonestar/reference/human/HISAT2/grch38_snp_tran/genome_snp_tran"
CHESSGTF="/work/06059/ozadamh/lonestar/projects/DBA/reference/chessmod.gtf"
NP=16

source activate dba
###############################################################################

mkdir -p  $OUTFOLDER/fastqc/${output1}
mkdir -p  $OUTFOLDER/clippedreads/${output1}
mkdir -p  $OUTFOLDER/alignment/${output1}
mkdir -p  $OUTFOLDER/sorted/${output1}
#echo $output1
#echo $output2

### C U T A D A P T  #########################

cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG           -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG           -A CAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTC           -u 6 -U 6 -q 22 -m 22           -o $OUTFOLDER/clippedreads/${output1}/clipped_${output1}_R1.fastq.gz           -p $OUTFOLDER/clippedreads/${output1}/clipped_${output2}_R2.fastq.gz           --cores ${NP}           $input1 $input2

#############################################

###### F A S T Q C  ########################
fastqc -o $OUTFOLDER/fastqc/${output1}/        $OUTFOLDER/clippedreads/${output1}/clipped_${output1}_R1.fastq.gz

fastqc -o $OUTFOLDER/fastqc/${output2}/        $OUTFOLDER/clippedreads/${output1}/clipped_${output2}_R2.fastq.gz
###########################################

#### H I S A T 2 ##########################
hisat2 -p ${NP} -x ${HISAT2INDEX}      -1 $OUTFOLDER/clippedreads/${output1}/clipped_${output1}_R1.fastq.gz      -2 $OUTFOLDER/clippedreads/${output1}/clipped_${output2}_R2.fastq.gz      -S $OUTFOLDER/alignment/${output1}/${output1}_align.sam      --un-conc-gz $OUTFOLDER/alignment/${output1}/${output1}_unaligned.gz

#### S A M T O O L S  #####################

samtools view -@ ${NP} -bS $OUTFOLDER/alignment/${output1}/${output1}_align.sam    | samtools sort -@ ${NP} -o $OUTFOLDER/sorted/${output1}/${output1}_sorted.bam

samtools index -@ ${NP} -b $OUTFOLDER/sorted/${output1}/${output1}_sorted.bam

#### S T R I N G T I E #####################
stringtie -M 0.75  -e -p ${NP} -G ${CHESSGTF}    -b $OUTFOLDER/ballgown/${output1}/    -A $OUTFOLDER/ballgown/${output1}/gene_abundance.tab    $OUTFOLDER/sorted/${output1}/${output1}_sorted.bam    -o $OUTFOLDER/misc/${output1}/n.gtf

