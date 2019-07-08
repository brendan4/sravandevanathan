
import glob
import os

#fastq_folder  = "/home/hakan/scratch/dba/DBA" 
fastq_folder  = "/work/06059/ozadamh/lonestar/projects/DBA/fastq/merged"
output_folder = "/scratch/06059/ozadamh/dba/run_july_4" 


## Not used for reference purposes we keep it here
slurm_options_example = \
"""
#SBATCH -J Job%j           # Job name
#SBATCH -o /PATH_to_be_determined/out.o%j       # Name of stdout output file
#SBATCH -e /PATH_to_be_determined/err.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 16              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 12:30:00        # Run time (hh:mm:ss)
#SBATCH -A Ceniklab-CPRIT
"""

slurm_footer = \
"""#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 16              # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 12:30:00        # Run time (hh:mm:ss)
#SBATCH -A Ceniklab-CPRIT
"""


def generate_slurm_commands(exp_name, output_folder):
  slurm_str   = "#!/bin/bash\n"
  slurm_str  += "#SBATCH -J {}".format(exp_name)
  outdir_line = "#SBATCH -o {}/slurm_out/{}.out%j".format(output_folder, exp_name)
  errdir_line = "#SBATCH -e {}/slurm_err/{}.err%j".format(output_folder, exp_name)

  return "\n".join( (slurm_str,   outdir_line,
                     errdir_line, slurm_footer) )



base_script = \
"""
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

cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
          -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
          -A CAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTC \
          -u 6 -U 6 -q 22 -m 22 \
          -o $OUTFOLDER/clippedreads/${output1}/clipped_${output1}_R1.fastq.gz \
          -p $OUTFOLDER/clippedreads/${output1}/clipped_${output2}_R2.fastq.gz \
          --cores ${NP} \
          $input1 $input2

#############################################

###### F A S T Q C  ########################
fastqc -o $OUTFOLDER/fastqc/${output1}/ \
       $OUTFOLDER/clippedreads/${output1}/clipped_${output1}_R1.fastq.gz

fastqc -o $OUTFOLDER/fastqc/${output2}/ \
       $OUTFOLDER/clippedreads/${output1}/clipped_${output2}_R2.fastq.gz
###########################################

#### H I S A T 2 ##########################
hisat2 -p ${NP} -x ${HISAT2INDEX} \
     -1 $OUTFOLDER/clippedreads/${output1}/clipped_${output1}_R1.fastq.gz \
     -2 $OUTFOLDER/clippedreads/${output1}/clipped_${output2}_R2.fastq.gz \
     -S $OUTFOLDER/alignment/${output1}/${output1}_align.sam \
     --un-conc-gz $OUTFOLDER/alignment/${output1}/${output1}_unaligned.gz

#### S A M T O O L S  #####################

samtools view -@ ${NP} -bS $OUTFOLDER/alignment/${output1}/${output1}_align.sam \
   | samtools sort -@ ${NP} -o $OUTFOLDER/sorted/${output1}/${output1}_sorted.bam

samtools index -@ ${NP} -b $OUTFOLDER/sorted/${output1}/${output1}_sorted.bam

#### S T R I N G T I E #####################
stringtie -M 0.75  -e -p ${NP} -G ${CHESSGTF} \
   -b $OUTFOLDER/ballgown/${output1}/ \
   -A $OUTFOLDER/ballgown/${output1}/gene_abundance.tab \
   $OUTFOLDER/sorted/${output1}/${output1}_sorted.bam \
   -o $OUTFOLDER/misc/${output1}/n.gtf
"""


def get_experiment_folders():
   lev_1_list   = os.listdir( fastq_folder ) 

   exp_folders = []
    
   for seq_lane in lev_1_list:
        sub_folder           = os.path.join( fastq_folder, seq_lane ) 
        seq_lane_experiments = os.listdir( sub_folder )

        for experiment in seq_lane_experiments:
    	    exp_folder = os.path.join(sub_folder, experiment)
    	    exp_folders.append(exp_folder)

   return exp_folders


def generate_script(paired_files):
    if len(paired_files) < 2:
        return 0

    exp_name = paired_files[0].split("/")[-2]
    
    with open(exp_name + ".sh", "w") as output_stream:

        slurm_str = generate_slurm_commands(exp_name, output_folder)
        print(slurm_str, file = output_stream)

        print( "input1=" + paired_files[0], 
                file = output_stream )
        print( "input2=" + paired_files[1] + "\n\n\n", 
                file = output_stream)

        print(base_script, file = output_stream)



def main():
    exp_folders = get_experiment_folders()
    print(exp_folders)

    for e in exp_folders:
        r1_contents = glob.glob( e + "/*_R1*gz" )
        r2_contents = glob.glob( e + "/*_R2*gz" )
        e_contents = r1_contents + r2_contents
        generate_script(e_contents)
    
    
if __name__ == '__main__':
  main()

