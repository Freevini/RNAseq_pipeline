##================
##subsetting files
##================

ROOT=/home/bioinf/bioinf_data/43_sovi/Projects/Metagenomics_pipeline/Model_data/e_coli/05_trial_05

for LIB in `ls $ROOT/FASTQ | awk -F"." '{print $1}'`
  do

      echo "${LIB}"
      less $ROOT/FASTQ/"${LIB}".fastq.gz| head -1000000 > $ROOT/FASTQ/"${LIB}".fastq
      gzip $ROOT/FASTQ/"${LIB}".fastq

done


ls $ROOT/FASTQ | awk -F"." '{print $1}'

less FASTQ/SRR2135666.fastq.gz| head -1000000 > FASTQ/SRR2135666.fastq
gzip FASTQ/SRR2135666.fastq
