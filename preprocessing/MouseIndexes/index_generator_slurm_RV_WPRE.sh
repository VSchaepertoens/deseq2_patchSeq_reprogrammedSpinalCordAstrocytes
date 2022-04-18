#!/bin/bash

#SBATCH --ntasks=4
#SBATCH --mem=128000

#fetch indexes
#wget ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz
#wget ftp://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz
#gunzip *.gz

module load ngs/STAR/2.7.1a



STAR --runThreadN 4 \
 --limitGenomeGenerateRAM=33527195690 \
 --runMode genomeGenerate \
 --genomeDir /home/vbednarova/mouseIndexes/STARmouse_RV_WPRE \
 --genomeFastaFiles /home/vbednarova/mouseIndexes/Mus_musculus.GRCm38.dna_sm.toplevel.fa /home/vbednarova/mouseIndexes/RV_WPRE/RV_WPRE.fa 
 --sjdbGTFfile /home/vbednarova/mouseIndexes/RV_WPRE/Mus_musculus.GRCm38.100_RV_WPRE.gtf \
 --sjdbOverhang 49

 