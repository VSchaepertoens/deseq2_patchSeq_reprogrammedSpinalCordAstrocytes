#!/bin/bash

#SBATCH --ntasks=2
#SBATCH --mem=80000

module load ngs/STAR/2.7.1a

STAR --runThreadN 2 \
	--genomeDir /home/vbednarova/mouseIndexes/STARmouse_AAV_1108 \
	--readFilesCommand gunzip -c \
	--readFilesIn /home/vbednarova/lafugaMay2020/rawdata/Goetz/200617_L183_0585_ACE5A3ANXX/Lane5_demultiplexed_reads/Sample${i}.txt.gz \
	--outFileNamePrefix /home/vbednarova/lafugaDecember2019/STARaligneddata/Sample${i} \
	--quantMode TranscriptomeSAM GeneCounts \