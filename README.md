## 
This readme.md file was generated on 15-04-2022 by Veronika Schäpertöns


GENERAL INFORMATION

1. Title of Dataset : Differential gene expression analysis of Ngn2 and Ascl1 reprogrammed spinal cord astrocytes

2. Author Information
	A. Principal Investigator Contact Information
		Name: Prof. Dr. Magdalena Götz 
		Institution: Ludwig-Maximilians-Universität München, Department of Physiological Genomics, BioMedical Center - BMC
		Address: Großhaderner Str. 9, D-82152 Planegg-Martinsried 
		
	B. Associate or Co-investigator Contact Information
		Name: Dr. Giacomo Masserdotti 
		Institution: Ludwig-Maximilians-Universität München, Department of Physiological Genomics, BioMedical Center - BMC
		Address: Großhaderner Str. 9, D-82152 Planegg-Martinsried 
		Email: Giacomo.Masserdotti@bmc.med.lmu.de

	C. Associate or Co-investigator Contact Information
		Name: Dr. Veronika Schäpertöns
		Institution: Ludwig-Maximilians-Universität München, Department of Physiological Genomics, BioMedical Center - BMC 
		Address: Großhaderner Str. 9, D-82152 Planegg-Martinsried 
		Email: vschaepertoens@protonmail.com

3. Date of data collection : March-April 2020

4. Geographic location of data collection : Munich, Germany

5. Information about funding sources that supported the collection of the data: German Research Foundation grants SFB 870 , SPP1757, the advanced ERC ChroNeuroRepair, ERA-Net neuron grant MICRONET, the EU consortium NSC Reconstruct, German Research Foundation grant TRR274


SHARING/ACCESS INFORMATION

1. Links to publications that cite or use the data : https://doi.org/10.1016/j.celrep.2021.109409

2. Links to other publicly accessible locations of the data : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173978


DATA & FILE OVERVIEW

1. Folders/files List: 
- __preprocessing__ 

  - __rawdata__ - folder containing raw fastq files e.g. lane5_R1.fastq.gz.fastqsanger.gz (folder empty)
  - __Demultiplex__ - folder containing barcodes *Indexes_08072020_2.txt* of each single patch-seq sample, script *script_to_demultiplex_08072020_6.sh* used to demultiplex dataset on lane5 and stats file after demultiplexing *jemultiplexer_out_stats.txt*
  - __MouseIndexes__ - folder containing script *index_generator_slurm_RV_WPRE.sh* used to build reference genome for Mus_musculus.GRCm38 and gene annotation table *geneInfo.tab*
  - __Alignment__ - folder containing script *mapToGenome_slurm_08072020.sh* to align patch-seq RNA samples to mouse reference genome  

- __data__ 

  - *counts_spinalCord.Rdata* contains raw GeneCounts and meta data
  - *normalized_counts.txt* contains normalized counts of GeneCounts using median of ratios method from DESeq2
  - *normalized_counts_geneIDs.txt* contains normalized counts of GeneCounts using median of ratios method from DESeq2 annotated by gene IDs

- *de_script.R* - R Studio script used to run DESeq2 (https://bioconductor.org/packages/release/workflows/html/RNAseq123.html)


- __figures__

  - folder contains plots created by running *de_script.R*
  

2. Relationship between files, if important:

  - *de_script.R* needs the following files in order to run smoothly
    - data/counts_spinalCord.Rdata
    - preprocessing/MouseIndexes/geneInfo.tab

 

METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 
Single-cell Ascl1, Ngn2 reprogrammed and not reprogrammed astrocytes were collected according to the patch-seq protocol published by Cadwell, et al. 2017 (https://doi.org/10.1038/nprot.2017.120). After the collection of single-cell RNA samples, we followed the protocol of SmartSeq2 for full-length RNA-seq published by Picelli, et al. 2014 (https://doi.org/10.1038/nprot.2014.006). Next generation sequencing was done using Illumina system at the collaborating sequencing facility Lafuga (https://www.genzentrum.uni-muenchen.de/research-groups/lafuga/index.html).

2. Methods for processing the data: 
Raw fastq files were demultiplexed and aligned using a high performance computing system available at the Bioinformatics Core Facility of BMC (https://github.com/bmc-CompBio/HPC_doc). Demultiplexing of the dataset was done with the help of Je_2.0.2 (https://github.com/gbcs-embl/Je/). The mouse reference genome and the alignment of the dataset was done with the help of STAR 2.7.1 (https://github.com/alexdobin/STAR).
