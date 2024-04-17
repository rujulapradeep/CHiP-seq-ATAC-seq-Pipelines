# CHiP-seq/ATAC-seq Pipelines

A collection of R, Python, and command line pipelines designed for analyzing CHiP-seq and ATAC-seq data. These pipelines facilitate the identification of group specific peaks, annotation and motif analysis of these peaks, and result visualization. 

# Pre-requirements
- R
- Python
- Conda/Pip
- Homer

# Instruction
There are 6 pipelines currently. 

Here is the order to run them:
- peaks.sh
- heatmaps.R
- annotations.R
- annotations_piecharts.py
- motifs.sh
- motifs_visualization.R

In order to execute these pipelines, you will need bed files containing peaks for each sample of CHiP-seq/ATAC-seq data representing two groups (e.g., female and male). Additionally, you will need list of hg38 chromosome sizes and hg38 genome fasta fai file. 

## Sample Data Files

### Bed file containing peaks

```
chr1	191330	191677	Peak_24915	74	.	4.06537	7.46282	4.67336	141
chr1	844331	844643	Peak_24977	74	.	4.06537	7.46282	4.67336	168
chr1	870013	870503	Peak_24980	74	.	4.06537	7.46282	4.67336	142
chr1	905158	905816	Peak_24981	74	.	4.06537	7.46282	4.67336	422
chr1	906459	907081	Peak_32333	67	.	3.87232	6.76931	4.10857	236
```
#### *Important*
- no header
- seperated by tab

### Chromosome Sizes List

```
chr1 ####
chr2 ####
chr3 ####
chr4 ####
chr5 ####
```
#### *Important*
- no header
- seperated by space

### Genome Fasta Fai File
```
chr1	248956422	112	70	71
chr2	242193529	252513167	70	71
chr3	198295559	498166716	70	71
chr4	190214555	699295181	70	71
chr5	181538259	892227221	70	71
chr6	170805979	1076358996	70	71
chr7	159345973	1249605173	70	71
chr8	145138636	1411227630	70	71
chr9	138394717	1558439788	70	71
chr10	133797422	1698811686	70	71
```
#### *Important*
- no header
- seperated by space

# peaks.sh
```
peaks.sh <path_to_input_file> <path_to_output_directory>
```
## Input File Format
1) Sample Type: Specify types of samples (e.g., histone, ATAC-seq)
2) Group 1: Paths to bed files for each sample in group 1
3) Group 2: Paths to bed files for each sample in group 2
4) Repeat for Each Sample Type Repeat the above structure for additional sample types, if any    

### Sample Input File
```
H3K4me1
Data/H3K4me1/Female/female_46.bed	Data/H3K4me1/Female/female_51.bed	Data/H3K4me1/Female/female_53.bed	Data/H3K4me1/Female/female_56.bed	Data/H3K4me1/Female/female_59.bed
Data/H3K4me1/Male/male_40.bed	Data/H3K4me1/Male/male_43.bed	Data/H3K4me1/Male/male_54.bed	Data/H3K4me1/Male/male_55.bed	Data/H3K4me1/Male/male_61.bed	Data/H3K4me1/Male/male_66.bed	Data/H3K4me1/Male/male_69.bed	Data/H3K4me1/Male/male_73.bed
H3K4me3
Data/H3K4me3/Female/female_51.bed	Data/H3K4me3/Female/female_53.bed	Data/H3K4me3/Female/female_56.bed	Data/H3K4me3/Female/female_59.bed
Data/H3K4me3/Male/male_43.bed	Data/H3K4me3/Male/male_54.bed	Data/H3K4me3/Male/male_55.bed	Data/H3K4me3/Male/male_61.bed	Data/H3K4me3/Male/male_66.bed	Data/H3K4me3/Male/male_69.bed	Data/H3K4me3/Male/male_73.bed
CTCF
Data/CTCF/Female/female_46.bed	Data/CTCF/Female/female_51.bed	Data/CTCF/Female/female_56.bed	Data/CTCF/Female/female_59.bed
Data/CTCF/Male/male_40.bed	Data/CTCF/Male/male_43.bed	Data/CTCF/Male/male_54.bed	Data/CTCF/Male/male_55.bed	Data/CTCF/Male/male_61.bed	Data/CTCF/Male/male_66.bed	Data/CTCF/Male/male_69.bed	Data/CTCF/Male/male_73.bed
ATACseq
Data/ATACseq/Female/female_41.bed	Data/ATACseq/Female/female_42.bed	Data/ATACseq/Female/female_47.bed	Data/ATACseq/Female/female_51.bed	Data/ATACseq/Female/female_53.bed	Data/ATACseq/Female/female_59.bed	Data/ATACseq/Female/female_66.bed
Data/ATACseq/Male/male_43.bed	Data/ATACseq/Male/male_54.bed	Data/ATACseq/Male/male_55.bed	Data/ATACseq/Male/male_61.bed	Data/ATACseq/Male/male_66.bed	Data/ATACseq/Male/male_69.bed	Data/ATACseq/Male/male_73.bed

```
#### *Important*
- no header
- pathways seperated by tab
- 3 seperate lines for each sample type
- add empty line at end

## Output
***Upon execution, the pipeline will generate a <output_directory> folder containing subfolders for each sample type. Within each sample type folder, there will be group specific peak bed files for group 1 and group 2. These bed files will be used to execute the remaining analyses.***

# heatmaps.R
```
heatmaps.R <path_to_input_file> <path_to_output_directory> <path_to_fasta.fai_genome_file> <path_to_chromosome_sizes_file>
```
## Input File Format
1) path_to_group1_specific_peaks_bed_file / path_to_group2_specific_peaks_bed_file / sample_type
2) Repeat for Each Sample Type Repeat the above structure for additional sample types, if any

### Sample Input File
```
Peaks/H3K4me1/female_specific_peaks.bed	Peaks/H3K4me1/male_specific_peaks.bed	H3K4me1
Peaks/H3K4me3/female_specific_peaks.bed	Peaks/H3K4me3/male_specific_peaks.bed	H3K4me3
Peaks/H3K9me3/female_specific_peaks.bed	Peaks/H3K9me3/male_specific_peaks.bed	H3K9me3
Peaks/H3K27ac/female_specific_peaks.bed	Peaks/H3K27ac/male_specific_peaks.bed	H3K27ac
Peaks/H3K27me3/female_specific_peaks.bed Peaks/H3K27me3/male_specific_peaks.bed	H3K27me3
Peaks/H3K36me3/female_specific_peaks.bed Peaks/H3K36me3/male_specific_peaks.bed	H3K36me3
Peaks/CTCF/female_specific_peaks.bed Peaks/CTCF/male_specific_peaks.bed	CTCF
Peaks/ATACseq/female_specific_peaks.bed	Peaks/ATACseq/male_specific_peaks.bed	ATACseq
```

#### *Important*
- no header
- seperated by tab

## Output
***Upon execution, the pipeline will generate a <output_directory> folder containing subfolders for each sample type. Within each sample type folder, there will be group_1_heatmap.png and group_2_heatmap.png which are enrichment heatmaps visualizaing group1 and group2 specific peaks. Right now, group 1 output file will be labeled female_heatmap.png and group 2 output will be labeled male_heatmap.png; this will be edited later.***

### Output Examples
![image](https://github.com/rujulapradeep/CHiP-seq-ATAC-seq-Pipelines/assets/132700660/2228a347-96bf-4aa1-8da0-e04d544997f5) ![image](https://github.com/rujulapradeep/CHiP-seq-ATAC-seq-Pipelines/assets/132700660/4088d9ac-3bdf-4a11-83c5-3549e8eb0c40)

# annotations.R
```
annotations.R <path_to_input_file> <path_to_output_directory>
```
## Input File Format
1) path_to_group1_specific_peaks_bed_file / path_to_group2_specific_peaks_bed_file / sample_type
2) Repeat for Each Sample Type Repeat the above structure for additional sample types, if any

### Sample Input File
```
Peaks/H3K4me1/female_specific_peaks.bed	Peaks/H3K4me1/male_specific_peaks.bed	H3K4me1
Peaks/H3K4me3/female_specific_peaks.bed	Peaks/H3K4me3/male_specific_peaks.bed	H3K4me3
Peaks/H3K9me3/female_specific_peaks.bed	Peaks/H3K9me3/male_specific_peaks.bed	H3K9me3
Peaks/H3K27ac/female_specific_peaks.bed	Peaks/H3K27ac/male_specific_peaks.bed	H3K27ac
Peaks/H3K27me3/female_specific_peaks.bed Peaks/H3K27me3/male_specific_peaks.bed	H3K27me3
Peaks/H3K36me3/female_specific_peaks.bed Peaks/H3K36me3/male_specific_peaks.bed	H3K36me3
Peaks/CTCF/female_specific_peaks.bed Peaks/CTCF/male_specific_peaks.bed	CTCF
Peaks/ATACseq/female_specific_peaks.bed	Peaks/ATACseq/male_specific_peaks.bed	ATACseq
```

#### *Important*
- no header
- seperated by tab

## Output
***Upon execution, the pipeline will generate a <output_directory> folder containing subfolders for each sample type. Within each sample type folder, there will sampletype_group1_annotations.txt and sampletype_group2_annotations.txt which will have the annotated group specific peaks for each group. Right now, group 1 output file will be labeled sampletype_female_annotations.txt and group 2 output will be labeledsampletype_male_annotations.txt; this will be edited later. These output files will be utilized when executing the annotations_piechart pipeline.***

### Output Example
```
seqnames	start	end	width	strand	name	score	thick.start	thick.end	thick.width	annotation	geneChr	geneStart	geneEnd	geneLength	geneStrand	geneId	distanceToTSS	ENSEMBL	SYMBOL	GENENAME
1	chr1	1058878	1059813	936	*	Peak_206711	64	483	1059813	1059331	Promoter (<=1kb)	1	1059708	1069355	9648	1	100288175	0	ENSG00000217801	LOC100288175	uncharacterized LOC100288175
2	chr1	1058878	1059813	936	*	Peak_24331	743	807	1059813	1059007	Promoter (<=1kb)	1	1059708	1069355	9648	1	100288175	0	ENSG00000217801	LOC100288175	uncharacterized LOC100288175
3	chr1	1349601	1349986	386	*	Peak_222312	58	134	1349986	1349853	Promoter (<=1kb)	1	1335276	1349418	14143	2	1855	-183	ENSG00000107404	DVL1	dishevelled segment polarity protein 1
4	chr1	1349601	1349986	386	*	Peak_219150	59	323	1349986	1349664	Promoter (<=1kb)	1	1335276	1349418	14143	2	1855	-183	ENSG00000107404	DVL1	dishevelled segment polarity protein 1
.....
```

# annotations_piechart.py
```
Python annotations_piechart.py <path_to_input_file> <path_to_output_directory>
```
## Input File Format
1) path_to_group1_specific_peaks_annotations_file / path_to_group2_specific_peaks_annotations_file / sample_type
2) Repeat for Each Sample Type Repeat the above structure for additional sample types, if any

### Sample Input File

```
Annotations/H3K4me1_annotations/H3K4me1_female_annotations.txt	Annotations/H3K4me1_annotations/H3K4me1_male_annotations.txt	H3K4me1
Annotations/H3K4me3_annotations/H3K4me3_female_annotations.txt	Annotations/H3K4me3_annotations/H3K4me3_male_annotations.txt	H3K4me3
Annotations/H3K36me3_annotations/H3K36me3_female_annotations.txt	Annotations/H3K36me3_annotations/H3K36me3_male_annotations.txt	H3K36me3
Annotations/H3K9me3_annotations/H3K9me3_female_annotations.txt	Annotations/H3K9me3_annotations/H3K9me3_male_annotations.txt	H3K9me3
Annotations/H3K27ac_annotations/H3K27ac_female_annotations.txt	Annotations/H3K27ac_annotations/H3K27ac_male_annotations.txt	H3K27ac
Annotations/H3K27me3_annotations/H3K27me3_female_annotations.txt	Annotations/H3K27me3_annotations/H3K27me3_male_annotations.txt	H3K27me3
Annotations/CTCF_annotations/CTCF_female_annotations.txt	Annotations/CTCF_annotations/CTCF_male_annotations.txt	CTCF
Annotations/ATACseq_annotations/ATACseq_female_annotations.txt	Annotations/ATACseq_annotations/ATACseq_male_annotations.txt	ATACseq
```
#### *Important*
- no header
- seperated by tab

## Output
***Upon execution, the pipeline will generate a <output_directory> folder containing subfolders for each sample type. Within each sample type folder, there will be sampletype_group1_piechart.png and sampletype_group2_piechart.png which will visualize the distribution of group1 and group2 specific peaks across promoters, intergenic regions, and gene bodies. Right now, group 1 output file will be labeled sampletype_female_piechart.png and group 2 output will be labeled sampletype_male_piehcart.png; this will be edited later.***

### Output Example
![image](https://github.com/rujulapradeep/CHiP-seq-ATAC-seq-Pipelines/assets/132700660/2fbcdc0d-c90c-4f79-8135-6931eb80dac7)

# motifs.sh
```
motifs.sh <path_to_input_file> <path_to_output_directory>
```
## Input File Format
1) path_to_group1_specific_peaks_bed_file / path_to_group2_specific_peaks_bed_file / sample_type
2) Repeat for Each Sample Type Repeat the above structure for additional sample types, if any

### Sample Input File
```
Peaks/H3K4me1/female_specific_peaks.bed   Peaks/H3K4me1/male_specific_peaks.bed   H3K4me1
Peaks/H3K4me3/female_specific_peaks.bed   Peaks/H3K4me3/male_specific_peaks.bed   H3K4me3
Peaks/H3K9me3/female_specific_peaks.bed   Peaks/H3K9me3/male_specific_peaks.bed   H3K9me3
Peaks/H3K27ac/female_specific_peaks.bed  Peaks/H3K27ac/male_specific_peaks.bed  H3K27ac
Peaks/H3K27me3/female_specific_peaks.bed Peaks/H3K27me3/male_specific_peaks.bed H3K27me3
Peaks/H3K36me3/female_specific_peaks.bed Peaks/H3K36me3/male_specific_peaks.bed H3K36me3
Peaks/CTCF/female_specific_peaks.bed     Peaks/CTCF/male_specific_peaks.bed     CTCF
Peaks/ATACseq/female_specific_peaks.bed Peaks/ATACseq/male_specific_peaks.bed ATACseq

```
#### *Important*
- no header
- seperated by space
- add empty line at end

## Output
***Upon execution, the pipeline will generate a <output_directory> folder containing subfolders for each sample type. Within each sample type folder, there will be a homer_results_group1 folder and homer_results_group2 folder. Within each of these folders, there will be a knownResults.html and knownResults.txt file with motifs for group1 and group2 specific peaks. Right now, group 1 outputs will be labeled female and group 2 outputs will be labeled male; this will be edited later. The knownResults.txt for each group will be utilized to execute motifs_visualization.R.***

### Output Examples

```
REST-NRSF(Zf)/Jurkat-NRSF-ChIP-Seq/Homer	GGMGCTGTCCATGGTGCTGA	1e-13	-3.065e+01	0.0000	138.0	1.87%	381.5	0.92%
OCT4-SOX2-TCF-NANOG(POU,Homeobox,HMG)/mES-Oct4-ChIP-Seq(GSE11431)/Homer	ATTTGCATAACAATG	1e-9	-2.182e+01	0.0000	509.0	6.89%	2156.6	5.22%
Dlx3(Homeobox)/Kerainocytes-Dlx3-ChIP-Seq(GSE89884)/Homer	NDGTAATTAC	1e-7	-1.811e+01	0.0000	1745.0	23.63%	8661.9	20.95%
TCFL2(HMG)/K562-TCF7L2-ChIP-Seq(GSE29196)/Homer	ACWTCAAAGG	1e-7	-1.682e+01	0.0000	277.0	3.75%	1110.3	2.69%
E2F(E2F)/Hela-CellCycle-Expression/Homer	TTSGCGCGAAAA	1e-6	-1.390e+01	0.0001	599.0	8.11%	2760.1	6.68%
.....
```
![image](https://github.com/rujulapradeep/CHiP-seq-ATAC-seq-Pipelines/assets/132700660/58a5c093-744b-4b7d-8167-08766257bd6b)

# motifs_visualization.R
```
motifs_visualization.R <path_to_input_file> <path_to_output_directory>
```
## Input File Format
1) sample_type / path_to_group1_specific_peaks_motifs_txt_file / path_to_group2_specific_peaks_motifs_txt_file
2) Repeat for Each Sample Type Repeat the above structure for additional sample types, if any

### Sample Input File
```
H3K4me1 Motifs/H3K4me1_motifs/homer_results_female/knownResults.txt Motifs/H3K4me1_motifs/homer_results_male/knownResults.txt
H3K4me3 Motifs/H3K4me3_motifs/homer_results_female/knownResults.txt Motifs/H3K4me3_motifs/homer_results_male/knownResults.txt
H3K9me3 Motifs/H3K9me3_motifs/homer_results_female/knownResults.txt Motifs/H3K9me3_motifs/homer_results_male/knownResults.txt
H3K27ac Motifs/H3K27ac_motifs/homer_results_female/knownResults.txt Motifs/H3K27ac_motifs/homer_results_male/knownResults.txt
H3K27me3 Motifs/H3K27me3_motifs/homer_results_female/knownResults.txt Motifs/H3K27me3_motifs/homer_results_male/knownResults.txt
H3K36me3 Motifs/H3K36me3_motifs/homer_results_female/knownResults.txt Motifs/H3K36me3_motifs/homer_results_male/knownResults.txt
CTCF Motifs/CTCF_motifs/homer_results_female/knownResults.txt Motifs/CTCF_motifs/homer_results_male/knownResults.txt
ATACseq Motifs/ATACseq_motifs/homer_results_female/knownResults.txt Motifs/ATACseq_motifs/homer_results_male/knownResults.txt
```
#### *Important*
- no header
- seperated by space

## Output
***Upon execution, the pipeline will generate a <output_directory> folder containing subfolders for each sample type. Within each sample type folder, there will be a motifs_table.html which will visualize group1 and group2 motifs with their relative p-values.***

### Output Example

![image](https://github.com/rujulapradeep/CHiP-seq-ATAC-seq-Pipelines/assets/132700660/89e6ff19-7bfa-447e-a24d-d03827c98523)


  

 




                                                                


