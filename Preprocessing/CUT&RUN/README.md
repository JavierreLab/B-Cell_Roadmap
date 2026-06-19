# CUT&RUN processing

Here, you can find a detailed description of how CUT&RUN data were processed in this project.

## Dependencies

* [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://www.htslib.org/download/)
* [sambamba](https://github.com/biod/sambamba)
* [deepTools](https://deeptools.readthedocs.io/en/latest/)
* [featureCounts](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html)
* [csaw (R package)](https://bioconductor.org/packages/release/bioc/html/csaw.html)
* [GenomicRanges (R package)](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
* [DESeq2 (R package)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* [ChromHMM](https://ernstlab.github.io/ChromHMM/)
* [ChromTime](https://github.com/ernstlab/ChromTime)
* [ROSE](http://younglab.wi.mit.edu/super_enhancer_code.html)



## Summary of the workflow

1. **Trimming and Quality**: Trim Galore and FastQC
2. **Alignment**: Bowtie2* (with CUT&RUN-specific alignment parameters)
3. **Filtering**: samtools and sambamba
4. **Coverage**: deepTools
5. **Promoter Counts**: csaw and GenomicRanges
6. **Background Matrix**: csaw and GenomicRanges
7. **Chromatin Segmentation**: ChromHMM and ChromTime
8. **Super Enhancer Identification**: ROSE

A detailed description of how steps 1-4 were performed can be found [here](https://github.com/JavierreLab/liCHiC/tree/main/3.ChIPseq%20Processing), the GitHub page associated with Tomás-Daza, L *et al.* Low input capture Hi-C (liCHi-C) identifies promoter-enhancer interactions at high-resolution. *Nature Communications* **14**, 268 (2023). [https://doi.org/10.1038/s41467-023-35911-8](https://doi.org/10.1038/s41467-023-35911-8)
We follow mostly the same logic as for a ChIP-seq experiment.

## 2. Alignment*

The alignment step is the main point of divergence from the ChIP-seq preprocessing pipeline. In particular, we enable the --dovetail option in Bowtie2 to retain paired-end reads in which mates overlap. Such overlapping pairs are expected in CUT&RUN data due to the generation of very short DNA fragments and would otherwise be discarded by default Bowtie2 filtering. Retaining these reads is therefore essential to avoid the systematic loss of valid CUT&RUN signal.

```bash
bowtie2 -x $INDEX -1 $FASTQ1 -2 $FASTQ2 --very-sensitive-local --no-unal --no-mixed --no-discordant -k 2 --phred33 -I 10 -X 700 --dovetail -p $THREADS -S $OUTDIR/$NAME.sam
```
## 5. Promoter Counts

This step and the following are performed to quantify the signal for a given histone mark (specifically H3K27ac and H3K27me3) across different regions of the genome, in this case at the promoter level.

First, we need to slightly modify the [promoter file](https://github.com/JavierreLab/B-Cell_Roadmap/blob/main/Data/Complementary_Data/promoter_regions_Homo_sapiens.GRCh38.104.bed) to fit the requirements of the featureCounts. 

```R
promoters <- read.table("promoter_regions_Homo_sapiens.GRCh38.104.bed",header = T)
### Change to SAF format to use with FeatureCount from subread
saf_prom <- promoters[,c(6,1:4)]
colnames(saf_prom) <- c("GeneID","Chr","Start","End","Strand")
saf_prom$Start[saf_prom$Start < 0] <- 1
saf_prom <- saf_prom[saf_prom$Chr %in% unique(saf_prom$Chr)[c(1:22,45:47)],]
write.table(saf_prom, file="promoter_saf.bed",col.names = T, row.names = F, quote = F, sep = "\t")
```
Once we have generated the promoter file, we simply need to run featureCounts for these regions in all the bams files for the given histone mark.

```bash
featureCounts -a promoter_saf.bed -F SAF -O -p -d 0 -D 10000 -T 30 -o CUTnRUN_H3K27ac_promoter_counts.tsv $BAMS
```

## 6. Background Matrix

A detailed description of how this step was performed can be found [here](https://github.com/JavierreLab/p53/tree/main/preprocessing/ChIP#6-generating-peakmatrix), the github page associated with Serra, F *et al.* p53 rapidly restructures 3D chromatin organization to trigger a transcriptional response. *Nature Communications* **15**, 2821 (2024). [https://doi.org/10.1038/s41467-024-46666-1](https://doi.org/10.1038/s41467-024-46666-1)

We did slight modifications; we used the ChromHMM binarization as a proxy of peaks to define the background. The explanation on how to obtain the bed from the binarization is in section [7.3.3.](https://github.com/JavierreLab/B-Cell_Roadmap/tree/main/Preprocessing/CUT%26RUN#733-generating-peaks-bed-files-from-binarization)

## 7. Chromatin Segmentation

### Overview

This section describes the complete chromatin segmentation workflow performed after CUT&RUN preprocessing, starting from filtered BAM files obtained after alignment and filtering steps.

---

### 7.1. Starting input files

We started from **filtered BAM files** generated after alignment and filtering steps.

Replicates strategy:

- **ChromTime (narrow marks)** → replicates were supermerged per cell type.
- **ChromHMM broad mark binarization** → all replicates were used individually.
- IgG was used as background control in all relevant steps.

---

### 7.2. ChromTime Peak Calling (Narrow Histone Marks)

For narrow histone marks **(H3K27ac and H3K4me3)**, we performed peak calling using ChromTime, a trajectory-aware peak caller designed to model temporal dynamics.

#### 7.2.1 Strategy

- **Supermerged BAM files** (per cell type) were used.
- Background: supermerged IgG BAM per cell type.
- Two differentiation trajectories were modeled:
  - Trajectory from HSC cell type to **memB** cell type
  - Trajectory from HSC cell type to **PC** cell type

After running both trajectories:

1. Peaks were extracted per cell type.
2. The **intersection of peaks** from both trajectories was computed.
3. The resulting intersected BED files were used for ChromHMM binarization.

#### ChromTime Command

```bash

python ChromTime.py -i $Path_to_ChromTime_input_table -o $Path_to_ChromTime_output -t 8 -g grc38 -p ct -q 0.01 --q-value-seed 0.01 --p-value-extend 0.01

```

---

### 7.3. ChromHMM Binarization (200 bp resolution)

We performed binarization using ChromHMM at the default 200 bp resolution, following the official manual: https://ernstlab.github.io/ChromHMM/ChromHMM_manual.pdf

#### 7.3.1. Broad Histone Marks

For broad marks **(H3K4me1, H3K27me3, H3K9me3)** we used `BinarizeBam` directly on filtered BAM files:

```
java -mx4000M -jar ChromHMM.jar BinarizeBam \
  $chr_length_file \
  $samples_directory \
  $cell_mark_file_broad \
  $output_binary_directory

```

#### 7.3.2. Narrow Histone Marks

For narrow marks **(H3K27ac and H3K4me3)**, we used `BinarizeBed -peaks` using **peak files** generated by ChromTime:

```
java -mx4000M -jar ChromHMM.jar BinarizeBed -peaks \
  $chr_length_file \
  $samples_directory_intersect_peaks_ChromTime \
  $cell_mark_file_narrow \
  $output_binary_directory

```

#### 7.3.3. Generating Peaks BED Files from Binarization

After binarization, we generated for each histone mark and cell type a corresponding **peaks BED file** for downstream analyses. This step was done in order to obtain the genomic coordinates of each histone mark as a proxy of "peaks" generated by a peak caller.

Briefly, we identified all bins marked as present (`1`) for a given histone mark and reconstructed their genomic coordinates based on the fixed ChromHMM bin size (200 bp).

For each positive bin:

- The start coordinate (0-based) was calculated as: start_position = (`bin_index` - 1) * 200
- The end coordinate (1-based) was calculated as: end_position = start_position + 200

Where:

- `bin_index` corresponds to the row number of the bin in the binarized file.


#### 7.3.4. Merging Narrow and Broad binarizations

For each chromosome and each cell type, we merged the binarization output files from Narrow and Broad histone marks into a unified dataset, which was used for model learning.

### 7.4. ChromHMM Model Learning

Using the merged binarized files, we trained a ChromHMM model with **15 states** using `LearnModel` function:

```
java -mx20000M -jar ChromHMM.jar LearnModel \
  -p 96 \
  $output_binary_directory \
  $output_model_directory \
  15 \
  GRCh38

```
Finally, the emission matrix obtained was annotated based on prior biological knowledge to relabel the 15 states into biologically meaningful chromatin state categories.

## 8. Super Enhancer Identification

### Overview

This section describes the workflow to define super enhancers, starting from H3K27ac and IgG filtered BAM files and active enhancers defined after ChromHMM processing. This pipeline was run in a conda environment using:
- python 2.7.18
- samtools 1.19.2
- R 4.3.2
  

### 8.1. Prepare input files

ROSE requires two types of input files:
- **BAM files**: we used **H3K27ac** and **IgG** merged for each cell type. Remember that ROSE requires that .bam files have chromosome IDs starting with "chr" and that they are sorted and indexed.
- **GFF of constituent enhancers**: regions defined as open active enhancers. In this case, we used active enhancers overlapping ATAC-seq peaks for each cell type. GFF file must be formatted as explained in the ROSE manual.
  

### 8.2. Run ROSE

We run ROSE for each cell type following this command:

```
python ROSE_main.py -g HG38 \
-i $Celltype_OpenActiveEnhancer.gff \
-r $Celltype_H3K27ac.bam \
-c $Celltype_IgG.bam \
-o $Celltype_output_directory \
-s 12500 \
-t 2500

```
