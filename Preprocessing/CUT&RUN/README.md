# CUT&RUN processing

Here, you can find a detailed description of how CUT&RUN data were processed in this project.

## Dependencies

* [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://www.htslib.org/download/)
* [sambamba](https://github.com/biod/sambamba)
* [deepTools](https://deeptools.readthedocs.io/en/latest/)
* [macs2](https://github.com/macs3-project/MACS)
* [csaw (R package)](https://bioconductor.org/packages/release/bioc/html/csaw.html)
* [GenomicRanges (R package)](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
* [DESeq2 (R package)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)


## Summary of the workflow

1. **Trimming and Quality**: Trim Galore and FastQC
2. **Alignment**: Bowtie2*
3. **Filtering**: samtools and sambamba
4. **Coverage**: deepTools
5. **Generating Peakmatrix**: csaw and GenomicRanges
6. **Generating Background Matrix**: csaw and GenomicRanges
* with CUT&RUN-specific alignment parameters

A detailed description of how steps 1-4 were performed can be found [here](https://github.com/JavierreLab/liCHiC/tree/main/3.ChIPseq%20Processing), the GitHub page associated to Tomás-Daza, L *et al.* Low input capture Hi-C (liCHi-C) identifies promoter-enhancer interactions at high-resolution. *Nature Communications* **14**, 268 (2023). [https://doi.org/10.1038/s41467-023-35911-8](https://doi.org/10.1038/s41467-023-35911-8)
We follow mostly the same logic as for a ChIP-seq experiment.

## 2. Alignment*

The alignment step is the main point of divergence from the ChIP-seq preprocessing pipeline. In particular, we enable the --dovetail option in Bowtie2 to retain paired-end reads in which mates overlap. Such overlapping pairs are expected in CUT&RUN data due to the generation of very short DNA fragments and would otherwise be discarded by default Bowtie2 filtering. Retaining these reads is therefore essential to avoid the systematic loss of valid CUT&RUN signal.

```bash
bowtie2 -x $INDEX -1 $FASTQ1 -2 $FASTQ2 --very-sensitive-local --no-unal --no-mixed --no-discordant -k 2 --phred33 -I 10 -X 700 --dovetail -p $THREADS -S $OUTDIR/$NAME.sam
```
