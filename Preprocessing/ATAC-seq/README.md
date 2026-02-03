# ATAC-seq processing

Here, you can find a detailed description of how ATAC-seq data were processed in this project.

## Dependencies

* [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://www.htslib.org/download/)
* [sambamba](https://github.com/biod/sambamba)
* [deepTools](https://deeptools.readthedocs.io/en/latest/)
* [HMMRATAC](https://macs3-project.github.io/MACS/docs/hmmratac.html)


## Summary of the workflow

1. **Trimming and Quality**: Trim Galore and FastQC
2. **Alignment**: Bowtie2
3. **Filtering**: samtools and sambamba
4. **Coverage**: deepTools
5. **Peak Calling**: HMMRATAC


The steps 1-4 were performed for the CUT&RUN data. A detailed description of these steps can be found [here](https://github.com/JavierreLab/B-Cell_Roadmap/tree/main/Preprocessing/CUT%26RUN).

## 5. Peak Calling


