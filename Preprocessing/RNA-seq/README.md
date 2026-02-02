# RNA-seq processing

Here, you can find a detailed description of how RNA-seq data were processed in this project.

## Dependencies

* [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [STAR](https://github.com/alexdobin/STAR/tree/master)
* [sambamba](https://github.com/biod/sambamba)
* [featureCounts](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html)
* [DESeq2 (R package)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)


## Summary of the workflow

1. **Trimming and Quality**: Trim Galore, samtools and FastQC
2. **Alignment**: Bowtie2 
3. **Counting**: sambamba and featureCounts
4. **Differential Analysis**: DESeq2


A detailed description of how steps 1-4 were performed can be found [here](https://github.com/JavierreLab/p53/tree/main/preprocessing/RNA), the github page associated to Serra, F *et al.* p53 rapidly restructures 3D chromatin organization to trigger a transcriptional response. *Nature Communications* **15**, 2821 (2024). [https://doi.org/10.1038/s41467-024-46666-1](https://doi.org/10.1038/s41467-024-46666-1)
