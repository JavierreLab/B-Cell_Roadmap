# PCHi-C processing

Here you can find a detailed description of how the liCHi-C data were processed in this project.

## Dependencies

* [HiCUP](https://www.bioinformatics.babraham.ac.uk/projects/hicup/)
* [Chicago (R package)](https://www.bioconductor.org/packages/release/bioc/html/Chicago.html)
* [bowtie2](https://github.com/BenLangmead/bowtie2)


## Summary of the workflow

1. **Mapping and Filtering**: HiCUP
2. **Capture Efficiency**: HiCUP miscellaneous script
3. **Interaction Calling**: Chicago


A detailed description of how steps 1-3 were performed can be found [here](https://github.com/JavierreLab/liCHiC/tree/main/1.liCHiC%20Processing), the github page associated to Tomás-Daza, L *et al.* Low input capture Hi-C (liCHi-C) identifies promoter-enhancer interactions at high-resolution. *Nature Communications* **14**, 268 (2023). [https://doi.org/10.1038/s41467-023-35911-8](https://doi.org/10.1038/s41467-023-35911-8)
