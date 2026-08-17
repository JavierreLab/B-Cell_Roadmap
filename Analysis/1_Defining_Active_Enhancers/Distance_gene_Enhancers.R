library(tidyverse)
library(GenomicRanges)

# Distance between Active Enhancers and Target gene or nearest gene

### Procedure:
### Compute the distance to the nearest gene and the
### target gene for each Active Enhancer in each
### cell type.

# 1. Data ----
## Locate data necessary for this analysis

## promoters ----
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/promoter_regions_Homo_sapiens.GRCh38.104.bed")

#2. Reading Data ----
## Promoters ----
promoter <- data.table::fread("promoter_regions_Homo_sapiens.GRCh38.104.bed") |>
    mutate(TSS = paste(seqnames, TSS, sep = ":")) |>
    dplyr::rename(Gene_ID = ID, Gene_Name = Name, Gene_type = type, Gene_biotype = biotype) |>
    dplyr::select(Gene_ID, Gene_Name, TSS, Gene_type, Gene_biotype)


# Gene - Active Enhancers pairs info ----
gene_enh <- data.table::fread(file = "Gene_ActiveEnhancer_pairs_per_celltype.tsv", na.strings = "") ## Generated in Defining_Enhancers_ATAC_liCHiC.R
gene_enh[gene_enh == ""] <- NA

#3. Computing distances ----
## To nearest gene ----
enhGR <- gene_enh %>%
    filter(!is.na(Active_Enhancer)) %>%
    dplyr::select(Active_Enhancer, cell) %>%
    unique() %>%
    separate(Active_Enhancer, into = c("seqnames", "ranges"), sep = ":") %>%
    separate(ranges, into = c("start", "end"), sep = "-") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

TSSGR <- promoter %>%
    dplyr::select(Gene_ID, TSS, Gene_biotype) %>%
    separate(TSS, into = c("seqnames", "start"), sep = ":") %>%
    mutate(end = as.numeric(start) + 1) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

near <- resize(TSSGR[nearest(enhGR, TSSGR)], width = 1)

ginear <- GenomicInteractions::GenomicInteractions(enhGR, near)
ginear$distance <- GenomicInteractions::calculateDistances(ginear)

ginear$group <- "Nearest Gene"

## To target gene ----
gene_enhGR <- gene_enh %>%
    filter(!is.na(Active_Enhancer) & Gene_ID != "non-annotated") %>%
    dplyr::select(Active_Enhancer, cell) %>%
    separate(Active_Enhancer, into = c("seqnames", "ranges"), sep = ":") %>%
    separate(ranges, into = c("start", "end"), sep = "-") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

enh_geneGR <- gene_enh %>%
    filter(!is.na(Active_Enhancer) & Gene_ID != "non-annotated") %>%
    dplyr::select(TSS, Gene_ID, Gene_biotype) %>%
    separate(TSS, into = c("seqnames", "start"), sep = ":") %>%
    mutate(end = start) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)

gitarget <- GenomicInteractions::GenomicInteractions(gene_enhGR, enh_geneGR)
gitarget$distance <- GenomicInteractions::calculateDistances(gitarget)

gitarget$group <- "Target Gene"


#4. Plot distances ----
toplot <- bind_rows(as_tibble(ginear)[, c("distance", "group", "anchor1.cell", "anchor2.Gene_ID", "anchor2.Gene_biotype")], as_tibble(gitarget)[, c("distance", "group", "anchor1.cell", "anchor2.Gene_ID", "anchor2.Gene_biotype")]) %>%
    dplyr::rename(cell = anchor1.cell, Gene_ID = anchor2.Gene_ID, Gene_biotype = anchor2.Gene_biotype) %>%
    mutate(cell = factor(cell, levels = c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC", "CMP", "Mon", "nCD8")))

toplot %>%
    group_by(group,cell) %>%
    summarise(med = median(distance, na.rm = TRUE)) |> 
    group_by(group) |> 
    reframe(mean = log10(mean(med)))

toplot %>%
    group_by(group,cell) %>%
    summarise(med = median(distance, na.rm = TRUE)) |> 
    group_by(group) |> 
    reframe(mean = mean(med))

### Figure 1B
ggplot(toplot, aes(x = log10(distance), fill = group, col = group)) +
    geom_density() +
    ggpubr::theme_classic2() +
    scale_fill_manual(values = c("#bebebe80", "#c1153ab2")) +
    scale_color_manual(values = c("#8f8f8fff", "#c1153aff")) +
    geom_vline(
        data = toplot %>%
            group_by(group,cell) %>%
            summarise(med = median(distance, na.rm = TRUE)) |> 
            group_by(group) |> 
            reframe(mean = log10(mean(med))),
        aes(xintercept = mean), col = "grey", linetype = "dashed"
    ) +
    geom_text(
        data = toplot %>%
            group_by(group,cell) %>%
            summarise(med = median(distance, na.rm = TRUE)) |> 
            group_by(group) |> 
            reframe(mean = log10(mean(med)),
                    label = paste0(round(mean(med)), " bp")
            ) %>%
            mutate(
                mean = ifelse(group == "Nearest Gene", mean - 1.5, mean + 1.5),
                y = ifelse(group == "Nearest Gene", 0.8, 1)
            ),
        aes(x = mean, y = y, label = label), inherit.aes = FALSE
    ) + 
    guides(col = "none") +
    labs(
        fill = NULL, x = "log10(distance in bp)", y = "Density",
        subtitle = "Distance from AE to genes"
    )