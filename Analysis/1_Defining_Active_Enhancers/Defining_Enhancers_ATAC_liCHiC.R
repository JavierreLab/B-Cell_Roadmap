library(HiCaptuRe)
library(tidyverse)

# Functional Enhancers and their target genes

### Functional Enhancer definition:
### Enhancers (active, primed, poised, super) that 
### are open and interacts with a gene.

### Procedure:
### Integrate enhancers from ChromHMM with
### ATAC-seq peaks so only the open (potentially 
### functional) part is kept. Then use liCHi-C data to
### keep the ones interacting with genes (potentially 
### regulating). Finally pair Enhancers with target genes.

# 1. Downloading data ----
## Download all processed public data and 
## complementary files needed for the analysis

## promoters ----
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/promoter_regions_Homo_sapiens.GRCh38.104.bed")

## lichic ----
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320102/suppl/GSE320102_liCHiC_CHiCAGO_peakmatrix_merge.txt.gz")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/liCHi-C/baits_coordinates_annotated_ensembl_gene_id_GRCh38_104.bed")

## ATAC ----
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_CMP_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_GCB_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_HSC_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_Mon_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_PC_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_PreB_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_PreProB_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_ProB_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_immtransB_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_memB_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_nB_ATAC_merged_peaks.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320104/suppl/GSE320104_nCD8_ATAC_merged_peaks.bed.gz")

ATAC_files <- list.files(pattern = "ATAC_merged")

## ChromHMM ----
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_CMP_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_GCB_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_HSC_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_Mon_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_PC_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_PreB_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_PreProB_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_ProB_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_immtransB_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_memB_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_nB_ChromHMM.bed.gz")
system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE320nnn/GSE320103/suppl/GSE320103_nCD8_ChromHMM.bed.gz")

chromHMM_files <- list.files(pattern = "ChromHMM")

## SuperEnhancers ----
### Defined as explained in section 8 of the tutorial https://github.com/JavierreLab/B-Cell_Roadmap/blob/main/Preprocessing/CUT%26RUN/README.md
### Note: this super-Enhancers are already filtered by ATAC-seq
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/CMP_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/GCB_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/HSC_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/Mon_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/PC_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/PreB_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/PreProB_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/ProB_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/immtransB_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/memB_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/nB_Super_Enhancers.bed")
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/SuperEnhancers/nCD8_Super_Enhancers.bed")

SE_files <- list.files(pattern = "Super_Enhancers")

#2. Reading Data ----
## Promoters ----
promoter <- data.table::fread("promoter_regions_Homo_sapiens.GRCh38.104.bed") |>
    mutate(TSS = paste(seqnames, TSS, sep = ":")) |>
    dplyr::rename(Gene_ID = ID, Gene_Name = Name, Gene_type = type, Gene_biotype = biotype) |>
    dplyr::select(Gene_ID, Gene_Name, TSS, Gene_type, Gene_biotype)

## lichic ----
pm <- load_interactions(file = "GSE320102_liCHiC_CHiCAGO_peakmatrix_merge.txt.gz") 
### Annotate the peakmatrix with Ensembl gene ID
pm <- annotate_interactions(interactions = pm, annotation = "baits_coordinates_annotated_ensembl_gene_id_GRCh38_104.bed")

### Split peakmatrix by cell type
pm_l <- peakmatrix2list(pm)
### Keep only the interactions bait with other-end since 
### we'll focus on enhancers potentially regulating genes
pm_l <- lapply(pm_l, function(x) x[x$int == "B_OE"])


# 3. Gene-Enhancer pairs definition ----

gene_AE_pairs <- data.frame()
gene_PrE_pairs <- data.frame()
gene_PoE_pairs <- data.frame()
gene_SE_pairs <- data.frame()

for (cell in names(pm_l)) 
{
    print(cell)
    chromhmm_file <- grep(paste0("_",cell), chromHMM_files, value = T)
    ATAC_file <- grep(paste0("_",cell), ATAC_files, value = T)
    SE_file <- grep(paste0("^",cell), SE_files, value = T)
    
    ints <- pm_l[[cell]]
    
    chromhmmGR <- data.table::fread(file = chromhmm_file)[,1:4] |>
        makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3",keep.extra.columns = T) |>
        sort()
    colnames(mcols(chromhmmGR)) <- "States"

    AE_GR <- chromhmmGR |> plyranges::filter(States %in% 4:5) |> GenomicRanges::reduce() ## Active Enhancers - states 4 and 5
    PrE_GR <- chromhmmGR |> plyranges::filter(States %in% 6:7) |> GenomicRanges::reduce() ## Primed Enhancers - states 6 and 7
    PoE_GR <- chromhmmGR |> plyranges::filter(States %in% 8:9) |> GenomicRanges::reduce() ## Poised Enhancers - states 8 and 9

    SE_ATAC <- data.table::fread(file = SE_file)[,1:3] |> 
        makeGRangesFromDataFrame() |> 
        sort()  
    
    ATAC_GR <- data.table::fread(file = ATAC_file) |>
        makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")
    
    ### Defining open region of each enhancer
    AE_ATAC <- intersect(AE_GR, ATAC_GR)
    PrE_ATAC <- intersect(PrE_GR, ATAC_GR)
    PoE_ATAC <- intersect(PoE_GR, ATAC_GR)

    ### Keeping interactions with Enhancers
    ints_AE_ATAC <- interactionsByRegions(ints, regions = AE_ATAC)
    ints_PrE_ATAC <- interactionsByRegions(ints, regions = PrE_ATAC)
    ints_PoE_ATAC <- interactionsByRegions(ints, regions = PoE_ATAC)
    ints_SE_ATAC <- interactionsByRegions(ints, regions = SE_ATAC)
    
    ### Defining and saving open part of Enhancers that interact with at least 1 gene
    ints_AE_ATAC_enhancers <- as_tibble(getByRegions(ints_AE_ATAC)[[1]]) |>
        filter(!is.na(Nfragment) & NfragmentOE != 0) |>
        dplyr::select(seqnames, start, end)
    data.table::fwrite(ints_AE_ATAC_enhancers, paste0(cell, "_Active_Enhancers.bed"), col.names = T, row.names = F, quote = F, sep = "\t")

    ints_PrE_ATAC_enhancers <- as_tibble(getByRegions(ints_PrE_ATAC)[[1]]) |>
        filter(!is.na(Nfragment) & NfragmentOE != 0) |>
        dplyr::select(seqnames, start, end)
    data.table::fwrite(ints_PrE_ATAC_enhancers, paste0(cell, "_Primed_Enhancers.bed"), col.names = T, row.names = F, quote = F, sep = "\t")

    ints_PoE_ATAC_enhancers <- as_tibble(getByRegions(ints_PoE_ATAC)[[1]]) |>
        filter(!is.na(Nfragment) & NfragmentOE != 0) |>
        dplyr::select(seqnames, start, end)
    data.table::fwrite(ints_PoE_ATAC_enhancers, paste0(cell, "_Poised_Enhancers.bed"), col.names = T, row.names = F, quote = F, sep = "\t")
    
    ints_SE_ATAC_enhancers <- as_tibble(getByRegions(ints_SE_ATAC)[[1]]) |> 
        filter(!is.na(Nfragment) & NfragmentOE != 0) |> 
        dplyr::select(seqnames, start, end)
    data.table::fwrite(ints_SE_ATAC_enhancers, paste0(cell, "_Super_Enhancers.bed"), col.names = T, row.names = F, quote = F, sep = "\t")
    
    ints_AE_ATAC_2 <- ints_AE_ATAC[ints_AE_ATAC$region_2]
    ints_PrE_ATAC_2 <- ints_PrE_ATAC[ints_PrE_ATAC$region_2]
    ints_PoE_ATAC_2 <- ints_PoE_ATAC[ints_PoE_ATAC$region_2]
    ints_SE_ATAC_2 <- ints_SE_ATAC[ints_SE_ATAC$region_2]
    
    ### bait centric table with all OE interacting with them
    bait_centric_all_OE <- as_tibble(ints) %>%
        mutate(
            bait_region = paste(seqnames1, paste(start1, end1, sep = "-"), sep = ":"),
            OE_region = paste(seqnames2, paste(start2, end2, sep = "-"), sep = ":")
        ) %>%
        dplyr::select(ID_1, bait_1, bait_region, ID_2, OE_region) %>%
        separate_rows(bait_1, sep = ",") %>%
        dplyr::rename(ID = ID_1, bait = bait_1, interactingID = ID_2)
    
    ### Active Enhancer centric table with all the OE that they overlap with
    state_centric_OE_AE_ATAC <- as_tibble(getByRegions(ints_AE_ATAC_2)[[1]]) |>
        separate_rows(fragmentID) |>
        mutate(State_region = paste(seqnames, paste(start, end, sep = "-"), sep = ":"),
               fragmentID = as.numeric(fragmentID)) |>
        filter(!is.na(Nfragment) & NfragmentOE != 0) |>
        dplyr::select(State_region, fragmentID) |>
        filter(!is.na(fragmentID)) |>
        dplyr::rename(interactingID = fragmentID)

    ### merge both info so for all interactions keep the info whether it has an Active Enhancer or not and which one
    state_centric_genes_AE_ATAC <- left_join(x = bait_centric_all_OE, y = state_centric_OE_AE_ATAC,
                                             by = "interactingID", relationship = "many-to-many")


    ### Primed Enhancer centric table with all the OE that they overlap with
    state_centric_OE_PrE_ATAC <- as_tibble(getByRegions(ints_PrE_ATAC_2)[[1]]) |>
        separate_rows(fragmentID) |>
        mutate(State_region = paste(seqnames, paste(start, end, sep = "-"), sep = ":"),
               fragmentID = as.numeric(fragmentID)) |>
        dplyr::select(State_region, fragmentID) |>
        filter(!is.na(fragmentID)) |>
        dplyr::rename(interactingID = fragmentID)

    ### merge both info so for all interactions keep the info whether it has an Primed Enhancer or not and which one
    state_centric_genes_PrE_ATAC <- left_join(x = bait_centric_all_OE, y = state_centric_OE_PrE_ATAC,
                                             by = "interactingID", relationship = "many-to-many")


    ### Poised Enhancer centric table with all the OE that they overlap with
    state_centric_OE_PoE_ATAC <- as_tibble(getByRegions(ints_PoE_ATAC_2)[[1]]) |>
        separate_rows(fragmentID) |>
        mutate(State_region = paste(seqnames, paste(start, end, sep = "-"), sep = ":"),
               fragmentID = as.numeric(fragmentID)) |>
        dplyr::select(State_region, fragmentID) |>
        filter(!is.na(fragmentID)) |>
        dplyr::rename(interactingID = fragmentID)

    ### merge both info so for all interactions keep the info whether it has an Poised Enhancer or not and which one
    state_centric_genes_PoE_ATAC <- left_join(x = bait_centric_all_OE, y = state_centric_OE_PoE_ATAC,
                                             by = "interactingID", relationship = "many-to-many")

    
    ### Super Enhancer centric table with all the OE that they overlap with
    state_centric_OE_SE_ATAC <- as_tibble(getByRegions(ints_SE_ATAC_2)[[1]]) |>
        separate_rows(fragmentID) |>
        mutate(State_region = paste(seqnames, paste(start, end, sep = "-"), sep = ":"),
               fragmentID = as.numeric(fragmentID)) |>
        dplyr::select(State_region, fragmentID) |>
        filter(!is.na(fragmentID)) |>
        dplyr::rename(interactingID = fragmentID)
    
    ### merge both info so for all interactions keep the info whether it has an Super Enhancer or not and which one
    state_centric_genes_SE_ATAC <- left_join(x = bait_centric_all_OE, y = state_centric_OE_SE_ATAC,
                                              by = "interactingID", relationship = "many-to-many")
    
    ### Rename columns, keep cell info and merge with previous cell type
    colnames(state_centric_genes_AE_ATAC) <- c("bait_ID", "Gene_ID", "bait_region", "OE_ID", "OE_region", "Active_Enhancer")
    state_centric_genes_AE_ATAC$cell <- cell
    gene_AE_pairs <- rbind(gene_AE_pairs, state_centric_genes_AE_ATAC)

    colnames(state_centric_genes_PrE_ATAC) <- c("bait_ID", "Gene_ID", "bait_region", "OE_ID", "OE_region", "Primed_Enhancer")
    state_centric_genes_PrE_ATAC$cell <- cell
    gene_PrE_pairs <- rbind(gene_PrE_pairs, state_centric_genes_PrE_ATAC)

    colnames(state_centric_genes_PoE_ATAC) <- c("bait_ID", "Gene_ID", "bait_region", "OE_ID", "OE_region", "Poised_Enhancer")
    state_centric_genes_PoE_ATAC$cell <- cell
    gene_PoE_pairs <- rbind(gene_PoE_pairs, state_centric_genes_PoE_ATAC)
    
    colnames(state_centric_genes_SE_ATAC) <- c("bait_ID", "Gene_ID", "bait_region", "OE_ID", "OE_region", "Super_Enhancer")
    state_centric_genes_SE_ATAC$cell <- cell
    gene_SE_pairs <- rbind(gene_SE_pairs, state_centric_genes_SE_ATAC)
    
}

### Add extra info for promoters and save gene-enhancers pairs
gene_AE_pairs <- left_join(gene_AE_pairs, promoter, by = "Gene_ID") |> 
    select(Gene_ID, Gene_Name, TSS, Gene_type, Gene_biotype, bait_ID, bait_region, OE_ID, OE_region, Active_Enhancer, cell)
data.table::fwrite(gene_AE_pairs, file = "Gene_ActiveEnhancer_pairs_per_celltype.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

gene_PrE_pairs <- left_join(gene_PrE_pairs, promoter, by = "Gene_ID") |> 
    select(Gene_ID, Gene_Name, TSS, Gene_type, Gene_biotype, bait_ID, bait_region, OE_ID, OE_region, Primed_Enhancer, cell)
data.table::fwrite(gene_PrE_pairs, file = "Gene_PrimedEnhancer_pairs_per_celltype.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

gene_PoE_pairs <- left_join(gene_PoE_pairs, promoter, by = "Gene_ID") |> 
    select(Gene_ID, Gene_Name, TSS, Gene_type, Gene_biotype, bait_ID, bait_region, OE_ID, OE_region, Poised_Enhancer, cell)
data.table::fwrite(gene_PoE_pairs, file = "Gene_PoisedEnhancer_pairs_per_celltype.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

gene_SE_pairs <- left_join(gene_SE_pairs, promoter, by = "Gene_ID") |> 
    select(Gene_ID, Gene_Name, TSS, Gene_type, Gene_biotype, bait_ID, bait_region, OE_ID, OE_region, Super_Enhancer, cell)
data.table::fwrite(gene_SE_pairs, file = "Gene_SuperEnhancer_pairs_per_celltype.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
