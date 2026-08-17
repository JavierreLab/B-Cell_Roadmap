library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(UpSetR)

# Categorization of Active Enhancer based on genomic position and targets

### Procedure:
### Compute the proportion of Active Enhancers in intragenic and 
### intergenic regions. And for intragenic ones split them if they
### interact with the same gene, a different one or both.

# 1. Data ----
## Locate data necessary for this analysis

## promoters ----
system("wget https://raw.githubusercontent.com/JavierreLab/B-Cell_Roadmap/refs/heads/main/Data/Complementary_Data/promoter_regions_Homo_sapiens.GRCh38.104.bed")

promoters <- data.table::fread("promoter_regions_Homo_sapiens.GRCh38.104.bed") |> 
    filter(biotype == "protein_coding")  

## Genome and gene annotation ----
system("wget https://ftp.ensembl.org/pub/release-104/gff3/homo_sapiens/Homo_sapiens.GRCh38.104.gff3.gz")

txdb <-  txdbmaker::makeTxDbFromGFF("Homo_sapiens.GRCh38.104.gff3.gz", format = "gff3",organism = "Homo sapiens",)
genome(txdb) <- "hg38"
seqlevels(txdb) <- paste0("chr",seqlevels(txdb))

tx2gene <- GenomicFeatures::transcripts(txdb, columns = c("tx_id", "gene_id"))

## Gene - Active Enhancers pairs ----
gene_enh <- data.table::fread(file = "Gene_ActiveEnhancer_pairs_per_celltype.tsv", na.strings = "")
gene_enh[gene_enh == ""] <- NA
gene_enh <- gene_enh |> filter(!is.na(Active_Enhancer)) 

## ENCODE4 Distal Enhancer-like ----
system("wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed")

dir.create("ENCODE_cell_types")
system("wget https://downloads.wenglab.org/Registry-V4/ENCFF429CVH_ENCFF061VAL_ENCFF573FQS.bed -O ENCODE_cell_types/Bcell_ENCFF429CVH_ENCFF061VAL_ENCFF573FQS.bed")
system("wget https://downloads.wenglab.org/Registry-V4/ENCFF641UWX_ENCFF507AEW_ENCFF272WIY.bed -O ENCODE_cell_types/CD8ab_ENCFF641UWX_ENCFF507AEW_ENCFF272WIY.bed")
system("wget https://downloads.wenglab.org/Registry-V4/ENCFF057LCQ_ENCFF579KFE_ENCFF810PFO.bed -O ENCODE_cell_types/effmemCD8_ENCFF057LCQ_ENCFF579KFE_ENCFF810PFO.bed")
system("wget https://downloads.wenglab.org/Registry-V4/ENCFF500QFS_ENCFF165CYH_ENCFF228CYK.bed -O ENCODE_cell_types/memB_ENCFF500QFS_ENCFF165CYH_ENCFF228CYK.bed")
system("wget https://downloads.wenglab.org/Registry-V4/ENCFF310ZPG_ENCFF837WOY_ENCFF926THG.bed -O ENCODE_cell_types/nB39_ENCFF310ZPG_ENCFF837WOY_ENCFF926THG.bed")
system("wget https://downloads.wenglab.org/Registry-V4/ENCFF207ATZ_ENCFF067FAK_ENCFF901BFR.bed -O ENCODE_cell_types/nB40_ENCFF207ATZ_ENCFF067FAK_ENCFF901BFR.bed")

ENC_cells_files <- list.files(path = "ENCODE_cell_types", pattern = "bed", full.names = T)


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

#2. Defining (intra)genic and coding regions ----
## Intragenic regions are defined by: exons, introns, 3'UTR and 5'UTR
## Coding regions are defined only by: exons

exons0 <- GenomicFeatures::exonsBy(txdb)
exons <- unlist(exons0)
exons$tx_id <- rep(names(exons0), elementNROWS(exons0))

introns0 <- GenomicFeatures::intronsByTranscript(txdb)
introns <- unlist(introns0)
introns$tx_id <- rep(names(introns0), elementNROWS(introns0))
introns$Feature <- "Intron"

UTR50 <- GenomicFeatures::fiveUTRsByTranscript(txdb)
UTR5 <-  unlist(UTR50)
UTR5$tx_id <- rep(names(UTR50), elementNROWS(UTR50))
UTR5$Feature <- "5' UTR"

UTR30 <- GenomicFeatures::threeUTRsByTranscript(txdb)
UTR3 <-  unlist(UTR30)
UTR3$tx_id <- rep(names(UTR30), elementNROWS(UTR30))
UTR3$Feature <- "3' UTR"

# codings <- c(exons) # Coding regions were only tested for the Reviewers point-by-point answers, everything involving coding regions analysis is commented.
genics <- c(exons, introns, UTR3, UTR5)

# codings$gene_id <- unlist(tx2gene$gene_id[match(codings$tx_id,tx2gene$tx_id)])
genics$gene_id <- unlist(tx2gene$gene_id[match(genics$tx_id,tx2gene$tx_id)])

## Add layer of only protein_coding genes
## This filter was only applied for the Reviewers point-by-point answers
# codings <- codings[codings$gene_id %in% promoters$ID,]
# genics <- genics[genics$gene_id %in% promoters$ID,]

# coding_l <- codings |> split(codings$gene_id)
genic_l <- genics |> split(genics$gene_id)

# coding_red_l <- parallel::mclapply(names(coding_l), function(x)
# {
#     red <- IRanges::reduce(coding_l[[x]])
#     red$gene_id <- x
#     return(red)
# }, mc.cores = 20)
#
# coding_red <- unlist(GRangesList(coding_red_l))
# coding_red <- sort(coding_red)
# seqlevelsStyle(coding_red) <- "NCBI"
#
# coding_prop <- coding_red %>% filter(seqnames %in% c(1:22,"X","Y"))
# strand(coding_prop) <- "*"

genic_red_l <- parallel::mclapply(names(genic_l), function(x)
{
    red2 <- IRanges::reduce(genic_l[[x]])
    red2$gene_id <- x
    return(red2)
}, mc.cores = 20)

genic_red <- unlist(GRangesList(genic_red_l))
genic_red <- sort(genic_red)
seqlevelsStyle(genic_red) <- "NCBI"

genic_prop <- genic_red %>% filter(seqnames %in% c(1:22,"X","Y"))
strand(genic_prop) <- "*"


#3. Types of Enhancers ----
## Active Enhancers - gene pairs (AE-ATAC-c) ----

active_enh <- GRanges(gene_enh$Active_Enhancer) |> mutate(cell=gene_enh$cell) |> 
    unique() |> sort() %>% mutate(AE_ID=1:length(.))

sum(width(reduce(active_enh))) ## 194602 elements - 50907275 bp


## All AE and AE-ATAC ----
cells <- unique(gene_enh$cell)

all_AE <- GRanges()
all_AE_ATAC <- GRanges()
for (cell in cells) 
{
    print(cell)
    chromhmm_file <- grep(paste0("_",cell), chromHMM_files, value = T)
    ATAC_file <- grep(paste0("_",cell), ATAC_files, value = T)
    
    chromhmmGR <- data.table::fread(file = chromhmm_file)[,1:4] |> 
        makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3",keep.extra.columns = T) |> 
        sort()    
    colnames(mcols(chromhmmGR)) <- "States"
    
    AE_GR <- chromhmmGR |> plyranges::filter(States %in% 4:5) |> GenomicRanges::reduce() %>% mutate(ID=paste(cell,"AE",1:length(.), sep = "_"))

    ATAC_GR <- data.table::fread(file = ATAC_file) |> 
        makeGRangesFromDataFrame(seqnames.field = "V1",start.field = "V2",end.field = "V3")
    
    AE_ATAC <- intersect(AE_GR, ATAC_GR)%>% mutate(ID=paste(cell,"AE_ATAC",1:length(.), sep = "_"))
    
    all_AE <- c(all_AE,AE_GR)
    all_AE_ATAC <- c(all_AE_ATAC,AE_ATAC)
    
}

sum(width(IRanges::reduce(all_AE_ATAC))) ## 714303 elements - 124956327 bp
sum(width(reduce(all_AE))) ## 716005 elements - 342537000 bp

## ENCODE4 enhancers (ENCODE-E) ----
enh_ENC <- data.table::fread("GRCh38-cCREs.ELS.bed") |> 
    makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T) %>% 
    mutate(AE_ID=1:length(.))
seqlevelsStyle(enh_ENC) <- "NCBI"

sum(width(reduce(enh_ENC))) ## 1718669 elements - 459966796 bp


enh_ENC_cell <- list()
for (file in ENC_cells_files)
{
    cell <- gsub("_.*","", basename(file))
    enh <- data.table::fread(file) |> filter(grepl("ELS",V10)) |> 
        makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T) %>%
        mutate(AE_ID=1:length(.))
    seqlevelsStyle(enh) <- "NCBI"
    
    enh_ENC_cell[[cell]] <- enh
}


# 4. Test intragenic proportion ----

df_perc <- data.frame()
# df_perc_coding <- data.frame()

## Active Enhancers (AE-ATAC-c)----

active_enh_genic <- active_enh |> join_overlap_left(genic_red) |> filter(!is.na(gene_id)) |> 
    as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)

perc <- sum(unique(active_enh$AE_ID) %in% unique(active_enh_genic$AE_ID))/length(unique(active_enh$AE_ID)) 
number <- length(unique(active_enh_genic$AE_ID)) 

df_perc <- rbind(df_perc,c("all",perc, number, "AE-ATAC-c"))

# active_enh_coding <- active_enh |> join_overlap_left(coding_red) |> filter(!is.na(gene_id)) |> 
#     as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
# 
# perc <- sum(unique(active_enh$AE_ID) %in% unique(active_enh_coding$AE_ID))/length(unique(active_enh$AE_ID))
# number <- length(unique(active_enh_coding$AE_ID))
# 
# df_perc_coding <- rbind(df_perc_coding,c("all",perc, number, "AE-ATAC-c"))

### by cell type ----

active_enh_cell <- active_enh %>% split(.$cell)
active_enh_cell_genic <- lapply(active_enh_cell, function(x){
    x |> join_overlap_left(genic_red) |> filter(!is.na(gene_id)) |> 
        as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
})
# active_enh_cell_coding <- lapply(active_enh_cell, function(x){
#     x |> join_overlap_left(coding_red) |> filter(!is.na(gene_id)) |> 
#         as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
# })

for (cell in names(active_enh_cell)) 
{
    perc <- sum(unique(active_enh_cell[[cell]]$AE_ID) %in% unique(active_enh_cell_genic[[cell]]$AE_ID))/length(unique(active_enh_cell[[cell]]$AE_ID))
    number <- length(unique(active_enh_cell[[cell]]$AE_ID)) 
    df_perc <- rbind(df_perc,c(cell,perc, number, "AE-ATAC-c"))
    
    # perc <- sum(unique(active_enh_cell[[cell]]$AE_ID) %in% unique(active_enh_cell_coding[[cell]]$AE_ID))/length(unique(active_enh_cell[[cell]]$AE_ID))
    # number <- length(unique(active_enh_cell[[cell]]$AE_ID)) 
    # df_perc_coding <- rbind(df_perc_coding,c(cell,perc, number, "AE-ATAC-c"))
}


## Act - ATAC - lichic (AE) ----

all_AE <- all_AE %>% 
    mutate(AE_ID=1:length(.))

all_AE_genic <- all_AE |> join_overlap_left(genic_red) |> filter(!is.na(gene_id)) |>
    as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)

perc <- sum(unique(all_AE$AE_ID) %in% unique(all_AE_genic$AE_ID))/length(unique(all_AE$AE_ID)) 
number <- length(unique(all_AE_genic$AE_ID)) 

df_perc <- rbind(df_perc,c("all",perc, number, "AE"))

# all_AE_coding <- all_AE |> join_overlap_left(coding_red) |> filter(!is.na(gene_id)) |> 
#     as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
# 
# perc <- sum(unique(all_AE$AE_ID) %in% unique(all_AE_coding$AE_ID))/length(unique(all_AE$AE_ID)) 
# number <- length(unique(all_AE_coding$AE_ID)) 
# 
# df_perc_coding <- rbind(df_perc_coding,c("all",perc, number, "AE"))

### by cell type ----

all_AE_cell <- all_AE |> mutate(cell = gsub("_.*E.*","",ID)) %>% split(.$cell)
all_AE_cell_genic <- lapply(all_AE_cell, function(x){
    x |> join_overlap_left(genic_red) |> filter(!is.na(gene_id)) |> 
        as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
})

# all_AE_cell_coding <- lapply(all_AE_cell, function(x){
#     x |> join_overlap_left(coding_red) |> filter(!is.na(gene_id)) |> 
#         as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
# })

for (cell in names(all_AE_cell)) 
{
    perc <- sum(unique(all_AE_cell[[cell]]$AE_ID) %in% unique(all_AE_cell_genic[[cell]]$AE_ID))/length(unique(all_AE_cell[[cell]]$AE_ID))
    number <- length(unique(all_AE_cell[[cell]]$AE_ID)) 
    df_perc <- rbind(df_perc,c(cell,perc, number, "AE"))
    
    # perc <- sum(unique(all_AE_cell[[cell]]$AE_ID) %in% unique(all_AE_cell_coding[[cell]]$AE_ID))/length(unique(all_AE_cell[[cell]]$AE_ID))
    # number <- length(unique(all_AE_cell[[cell]]$AE_ID)) 
    # df_perc_coding <- rbind(df_perc_coding,c(cell,perc, number, "AE"))
}

## Act + ATAC - lichic (AE-ATAC) ----

all_AE_ATAC <- all_AE_ATAC %>% 
    mutate(AE_ID=1:length(.))

all_AE_ATAC_genic <- all_AE_ATAC |> join_overlap_left(genic_red) |> filter(!is.na(gene_id)) |> 
    as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)

perc <- sum(unique(all_AE_ATAC$AE_ID) %in% unique(all_AE_ATAC_genic$AE_ID))/length(unique(all_AE_ATAC$AE_ID)) 
number <- length(unique(all_AE_ATAC_genic$AE_ID)) 

df_perc <- rbind(df_perc,c("all",perc, number, "AE-ATAC"))

# all_AE_ATAC_coding <- all_AE_ATAC |> join_overlap_left(coding_red) |> filter(!is.na(gene_id)) |> 
#     as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
# 
# perc <- sum(unique(all_AE_ATAC$AE_ID) %in% unique(all_AE_ATAC_coding$AE_ID))/length(unique(all_AE_ATAC$AE_ID)) 
# number <- length(unique(all_AE_ATAC_coding$AE_ID)) 
# 
# df_perc_coding <- rbind(df_perc_coding,c("all",perc, number, "AE-ATAC"))

### by cell type ----

all_AE_ATAC_cell <- all_AE_ATAC |> mutate(cell = gsub("_.*E.*","",ID)) %>% split(.$cell)
all_AE_ATAC_cell_genic <- lapply(all_AE_ATAC_cell, function(x){
    x |> join_overlap_left(genic_red) |> filter(!is.na(gene_id)) |> 
        as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
})

# all_AE_ATAC_cell_coding <- lapply(all_AE_ATAC_cell, function(x){
#     x |> join_overlap_left(coding_red) |> filter(!is.na(gene_id)) |> 
#         as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
# })

for (cell in names(all_AE_ATAC_cell)) 
{
    perc <- sum(unique(all_AE_ATAC_cell[[cell]]$AE_ID) %in% unique(all_AE_ATAC_cell_genic[[cell]]$AE_ID))/length(unique(all_AE_ATAC_cell[[cell]]$AE_ID))
    number <- length(unique(all_AE_ATAC_cell[[cell]]$AE_ID)) 
    df_perc <- rbind(df_perc,c(cell,perc, number, "AE-ATAC"))
    
    # perc <- sum(unique(all_AE_ATAC_cell[[cell]]$AE_ID) %in% unique(all_AE_ATAC_cell_coding[[cell]]$AE_ID))/length(unique(all_AE_ATAC_cell[[cell]]$AE_ID))
    # number <- length(unique(all_AE_ATAC_cell[[cell]]$AE_ID)) 
    # df_perc_coding <- rbind(df_perc_coding,c(cell,perc, number, "AE-ATAC"))
}

## ENCODE (ENCODE-E) ----
enh_ENC_genic <- enh_ENC |> join_overlap_left(genic_red) |> filter(!is.na(gene_id)) |> 
    as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)

perc <- sum(unique(enh_ENC$AE_ID) %in% unique(enh_ENC_genic$AE_ID))/length(unique(enh_ENC$AE_ID)) 
number <- length(unique(enh_ENC_genic$AE_ID)) 

df_perc <- rbind(df_perc,c("all",perc, number, "ENCODE-E"))

# enh_ENC_coding <- enh_ENC |> join_overlap_left(coding_red) |> filter(!is.na(gene_id)) |> 
#     as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
# 
# perc <- sum(unique(enh_ENC$AE_ID) %in% unique(enh_ENC_coding$AE_ID))/length(unique(enh_ENC$AE_ID)) 
# number <- length(unique(enh_ENC_coding$AE_ID)) 
# 
# df_perc_coding <- rbind(df_perc_coding,c("all",perc, number, "ENCODE-E"))

### by cell type ----

enh_ENC_cell_genic <- lapply(enh_ENC_cell, function(x){
    x |> join_overlap_left(genic_red) |> filter(!is.na(gene_id)) |> 
        as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
})

# enh_ENC_cell_coding <- lapply(enh_ENC_cell, function(x){
#     x |> join_overlap_left(coding_red) |> filter(!is.na(gene_id)) |> 
#         as_tibble() |> dplyr::select(AE_ID, gene_id) |> dplyr::rename(Gene_overlap = gene_id)
# })

for (cell in names(enh_ENC_cell)) 
{
    perc <- sum(unique(enh_ENC_cell[[cell]]$AE_ID) %in% unique(enh_ENC_cell_genic[[cell]]$AE_ID))/length(unique(enh_ENC_cell[[cell]]$AE_ID))
    number <- length(unique(enh_ENC_cell[[cell]]$AE_ID)) 
    df_perc <- rbind(df_perc,c(cell,perc, number, "ENCODE-E"))
    
    # perc <- sum(unique(enh_ENC_cell[[cell]]$AE_ID) %in% unique(enh_ENC_cell_coding[[cell]]$AE_ID))/length(unique(enh_ENC_cell[[cell]]$AE_ID))
    # number <- length(unique(enh_ENC_cell[[cell]]$AE_ID)) 
    # df_perc_coding <- rbind(df_perc_coding,c(cell,perc, number, "ENCODE-E"))
}

#5. Check which gene is regulating ----
active_enh <- GRanges(gene_enh$Active_Enhancer) |> 
    unique() |> sort() %>% mutate(AE_ID=1:length(.))

gene_enh$AE_ID <- active_enh$AE_ID[match(GRanges(gene_enh$Active_Enhancer),active_enh)]

gene_enh <- gene_enh |> 
    filter(Gene_biotype == "protein_coding")

## check AEs that target a gene different that the one that overlap
active_enh_genic_prot <- active_enh_genic |> filter(AE_ID %in% unique(gene_enh$AE_ID))
active_enh_prot <- active_enh |> filter(AE_ID %in% unique(gene_enh$AE_ID))


enh_target_other <- gene_enh |> dplyr::select(AE_ID,Gene_ID, Gene_Name) |> unique() |> right_join(active_enh_genic_prot) |> 
    group_by(AE_ID) |> reframe(other_gene=any(!Gene_ID %in% Gene_overlap),
                               same_gene=any(Gene_ID %in% Gene_overlap))

l <- list(intergenic=unique(active_enh_prot$AE_ID[!active_enh_prot$AE_ID %in% unique(active_enh_genic_prot$AE_ID)]),
          intragenic=unique(active_enh_genic_prot$AE_ID),
          target_same_gene=unique(enh_target_other$AE_ID[enh_target_other$same_gene]),
          target_other_gene=unique(enh_target_other$AE_ID[enh_target_other$other_gene]))

df <- t(as.data.frame(lapply(attr(gplots::venn(l),"intersections"), length))) |> as.data.frame() |> 
    mutate(genic=factor(c("Intergenic","Intragenic","Intragenic","Intragenic"), levels=c("Intragenic","Intergenic")),
           target=factor(c("","Host","Other","Both"), levels=c("Other","Both","Host",""))) |> 
    dplyr::rename(value=V1)


#6. Plotting ----
## Intragenic ----

## Figure 1C
df |> ggplot(aes(x="AE", y=value, fill=interaction(target,genic))) + 
    geom_bar(stat = "identity", position = "fill", col="white") + 
    ggpubr::theme_classic2() + 
    labs(x=NULL, y = "Proportion of AE", fill=NULL) + 
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_fill_manual(values = c(Other.Intragenic="#c4053cff",
                                 Both.Intragenic="#8c476bff",
                                 Host.Intragenic="#b55d44ff",
                                 .Intergenic="#e2d1bfff"))


colnames(df_perc) <- c("cell","perc","number","type")
df_perc <- df_perc |> mutate(perc=as.numeric(perc),
                             number=as.numeric(number),
                             type=factor(type, levels = c("AE-ATAC-c","AE-ATAC","AE","ENCODE-E")))

## Extended Data Figure 2D (subset of point-by-point figure IVa)
ggplot(rbind(df_perc |> filter(cell == "all") |> mutate(type2="intragenic"),
             df_perc |> filter(cell == "all") |> mutate(perc = 1-perc,
                                                        type2="intergenic")) |> 
           mutate(type2=factor(type2, levels = c("intragenic","intergenic"))), 
       aes(x=type,y=perc, fill=type2)) + 
    geom_bar(stat="identity", position = "stack") + 
    labs(x=NULL,y="Proportion of AE", fill=NULL) + 
    ggpubr::theme_classic2() + 
    theme(axis.text.x = element_text(angle = 90)) +
    ylim(c(0,1)) + 
    scale_fill_manual(values = c(intragenic="#c4053cff",
                                 intergenic="#e2d1bfff"))

## Point-by-point figure IVc
ggplot(df_perc |> filter(cell != "all"), aes(x=type,y=perc, label=cell)) + 
    geom_boxplot(outlier.shape = NA) + ggbeeswarm::geom_beeswarm(size=1) + 
    labs(x=NULL,y="Proportion of intragenic AE") + 
    ggpubr::theme_classic2() + 
    theme(axis.text.x = element_text(angle = 90)) +
    ylim(c(0,1))


## coding ----
# colnames(df_perc_coding) <- c("cell","perc","number","type")
# df_perc_coding <- df_perc_coding |> mutate(perc=as.numeric(perc),
#                                            number=as.numeric(number),
#                                            type=factor(type, levels = c("AE-ATAC-c","AE-ATAC","AE","ENCODE-E")))

## Point-by-point figure IVb
# ggplot(rbind(df_perc_coding |> filter(cell == "all") |> mutate(type2="coding"),
#              df_perc_coding |> filter(cell == "all") |> mutate(perc = 1-perc,
#                                                                type2="non_coding")) |> 
#            mutate(type2=factor(type2, levels = c("coding","non_coding"))), 
#        aes(x=type,y=perc, fill=type2)) + 
#     geom_bar(stat="identity", position = "stack") + 
#     labs(x=NULL,y="Proportion of AE", fill=NULL) + 
#     ggpubr::theme_classic2() + 
#     theme(axis.text.x = element_text(angle = 90)) +
#     ylim(c(0,1)) + 
#     scale_fill_manual(values = c(coding="#c4053cff",
#                                  non_coding="#c98467ff"))


