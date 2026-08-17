library(tidyverse)

# Basic statistics of Enhancers (number and coverage)

# 1. Data ----
## Locate data necessary for this analysis

AE_files <- list.files(pattern = "Active_Enhancers.bed", full.names = T, recursive = F)
PrE_files <- list.files(pattern = "Primed_Enhancers.bed", full.names = T, recursive = F)
PoE_files <- list.files(pattern = "Poised_Enhancers.bed", full.names = T, recursive = F)
SE_files <- list.files(pattern = "Super_Enhancers.bed", full.names = T, recursive = F)

# 2. Get the number of each type of Enhancer in each cell type and plot ----

cells <- gsub("_.*","",basename(AE_files))

df_final <- tibble()
for (cell in cells) 
{
    active <- data.table::fread(grep(paste0("/",cell),AE_files,value = T))
    primed <- data.table::fread(grep(paste0("/",cell),PrE_files,value = T))
    poised <- data.table::fread(grep(paste0("/",cell),PoE_files,value = T))
    super <- data.table::fread(grep(paste0("/",cell),SE_files,value = T))

    df <- tibble(cell=cell,AE=nrow(active),SE=nrow(super),PrE=nrow(primed),PoE=nrow(poised))
    df_final <- rbind(df_final,df)
}

colors <- c(AE="#c1153bff",
            SE="#e38814ff",
            PrE="#b55d44ff",
            PoE="#8a486cff")

### Figure 1D
df_final %>% pivot_longer(cols = -cell) %>% mutate(name=factor(name,levels = c("AE","SE","PrE","PoE")),
                                                   cell=factor(cell,levels = c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC", "CMP", "Mon", "nCD8"))) %>%
    ggplot(aes(x=name,y=value/1000, fill=name)) + geom_bar(stat = "identity", position = "dodge") + 
    facet_grid(~cell) + 
    ggpubr::theme_classic2(base_size = 5) + 
    labs(x=NULL,y="Number of Element (k)",fill="Type") +
    theme(axis.text.x = element_blank()) + 
    scale_fill_manual(values = colors)

# 3. Compute the total coverage of each type of Enhancer in each cell type and plot ----

cov_final <- tibble()
for (cell in cells) 
{
    print(cell)
    active <- GenomicRanges::makeGRangesFromDataFrame(data.table::fread(grep(paste0("/",cell),AE_files,value = T)))
    primed <- GenomicRanges::makeGRangesFromDataFrame(data.table::fread(grep(paste0("/",cell),PrE_files,value = T)))
    poised <- GenomicRanges::makeGRangesFromDataFrame(data.table::fread(grep(paste0("/",cell),PoE_files,value = T)))
    super <- GenomicRanges::makeGRangesFromDataFrame(data.table::fread(grep(paste0("/",cell),SE_files,value = T)))
    
    cov <- tibble(cell=cell,AE=sum(width(active)),SE=sum(width(super)),PrE=sum(width(primed)),PoE=sum(width(poised)))
    cov_final <- rbind(cov_final,cov)
}

### Extended Data Figure 2C
cov_final %>% pivot_longer(cols = -cell) %>% mutate(name=factor(name,levels = c("AE","SE","PrE","PoE")),
                                                    cell=factor(cell,levels = c("HSC", "PreProB", "ProB", "PreB", "immtransB", "nB", "GCB", "memB", "PC", "CMP", "Mon", "nCD8"))) %>%
    ggplot(aes(x=name,y=value/1000000, fill=name)) + geom_bar(stat = "identity", position = "dodge") + 
    facet_grid(~cell) + 
    ggpubr::theme_classic2(base_size = 5) + 
    labs(x=NULL,y="Coverage (Mb)",fill="Type") +
    theme(axis.text.x = element_blank()) + 
    scale_fill_manual(values = colors)