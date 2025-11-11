library(umap)
library(Rtsne)
library(data.table)
library(ggplot2)
library(pheatmap)
library(ggplotify)
library(grid)
library(preprocessCore)
library(readxl)
library(tidyverse)
library(edgeR)

#specify counts file
counts <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/unfiltered_analysis/batch_correction/batch_effect/quant_06_16_25/counts_w_neg.csv"

#specify metadata file
metadata <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/metadata/skin_layer/all_samples.csv"

#specify output directory
output_dir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/final/NG2_10_09_25/UMAP_heat_quant"


#specify whether I want to add voom normalization and/or quantile normalization
perform_voom <- FALSE
perform_quant <- TRUE

#specify whether I want to filter down to immune-related genes
immune_genes <- TRUE

#set working directory to the output folder, create it if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)}

setwd(output_dir)

#clean the counts data up
counts <- as.data.frame(fread(counts))
rownames(counts) <- counts$V1
counts$V1 <- NULL
dcc_kept <- sub("\\.dcc$", "", colnames(counts))
colnames(counts) <- dcc_kept

#subset the metadata down to the ROIs kept after QC
metadata <- as.data.frame(fread(metadata))
metadata <- subset(metadata, metadata$Sample_ID %in% dcc_kept)
metadata <- metadata[match(dcc_kept, metadata$Sample_ID), ]
metadata$Sample_ID == dcc_kept
metadata$V1 <- NULL

#remove the samples that do not have custom diagnosis assigned 
metadata <- subset(metadata, !(is.na(Custom_Diagnosis)))
counts <- counts[, metadata$Sample_ID]

#remove negative controls from the counts data
counts <- counts[!grepl("^RTS", rownames(counts)), ]

probes <- rownames(counts)
samples <- colnames(counts)

################# Normalization Steps #################

if (immune_genes == TRUE){
    genes_to_annotate <- probes

    genes_immune <- read_xls("/Users/elise/Downloads/innatedb_curated_genes.xls", col_names = TRUE)
    genes_imm <- read_xlsx("/Users/elise/Downloads/InnateDB_genes.xlsx")
    septic_genes <- read_xlsx("/Users/elise/Downloads/septic_shock_genes.xlsx")

    genes_imm <- genes_imm$name
    septic_genes <- septic_genes$name

    genes_immune <- genes_immune[, "Gene Symbol"]
    genes_immune <- genes_immune %>% as.data.frame() %>% unique()
    symbols_immune <- genes_immune[, "Gene Symbol"]

    #which ones are in my list of genes?
    all_genes_immune <- c(symbols_immune, genes_imm, septic_genes)
    all_genes_immune <- toupper(all_genes_immune)
    all_genes_immune <- unique(all_genes_immune)

    my_genes_immune <- intersect(genes_to_annotate, all_genes_immune)
    length(my_genes_immune) #down to 2451 genes that are in the innate immunity, immunity, and septic shock databases 

   #counts <- counts[my_genes_immune, ]
}

#voom log cpm counts
if (perform_voom == TRUE){
   counts_dge <- DGEList(counts)
   v <- calcNormFactors(counts_dge)
   v <- voom(v)
    #corfit <- duplicateCorrelation(v, block = metadata$Slide.Name)
    #v <- voom(counts, block = metadata$Slide.Name, correlation = corfit$consensus)

   voom_norm <- v$E %>% as.data.frame()
   counts <- voom_norm
}

#quantile normalize using the log cpm counts
if (perform_quant == TRUE){
    prep_counts <- as.matrix(counts)
    quant_norm <- normalize.quantiles(prep_counts) %>% as.data.frame()
    counts <- quant_norm
}

colnames(counts) <- samples
rownames(counts) <- probes

#subset down to immune genes only
counts <- counts[my_genes_immune, ]

#save these counts and metadata
counts_save <- counts
metadata_save <- metadata

################# Plotting ##############################

#create UMAP, tSNE, and gene expression heatmap for one cell type
UMAP_tSNE_heatmap <- function(metadata, counts, factor, filter) {

    dir.create(filter)

    custom_umap <- umap::umap.defaults
    custom_umap$random_state <- 42
    custom_umap$n_neighbors <- 5

    umap_out <-
        umap(t(counts),
            config = custom_umap)
    metadata[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

    #########UMAP##########
    UMAP <- ggplot(metadata,
       aes(x = UMAP1, y = UMAP2, color = metadata[[factor]])) +
    geom_point(size = 3) + 
    theme_bw()

    ggsave(paste(filter,"/UMAP_",factor,".png", sep = ""),
       UMAP, width = 10, height = 8, dpi = 300)

    #########t-SNE##########
    tsne_out <-
       Rtsne(t(counts),
           perplexity = ncol(counts)*.15)
    metadata[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
    tSNE <- ggplot(metadata, aes(x = tSNE1, y = tSNE2, color = metadata[[factor]])) +
       geom_point(size = 3, alpha = 0.6) +  # Slight transparency for clarity 
       theme_bw()

    ggsave(paste(filter, "/tSNE_",factor,".png", sep = ""),
       tSNE, width = 10, height = 8, dpi = 300)

    #########heatmap##########
    annotation <- metadata[, c("Sample_ID", factor), drop = FALSE]
    annotation <- annotation[order(annotation[, factor]), ]
    rownames(annotation) <- annotation[, "Sample_ID"]
    annotation[, "Sample_ID"] <- NULL
    counts <- counts[ , rownames(annotation)]

    calc_CV <- function(x) {sd(x) / mean(x)}
    CV_dat <- apply(counts, MARGIN = 1, FUN = calc_CV)
    sort(CV_dat, decreasing = TRUE)[1:30]

    GOI <- names(CV_dat)[CV_dat >= quantile(CV_dat, 0.9)]
    counts <- counts[GOI, ]

    heatmap_plot <- as.ggplot(pheatmap(counts,
                                   scale = "row",
                                   show_rownames = FALSE, show_colnames = FALSE,
                                   border_color = NA,
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   cluster_cols = FALSE,
                                   breaks = seq(-3, 3, 0.05),
                                   bg = "white",
                                   color = colorRampPalette(c("#0a71cb", "#FFFFFF", "#D1495B"))(120),
                                   annotation_col = annotation, 
                                   legend = TRUE,
                                   legend.title = list(gpar(fontsize = 14)),
                                   legend.text = list(gpar(fontsize = 12))))
    ggsave(paste(filter, "/ord_heatmap_",factor,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)

    heatmap_plot <- as.ggplot(pheatmap(log2(counts),
                                   scale = "row",
                                   show_rownames = FALSE, show_colnames = FALSE,
                                   border_color = NA,
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   cluster_cols = TRUE,
                                   breaks = seq(-3, 3, 0.05),
                                   bg = "white",
                                   color = colorRampPalette(c("#0a71cb", "#FFFFFF", "#D1495B"))(120),
                                   annotation_col = annotation, 
                                   legend = TRUE,
                                   legend.title = list(gpar(fontsize = 14)),
                                   legend.text = list(gpar(fontsize = 12))))
    ggsave(paste(filter, "/notord_heatmap_",factor,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)
}

#create UMAP, tSNE, and gene expression heatmap for only two diagnosis groups, from one cell type
#this includes confidence interval ellipses around each group
UMAP_tSNE_heatmap_2groups <- function(metadata, counts, factor, filter, groups, plot_name) {

    dir.create(paste0(output_dir, "/", filter))

    custom_umap <- umap::umap.defaults
    custom_umap$random_state <- 42
    umap_out <-
        umap(t(counts),
            config = custom_umap)
    metadata[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

    #########UMAP##########
    UMAP <- ggplot(metadata, aes(x = UMAP1, y = UMAP2, color = .data[[factor]])) +
        geom_point(size = 3) +
        stat_ellipse(data = metadata[metadata[[factor]] %in% groups, ],
                 aes(fill = .data[[factor]]),
                 alpha = 0.2, geom = "polygon") +
        theme_bw()

    ggsave(paste(filter,"/UMAP_",plot_name,".png", sep = ""),
       UMAP, width = 10, height = 8, dpi = 300)

    #########t-SNE##########
    tsne_out <-
       Rtsne(t(counts),
           perplexity = ncol(counts)*.15)
    metadata[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
    tSNE <- ggplot(metadata, aes(x = tSNE1, y = tSNE2, color = .data[[factor]])) +
        geom_point(size = 3, alpha = 0.6) +
        stat_ellipse(data = metadata[metadata[[factor]] %in% groups, ],
                 aes(fill = .data[[factor]]),
                 alpha = 0.2, geom = "polygon") +
        theme_bw()

    ggsave(paste(filter, "/tSNE_",plot_name,".png", sep = ""),
       tSNE, width = 10, height = 8, dpi = 300)

    #########heatmap##########
    md <- subset(metadata, metadata[[factor]] %in% groups)
    counts <- counts[, colnames(counts) %in% md$Sample_ID]

    annotation <- md[, c("Sample_ID", factor), drop = FALSE]
    annotation <- annotation[order(annotation[, factor]), ]
    rownames(annotation) <- annotation[, "Sample_ID"]
    annotation[, "Sample_ID"] <- NULL
    counts <- counts[ , rownames(annotation)]

    calc_CV <- function(x) {sd(x) / mean(x)}
    CV_dat <- apply(counts, MARGIN = 1, FUN = calc_CV)
    sort(CV_dat, decreasing = TRUE)[1:30]

    GOI <- names(CV_dat)[CV_dat >= quantile(CV_dat, 0.9)]
    counts <- counts[GOI, ]

    heatmap_plot <- as.ggplot(pheatmap(counts,
                                   scale = "row",
                                   show_rownames = FALSE, show_colnames = FALSE,
                                   border_color = NA,
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   cluster_cols = FALSE,
                                   breaks = seq(-3, 3, 0.05),
                                   bg = "white",
                                   color = colorRampPalette(c("#0a71cb", "#FFFFFF", "#D1495B"))(120),
                                   annotation_col = annotation, 
                                   legend = TRUE,
                                   legend.title = list(gpar(fontsize = 14)),
                                   legend.text = list(gpar(fontsize = 12))))
    ggsave(paste(filter, "/ord_heatmap_",plot_name,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)

    heatmap_plot <- as.ggplot(pheatmap(counts,
                                   scale = "row",
                                   show_rownames = FALSE, show_colnames = FALSE,
                                   border_color = NA,
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   cluster_cols = TRUE,
                                   breaks = seq(-3, 3, 0.05),
                                   bg = "white",
                                   color = colorRampPalette(c("#0a71cb", "#FFFFFF", "#D1495B"))(120),
                                   annotation_col = annotation, 
                                   legend = TRUE,
                                   legend.title = list(gpar(fontsize = 14)),
                                   legend.text = list(gpar(fontsize = 12))))
    ggsave(paste(filter, "/notord_heatmap_",plot_name,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)
}

#example implementations: 

##### dermal stroma ###### 2 groups
metadata_ds <- subset(metadata, metadata$Custom_CellType == "Dermal Stroma") #subset metadata down to dermal stroma only
counts_ds <- counts[, metadata_ds$Sample_ID] #subset counts down to dermal stroma only
groups <- c("CTCL", "HC", "AD") #groups I want outlined with a confidence interval
plot_name <- "CTCL_HC_AD" #how I want the plot to be named
UMAP_tSNE_heatmap_2groups(metadata_ds, counts_ds, "Custom_Diagnosis", "Dermal_Stroma", groups, plot_name) #use UMAP_tSNE_heatmap_2groups to show confidence intervals

#we're still doing dermal stroma, but now I just want CTCL and HC outlined
groups <- c("CTCL", "HC") #use CTCL and HC as groups
plot_name <- "CTCL_HC" #name the plot
UMAP_tSNE_heatmap_2groups(metadata_ds, counts_ds, "Custom_Diagnosis", "Dermal_Stroma", groups, plot_name) #plot with confidence intervals

##### keratinocytes ######
metadata_kt <- subset(metadata, metadata$Custom_CellType == "Keratinocyte") #start with the original metadata and subset down to keratinocytes
counts_kt <- counts[, metadata_kt$Sample_ID] #subset counts down to match
groups <- c("DR", "CTCL") #just want confidence intervals on Drug Rash and CTCL
plot_name <- "DR_CTCL" #name the plot
UMAP_tSNE_heatmap_2groups(metadata_kt, counts_kt, "Custom_Diagnosis", "Keratinocyte", groups, plot_name) #plot with confidence intervals

groups <- c("AD", "CTCL") #just AD and CTCL
plot_name <- "AD_CTCL" #name the plot
UMAP_tSNE_heatmap_2groups(metadata_kt, counts_kt, "Custom_Diagnosis", "Keratinocyte", groups, plot_name) #plot with confidence intervals

groups <- c("PSO", "CTCL") #just PSO and CTCL
plot_name <- "PSO_CTCL" #name the plot
UMAP_tSNE_heatmap_2groups(metadata_kt, counts_kt, "Custom_Diagnosis", "Keratinocyte", groups, plot_name) #plot with confidence intervals 


##### dermal stroma ###### all groups
metadata_ds <- subset(metadata, metadata$Custom_CellType == "Dermal Stroma")
counts_ds <- counts[, metadata_ds$Sample_ID]
UMAP_tSNE_heatmap(metadata_ds, counts_ds, "Custom_Diagnosis", "Dermal_Stroma")

##### keratinocytes ###### all groups
metadata_kt <- subset(metadata, metadata$Custom_CellType == "Keratinocyte")
counts_kt <- impute_zero[, metadata_kt$Sample_ID]
UMAP_tSNE_heatmap(metadata_kt, counts_kt, "Custom_Diagnosis", "Keratinocyte")

#now, I have a function to plot UMAP by other factors
#this helps us check for batch effects or other technical bias
umap_by_factor <- function(metadata, counts, factor, plot_add){
   custom_umap <- umap::umap.defaults
   custom_umap$random_state <- 42
   custom_umap$n_neighbors <- 5

   umap_out <-
      umap(t(counts),
        config = custom_umap)
   metadata[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

   #########UMAP##########
   UMAP <- ggplot(metadata, aes(x = UMAP1, y = UMAP2, color = .data[[factor]])) +
      geom_point(size = 3) + 
      theme_bw()
   
   dir.create(paste0(output_dir, "/", factor))

   ggsave(paste(factor,"/UMAP", plot_add,".png", sep = ""),
       UMAP, width = 10, height = 8, dpi = 300)
}

#Slide by Slide UMAPs
umap_by_factor(metadata_kt, counts_kt, factor = "Slide.Name", plot_add = "ker")
umap_by_factor(metadata_ds, counts_ds, factor = "Slide.Name", plot_add = "ds")

#Patient by Patient UMAPs
umap_by_factor(metadata_kt, counts_kt, factor = "DP", plot_add = "ker")
umap_by_factor(metadata_ds, counts_ds, factor = "DP", plot_add = "ds")


####################################################
#heatmap of just the CTCLvsNonCTCL differentially expressed genes

#before this, you must run some sort of differential expression testing using DESeq2, limma-voom, or other tool/test

#load in differential expression results file that I want to analyze 
CTCLvsOther_ker <- fread("/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/final/NG2_comm_meeting/DEG_analysis/CTCL_Other_Keratinocyte_TRUE.csv") %>% as.data.frame()

genes_to_plot <- as.vector(CTCLvsOther_ker$V1)

counts_kt_CTCLvsOther <- counts_kt[genes_to_plot, ]

factor = "CTCLvsOther"

all_diag <- as.vector(metadata_kt$Custom_Diagnosis)
CTCLvsOther <- all_diag == "CTCL"
CTCLvsOther <- ifelse(CTCLvsOther, "CTCL", "Other")
metadata_kt$CTCLvsOther <- CTCLvsOther

md <- metadata_kt

annotation <- md[, c("Sample_ID", factor), drop = FALSE]
annotation <- annotation[order(annotation[, factor]), ]
rownames(annotation) <- annotation[, "Sample_ID"]
annotation[, "Sample_ID"] <- NULL
counts_kt_CTCLvsOther <- counts_kt_CTCLvsOther[ , rownames(annotation)]

calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- apply(counts_kt_CTCLvsOther, MARGIN = 1, FUN = calc_CV)
sort(CV_dat, decreasing = TRUE)[1:30]

GOI <- names(CV_dat)[CV_dat >= quantile(CV_dat, 0.9)]
counts <- counts[GOI, ]

heatmap_plot <- as.ggplot(pheatmap(counts_kt_CTCLvsOther,
                                   scale = "row",
                                   show_rownames = TRUE, show_colnames = FALSE,
                                   border_color = NA,
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   cluster_cols = FALSE,
                                   breaks = seq(-3, 3, 0.05),
                                   bg = "white",
                                   color = colorRampPalette(c("#0a71cb", "#FFFFFF", "#D1495B"))(120),
                                   annotation_col = annotation, 
                                   legend = TRUE,
                                   legend.title = list(gpar(fontsize = 14)),
                                   legend.text = list(gpar(fontsize = 12))))
    ggsave(paste(filter, "/ord_heatmap_",plot_name,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)

    heatmap_plot <- as.ggplot(pheatmap(counts,
                                   scale = "row",
                                   show_rownames = FALSE, show_colnames = FALSE,
                                   border_color = NA,
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   cluster_cols = TRUE,
                                   breaks = seq(-3, 3, 0.05),
                                   bg = "white",
                                   color = colorRampPalette(c("#0a71cb", "#FFFFFF", "#D1495B"))(120),
                                   annotation_col = annotation, 
                                   legend = TRUE,
                                   legend.title = list(gpar(fontsize = 14)),
                                   legend.text = list(gpar(fontsize = 12))))
    ggsave(paste(filter, "/notord_heatmap_",plot_name,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)

#all done! Yay!
