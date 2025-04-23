library(umap)
library(Rtsne)
library(data.table)
library(ggplot2)
library(pheatmap)
library(ggplotify)
library(grid)

output_dir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/basic_pipeline"
plots_dir <- "plots_limited"
setwd(output_dir)
dir.create(plots_dir)
setwd(plots_dir)

num_cutoff <- 30 #for now, selecting top 30 most variable genes

####### need to include code to import the data. 
    #For now, it's imported in GeoMx_custom.R and still loaded.

#find the most highly variable genes first
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- apply(log2(impute_zero), MARGIN = 1, FUN = calc_CV)
top <- sort(CV_dat, decreasing = TRUE)[1:num_cutoff]

genes_to_keep <- names(top)

impute_zero <- impute_zero[genes_to_keep, ]

#start UMAP
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

#####by segment
umap_out <-
    umap(t(log(impute_zero)),
         config = custom_umap)
metadata[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
UMAP_segment <- ggplot(metadata,
       aes(x = UMAP1, y = UMAP2, color = Custom_CellType)) +
    geom_point(size = 3) + 
    stat_ellipse(data = subset(metadata, Custom_CellType %in% c("Keratinocyte", "Dermal Stroma")),
                 aes(fill = Custom_CellType), alpha = 0.2, geom = "polygon") +
    theme_bw()

ggsave("UMAP_celltype.png",
       UMAP_segment, width = 10, height = 8, dpi = 300)

#####by diagnosis
UMAP_diagnosis <- ggplot(metadata,
       aes(x = UMAP1, y = UMAP2, color = Custom_Diagnosis)) +
    geom_point(size = 3) + 
    theme_bw()

ggsave("UMAP_diagnosis.png",
       UMAP_diagnosis, width = 10, height = 8, dpi = 300)

########################run tSNE################################
set.seed(42) # set the seed for tSNE as well

#####by cell type
tsne_out <-
    Rtsne(t(log2(impute_zero)),
          perplexity = ncol(impute_zero)*.15)
metadata[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
tSNE_segment <- ggplot(metadata, aes(x = tSNE1, y = tSNE2, color = Custom_CellType)) +
    geom_point(size = 3, alpha = 0.6) +  # Slight transparency for clarity
    stat_ellipse(data = subset(metadata, Custom_CellType %in% c("Keratinocyte", "Dermal Stroma")),
                 aes(fill = Custom_CellType), alpha = 0.2, geom = "polygon") + 
    theme_bw()

ggsave("tSNE_celltype.png",
       tSNE_segment, width = 10, height = 8, dpi = 300)

#####by diagnosis
tsne_out <-
    Rtsne(t(log2(impute_zero)),
          perplexity = ncol(impute_zero)*.15)
metadata[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
tSNE_diagnosis <- ggplot(metadata, aes(x = tSNE1, y = tSNE2, color = Custom_Diagnosis)) +
    geom_point(size = 3, alpha = 0.6) +  # Slight transparency for clarity 
    theme_bw()

ggsave("tSNE_diagnosis.png",
       tSNE_diagnosis, width = 10, height = 8, dpi = 300)

#####other option: clustering high CV genes #####################################

#####for segment
annotation <- metadata[, c("Sample_ID", "Custom_CellType"), drop = FALSE]
annotation <- annotation[order(annotation$Custom_CellType), ]
rownames(annotation) <- annotation[, "Sample_ID"]
annotation[, "Sample_ID"] <- NULL

impute_zero <- impute_zero[ , rownames(annotation)]

calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- apply(log2(impute_zero), MARGIN = 1, FUN = calc_CV)
sort(CV_dat, decreasing = TRUE)[1:30]

# Identify genes in the top 3rd of the CV values
#GOI <- names(CV_dat)[CV_dat >= quantile(CV_dat, 0.9)]
heatmap_plot <- as.ggplot(pheatmap(log2(impute_zero),
                                   scale = "row",
                                   show_rownames = FALSE, show_colnames = FALSE,
                                   border_color = NA,
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   cluster_cols = FALSE,
                                   breaks = seq(-3, 3, 0.05),
                                   bg = "white",
                                   color = colorRampPalette(c("#0a71cb", "#FFFFFF", "#D1495B"))(120),
                                   annotation_col = annotation))

ggsave("celltype_heatmap_all.png",
       heatmap_plot, width = 10, height = 8, dpi = 300)

heatmap_plot <- as.ggplot(pheatmap(log2(impute_zero),
                                   scale = "row",
                                   show_rownames = FALSE, show_colnames = FALSE,
                                   border_color = NA,
                                   clustering_method = "average",
                                   clustering_distance_rows = "correlation",
                                   cluster_cols = TRUE,
                                   breaks = seq(-3, 3, 0.05),
                                   bg = "white",
                                   color = colorRampPalette(c("#0a71cb", "#FFFFFF", "#D1495B"))(120),
                                   annotation_col = annotation))

ggsave("notord_celltype_heatmap_all.png",
       heatmap_plot, width = 10, height = 8, dpi = 300)

#####by diagnosis
annotation <- metadata[, c("Sample_ID", "Custom_Diagnosis"), drop = FALSE]
annotation <- annotation[order(annotation$Custom_Diagnosis), ]
rownames(annotation) <- annotation[, "Sample_ID"]
annotation[, "Sample_ID"] <- NULL

impute_zero <- impute_zero[ , rownames(annotation)]

# Identify genes in the top 3rd of the CV values
heatmap_plot <- as.ggplot(pheatmap(log2(impute_zero),
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

ggsave("diagnosis_heatmap.png",
       heatmap_plot, width = 10, height = 8, dpi = 300)

heatmap_plot <- as.ggplot(pheatmap(log2(impute_zero),
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

ggsave("notord_diagnosis_heatmap.png",
       heatmap_plot, width = 10, height = 8, dpi = 300)

#####subset to one cell type#####################################

UMAP_tSNE_heatmap <- function(metadata, counts, factor, filter) {

    dir.create(filter)

    #########UMAP##########
    UMAP <- ggplot(metadata,
       aes(x = UMAP1, y = UMAP2, color = metadata[[factor]])) +
    geom_point(size = 3) + 
    theme_bw()

    ggsave(paste(filter,"/UMAP_",factor,".png", sep = ""),
       UMAP, width = 10, height = 8, dpi = 300)

    #########t-SNE##########
    tsne_out <-
       Rtsne(t(log2(counts)),
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

    heatmap_plot <- as.ggplot(pheatmap(log2(counts),
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

##### dermal stroma ######
metadata_ds <- subset(metadata, metadata$Custom_CellType == "Dermal Stroma")
counts_ds <- impute_zero[, metadata_ds$Sample_ID]
UMAP_tSNE_heatmap(metadata_ds, counts_ds, "Custom_Diagnosis", "Dermal_Stroma")

##### keratinocytes ######
metadata_kt <- subset(metadata, metadata$Custom_CellType == "Keratinocyte")
counts_kt <- impute_zero[, metadata_kt$Sample_ID]
UMAP_tSNE_heatmap(metadata_kt, counts_kt, "Custom_Diagnosis", "Keratinocyte")
