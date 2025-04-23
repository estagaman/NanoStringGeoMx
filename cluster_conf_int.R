library(umap)
library(Rtsne)
library(data.table)
library(ggplot2)
library(pheatmap)
library(ggplotify)
library(grid)

counts <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/basic_pipeline/QuantNorm_geneQC_10.csv"
output_dir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/basic_pipeline/plots_conf_int"
metadata <- as.data.frame(fread("/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/metadata/skin_layer/all_samples.csv"))

counts <- as.data.frame(fread(counts))

dir.create(output_dir)
setwd(output_dir)

#clean the data
rownames(counts) <- counts$V1
counts$V1 <- NULL
dcc_kept <- sub("\\.dcc$", "", colnames(counts))

#subset the metadata
dcc_kept <- sub("\\.dcc$", "", colnames(counts))
metadata <- subset(metadata, metadata$Sample_ID %in% dcc_kept)
metadata <- metadata[match(dcc_kept, metadata$Sample_ID), ]
metadata$Sample_ID == dcc_kept
metadata$V1 <- NULL
colnames(counts) <- dcc_kept

################# plotting ##############################

UMAP_tSNE_heatmap_2groups <- function(metadata, counts, factor, filter, groups, plot_name) {

    dir.create(filter)

    custom_umap <- umap::umap.defaults
    custom_umap$random_state <- 42
    umap_out <-
        umap(t(log(counts)),
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
       Rtsne(t(log2(counts)),
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
    annotation <- metadata[, c("Sample_ID", factor), drop = FALSE]
    annotation <- annotation[order(annotation[, factor]), ]
    rownames(annotation) <- annotation[, "Sample_ID"]
    annotation[, "Sample_ID"] <- NULL
    counts <- counts[ , rownames(annotation)]

    calc_CV <- function(x) {sd(x) / mean(x)}
    CV_dat <- apply(log2(counts), MARGIN = 1, FUN = calc_CV)
    sort(CV_dat, decreasing = TRUE)[1:30]

    GOI <- names(CV_dat)[CV_dat >= quantile(CV_dat, 0.9)]
    counts <- counts[GOI, ]

    md <- subset(metadata, metadata$filter %in% groups)
    counts <- counts[, colnames(counts) %in% md$Sample_ID]

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
    ggsave(paste(filter, "/ord_heatmap_",plot_name,".png", sep = ""),
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
    ggsave(paste(filter, "/notord_heatmap_",plot_name,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)
}

##### dermal stroma ######
metadata_ds <- subset(metadata, metadata$Custom_CellType == "Dermal Stroma")
counts_ds <- counts[, metadata_ds$Sample_ID]
groups <- c("CTCL", "HC", "AD")
plot_name <- "CTCL_HC_AD"
UMAP_tSNE_heatmap_2groups(metadata_ds, counts_ds, "Custom_Diagnosis", "Dermal_Stroma", groups, plot_name)

groups <- c("CTCL", "HC")
plot_name <- "CTCL_HC"
UMAP_tSNE_heatmap_2groups(metadata_ds, counts_ds, "Custom_Diagnosis", "Dermal_Stroma", groups, plot_name)


##### keratinocytes ######
metadata_kt <- subset(metadata, metadata$Custom_CellType == "Keratinocyte")
counts_kt <- counts[, metadata_kt$Sample_ID]
groups <- c("DR", "CTCL")
plot_name <- "DR_CTCL"
UMAP_tSNE_heatmap_2groups(metadata_kt, counts_kt, "Custom_Diagnosis", "Keratinocyte", groups, plot_name)

groups <- c("AD", "CTCL")
plot_name <- "AD_CTCL"
UMAP_tSNE_heatmap_2groups(metadata_kt, counts_kt, "Custom_Diagnosis", "Keratinocyte", groups, plot_name)

groups <- c("PSO", "CTCL")
plot_name <- "PSO_CTCL"
UMAP_tSNE_heatmap_2groups(metadata_kt, counts_kt, "Custom_Diagnosis", "Keratinocyte", groups, plot_name)



