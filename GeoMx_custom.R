library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

#name directory I would like output expression files to save to
output_dir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/basic_pipeline"
plots_dir <- "plots"

datadir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/data"
DCCFiles <- list.files(datadir, pattern = "\\.dcc$", 
                       full.names = TRUE, recursive = FALSE)
PKCFiles <- list.files(datadir, pattern = "\\.pkc$", 
                       full.names = TRUE, recursive = FALSE)
PKC_data <- readPKCFile(PKCFiles, default_pkc_vers=NULL)
SampleAnnotationFile <- list.files(datadir, pattern = "\\.xlsx$", 
                       full.names = TRUE, recursive = FALSE)

#Create GeoMx set object - yay!
demoData <- readNanoStringGeoMxSet(dccFiles = DCCFiles, #count/quality data
    pkcFiles = PKCFiles, #path to probes, gene targets used
    phenoDataFile = SampleAnnotationFile, #excel: regions/cell types
    phenoDataSheet = "Sheet1",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("Aoi", "Roi"),
    experimentDataColNames = c("Panel"))

#check the annotation file was loaded
library(knitr)
pkcs <- annotation(demoData) #use annotation() to extract tissue info
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules)) #Yay!

#Shift any zero counts to one
demoData <- shiftCountsOne(demoData, useDALogic = TRUE) #shifting the expression counts, the rna counts 

######check QC########
QC_params <-
    list(minSegmentReads = 1000, # Minimum number of reads (1000)
         percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
         percentStitched = 80,   # Minimum % of reads stitched (80%)
         percentAligned = 80,    # Minimum % of reads aligned (80%) - Zi used 80
         percentSaturation = 50, # Minimum sequencing saturation (50%)
         minNegativeCount = 2,   # Minimum negative control counts (10)
         maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
         minNuclei = 100,         # Minimum # of nuclei estimated (100)
         minArea = 5000        # Minimum segment area (5000)
         )
demoData <-
    setSegmentQCFlags(demoData, 
                      qcCutoffs = QC_params)   #actually do the filtering     
# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
    ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
    c(sum(QCResults[, "QCStatus"] == "PASS"),
      sum(QCResults[, "QCStatus"] == "WARNING"))

#to check: 
QC_Summary["TOTAL FLAGS", ] #subsetting down to 113 based on this

#make me a subset of samples that pass this neg control QC
QC_passOnly <- subset(QCResults, QCResults$QCStatus == "PASS")
samp_pass <- rownames(QC_passOnly)


##########QC done: samples unedited#########

##########NegGeoMean stuff#################
negativeGeoMeans <- 
    esBy(negativeControlSubset(demoData), 
         GROUP = "Module", 
         FUN = function(x) { 
             assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

#save NegGeoMeans as csv file
write.csv(protocolData(demoData)[["NegGeoMean"]], "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/NegGeoMean/NegGeoMean.csv")

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
for(ann in negCols) {
    plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
    print(plt)
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(demoData)$NegGeoMean),
      col.names = c("NTC Count", "# of Segments"))

#REMOVE SEGMENTS NOT PASSING FILTERS
demoData <- demoData[, QCResults$QCStatus == "PASS"]

#################### do PROBE-level QC ################
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)

exprs(target_demoData)[1:5, 1:2]

#subset to group that I preserved for NegGeoMean threshold
library('data.table')
kept_nGeoMean <- fread("/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/NegGeoMean/CANX_ACTBL2_PPIA.csv")
kept <- colnames(kept_nGeoMean)

target_demoData_sub <- target_demoData[, colnames(target_demoData) %in% kept]

t_f <- data.frame(sample_name = colnames(target_demoData), in_kept = colnames(target_demoData) %in% kept)
t_f <- subset(t_f, t_f$in_kept == FALSE)

t_f_kept <- data.frame(sample_name = kept, in_data = kept %in% colnames(target_demoData))
#filter down samples to ones containing housekeeping genes at high levels
cutoff <- 1
minLOQ <- 1

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                   module)
    if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
        LOQ[, module] <-
            pmax(minLOQ,
                 pData(target_demoData)[, vars[1]] * 
                     pData(target_demoData)[, vars[2]] ^ cutoff)
    }
}
pData(target_demoData)$LOQ <- LOQ

LOQ_Mat <- c()
for(module in modules) {
    ind <- fData(target_demoData)$Module == module
    Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                       FUN = function(x) {
                           x > LOQ[, module]
                       }))
    LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]

pData(target_demoData)$GenesDetected <- 
    colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
    pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
    cut(pData(target_demoData)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = Segment)) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment")

#can also review this as a table 
kable(table(pData(target_demoData)$DetectionThreshold,
            pData(target_demoData)$Segment))

#choose your filtering threshold: use the graph to do this depending on YOUR results 
target_demoData <-
    target_demoData[, pData(target_demoData)$GeneDetectionRate >= .05]

dim(target_demoData) #down to 110 samples

#what about a more strict threshold?
target_demoData_strict <- target_demoData[, pData(target_demoData)$GeneDetectionRate >= .075] #88 samples



#GENE DETECTION RATE across the study 
    #creating a gene list (goi) to review

library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
    fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

#calculating percent of segments in which a certain gene is detected 

# Gene of interest detection table: CHOOSE FOR YOUR STUDY WHAT IS INTERESTING
goi <- c("LACTB2", "LACTBL1", "ACTBL2", "LACTB", "TBP", 
"CANX", "B2M", "YWHAZ", "SDHA", "CYC1",
 "IPO8", "PUM1", "MRPL19", "PSMC4", "TFRC",
  "HPRT1", "GAPDHS", "PGK1", "PPIA", "TUBB", "APRT", "POLR2A", "EEF1A2") #choose the specific genes that are most of interest to you 
goi_df <- data.frame(
    Gene = goi,
    Number = fData(target_demoData)[goi, "DetectedSegments"],
    DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))

#GENE FILTERING
    #graph the total number of genes detected in different percentages of segments 
    #select how many low detected genes to filter out of the dataset 
    #gene filtering overall increases performance downstream 

plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
              vjust = 1.6, color = "black", size = 4) +
    scale_fill_gradient2(low = "orange2", mid = "lightblue",
                         high = "dodgerblue3", midpoint = 0.65,
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "% of Segments",
         y = "Genes Detected, % of Panel > LOQ")

#let's do a 5% filter and see what happens
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_demoData_detected <- 
    target_demoData[fData(target_demoData)$DetectionRate >= 0.05 |
                        fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData_detected)
#12788 features

#ok what about a 10% filter
target_demoData_detected_10 <- 
    target_demoData[fData(target_demoData)$DetectionRate >= 0.1 |
                        fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData_detected_10)
#7216 features

#look at which housekeeping genes were preserved
goi_05 <- goi %in% rownames(target_demoData_detected)
goi_10 <- goi %in% rownames(target_demoData_detected_10)

#all the ones I need are kept either way so it's fine
goi_kept_df <- data.frame(gene = goi, goi_05 = goi_05, goi_10 = goi_10)

####STOP: save some intermediate files
setwd(output_dir)

write.csv(exprs(target_demoData), "raw_counts_sampQC.csv")
write.csv(exprs(target_demoData_detected), "raw_counts_geneQC_05.csv")
write.csv(exprs(target_demoData_detected_10), "raw_counts_geneQC_10.csv")

############ 3rd quartile normalization ########################################
library(reshape2)  # for melt
library(cowplot)   # for plot_grid

target_demo <- target_demoData_detected_10 #set to the data you would like to use for normalization
    #depends on gene detection rate you subsetted to

# Graph Q3 value vs negGeoMean of Negatives: Just checking for separation 
#if there are segments with very little separation, I should remove them
ann_of_interest <- "Segment"
Stat_data <- 
    data.frame(row.names = colnames(exprs(target_demo)),
               Segment = colnames(exprs(target_demo)),
               Annotation = pData(target_demo)[, ann_of_interest],
               Q3 = unlist(apply(exprs(target_demo), 2,
                                 quantile, 0.75, na.rm = TRUE)),
               NegProbe = exprs(target_demo)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
    geom_histogram(bins = 40) + theme_bw() +
    scale_x_continuous(trans = "log2") +
    facet_wrap(~Annotation, nrow = 1) + 
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point() + guides(color = "none") + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                    rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

### perform the Q3 normalization
target_demo <- normalize(target_demo ,
                    norm_method = "quant", 
                    desiredQuantile = .75,
                    toElt = "q_norm")

#visualize what has changed
boxplot(exprs(target_demo)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

boxplot(assayDataElement(target_demo[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")

#now I need to save these results and check for housekeeping gene coexpression among my goi's
write.csv(assayDataElement(target_demo, elt = "q_norm"), "Q3norm_geneQC_10.csv")

###quantile normalization instead
library("preprocessCore")
norm.quantile <- normalize.quantiles(as.matrix(exprs(target_demo)))
dimnames(norm.quantile) <- dimnames(exprs(target_demo))

#save results of the quantile normalization as well
write.csv(norm.quantile, "QuantNorm_geneQC_10.csv")

##### try the coefficients of variance again
#####for segment
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- apply(log2(norm.quantile), MARGIN = 1, FUN = calc_CV)
sort(CV_dat, decreasing = TRUE)[(length(CV_dat)- 30):length(CV_dat)]

#####try geNorm on the normalized data to see which housekeeping genes pop out
library("ctrlGene")
count_matrix_limited <- norm.quantile[rownames(norm.quantile) %in% goi, ]

ge_results <- geNorm(t(count_matrix_limited), genes = data.frame(Genes = character(0), Avg.M =
  numeric(0)), ctVal = TRUE)

View(ge_results)

#try subsetting samples to those where housekeeping genes are above 1 in the original counts 
library("data.table")
norm_counts <- as.data.frame(fread("/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/imputation/samples_TFRC_82.csv"))
samples_kept <- colnames(norm_counts)[2:ncol(norm_counts)]

count_matrix_limited <- count_matrix_limited[, colnames(count_matrix_limited) %in% samples_kept]

ge_results2 <- geNorm(t(count_matrix_limited), genes = data.frame(Genes = character(0), Avg.M =
  numeric(0)), ctVal = TRUE)

View(ge_results2)

#alright, just try plotting the plots I want to plot and see what happens
############ PLOTTING YAY ###############

# Prep the data
impute_zero <- norm.quantile
metadata <- as.data.frame(fread("/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/metadata/skin_layer/all_samples.csv"))

dir.create(plots_dir)
setwd(plots_dir)

#subset the metadata - clean the counts data
dcc_kept <- sub("\\.dcc$", "", colnames(impute_zero))
metadata <- subset(metadata, metadata$Sample_ID %in% dcc_kept)
metadata <- metadata[match(dcc_kept, metadata$Sample_ID), ]
metadata$Sample_ID == dcc_kept
metadata$V1 <- NULL
colnames(impute_zero) <- dcc_kept

########################run UMAP################################
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
GOI <- names(CV_dat)[CV_dat >= quantile(CV_dat, 0.9)]
heatmap_plot <- as.ggplot(pheatmap(log2(impute_zero[GOI,]),
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

heatmap_plot <- as.ggplot(pheatmap(log2(impute_zero[GOI,]),
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
heatmap_plot <- as.ggplot(pheatmap(log2(impute_zero[GOI, ]),
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

heatmap_plot <- as.ggplot(pheatmap(log2(impute_zero[GOI, ]),
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

    heatmap_plot <- as.ggplot(pheatmap(log2(counts[GOI, ]),
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

    heatmap_plot <- as.ggplot(pheatmap(log2(counts[GOI, ]),
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

    calc_CV <- function(x) {sd(x) / mean(x)}
    CV_dat <- apply(log2(counts), MARGIN = 1, FUN = calc_CV)
    lower <- quantile(CV_dat, 0.8)
    upper <- quantile(CV_dat, 0.9)

    GOI_adj <- names(CV_dat)[CV_dat > lower & CV_dat <= upper]
    heatmap_plot <- as.ggplot(pheatmap(log2(counts[GOI, ]),
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
    ggsave(paste(filter, "/notord_heatmap_80_90_",factor,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)

    lower <- quantile(CV_dat, 0.7)
    upper <- quantile(CV_dat, 0.8)

    GOI_adj <- names(CV_dat)[CV_dat > lower & CV_dat <= upper]
    heatmap_plot <- as.ggplot(pheatmap(log2(counts[GOI, ]),
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
    ggsave(paste(filter, "/notord_heatmap_70_80_",factor,".png", sep = ""),
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

heatmap_adjust_quantiles <- function(metadata, counts, factor, out_dir){
    annotation <- metadata[, c("Sample_ID", factor), drop = FALSE]
    annotation <- annotation[order(annotation[, factor]), ]
    rownames(annotation) <- annotation[, "Sample_ID"]
    annotation[, "Sample_ID"] <- NULL
    counts <- counts[ , rownames(annotation)]


    calc_CV <- function(x) {sd(x) / mean(x)}
    CV_dat <- apply(log2(counts), MARGIN = 1, FUN = calc_CV)
    lower <- quantile(CV_dat, 0.8)
    upper <- quantile(CV_dat, 0.9)

    GOI_adj <- names(CV_dat)[CV_dat > lower & CV_dat <= upper]
    heatmap_plot <- as.ggplot(pheatmap(log2(counts[GOI, ]),
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
    ggsave(paste(filter, "/notord_heatmap_80_90_",factor,".png", sep = ""),
       heatmap_plot, width = 10, height = 8, dpi = 300)
}

out_dir_80_90 <- "quant_80_90"
dir.create(out_dir_80_90)

heatmap_adjust_quantiles(metadata, impute_zero, "Custom_CellType", out_dir_80_90)
heatmap_adjust_quantiles(metadata, impute_zero, "Custom_Diagnosis", out_dir_80_90)
