#QC_pipeline.R: pipeline for sample and gene-level quality control of raw DCC files from Nanostring GeoMx sequencing

#INPUT: 
    #output_dir: directory you would like output to save to
    #datadir: path to folder where DCC files, PKC files, and metadata (as excel .xlsx) are saved to
        #metadata should contain column "Sample_ID" which contains the name of the matching DCC file
    
    #along the script, there are options to change the QC parameters as you see fit for your data. 
        #the parameters input here are what I chose to use for my analysis. 

#load libraries
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

#name directory I would like output files to save to
output_dir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/final/NG1_gene15"

#name file that DCC files, PKC files, and metadata are stored to
datadir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/data"

dir.create(output_dir)
plots_dir <- "plots"

#load DCC files
DCCFiles <- list.files(datadir, pattern = "\\.dcc$", 
                       full.names = TRUE, recursive = FALSE)

#load PKC files                       
PKCFiles <- list.files(datadir, pattern = "\\.pkc$", 
                       full.names = TRUE, recursive = FALSE)
PKC_data <- readPKCFile(PKCFiles, default_pkc_vers=NULL)

#load metadata
SampleAnnotationFile <- list.files(datadir, pattern = "\\.xlsx$", 
                       full.names = TRUE, recursive = FALSE)

#Create GeoMx set object - yay!
demoData <- readNanoStringGeoMxSet(dccFiles = DCCFiles, #count/quality data
    pkcFiles = PKCFiles, #path to probes, gene targets used
    phenoDataFile = SampleAnnotationFile, #excel: regions/cell types
    phenoDataSheet = "Sheet1",
    phenoDataDccColName = "Sample_ID", #column in metadata that matches the DCC file name
    protocolDataColNames = c("Aoi", "Roi"),
    experimentDataColNames = c("Panel"))

#check the annotation file was loaded
library(knitr)
pkcs <- annotation(demoData) #use annotation() to extract tissue info
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules)) #Yay!

#Shift any zero counts to one
demoData <- shiftCountsOne(demoData, useDALogic = TRUE) #shifting the expression counts, the rna counts 

#extract the original list of features before any filtering
features_og <- pData(featureData(demoData))

#extract just the negative counts for inspection
save_neg_counts <- subset(features_og, features_og$CodeClass == "Negative")
neg_counts <- exprs(demoData)[rownames(save_neg_counts), ]

##inspecting total reads per segment
counts_matrix <- exprs(demoData)
colSums(counts_matrix)

######check QC########
QC_params <-
    list(minSegmentReads = 1000, # Minimum number of reads (1000)
         percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
         percentStitched = 80,   # Minimum % of reads stitched (80%)
         percentAligned = 80,    # Minimum % of reads aligned (80%) - Zi used 80
         percentSaturation = 50, # Minimum sequencing saturation (50%)
         minNegativeCount = 2,   # Minimum negative control counts (10)
         maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
         minNuclei = 100,         # Minimum # of nuclei estimated (100)
         minArea = 5000        # Minimum segment area (5000)
         )

#set QC flags according to the above parameters
demoData <-
    setSegmentQCFlags(demoData, 
                      qcCutoffs = QC_params)

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

#to check the number of ROIs we will keep after filtering: 
QC_Summary["TOTAL FLAGS", ]

#check Negative Control counts
negativeGeoMeans <- 
    esBy(negativeControlSubset(demoData), 
         GROUP = "Module", 
         FUN = function(x) { 
             assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

#save NegGeoMeans as csv file
write.csv(protocolData(demoData)[["NegGeoMean"]], paste0(output_dir, "/NegGeoMean.csv"))

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

#save the raw counts
write.csv(exprs(target_demoData), paste0(output_dir, "/raw_counts.csv"))

##now, we do gene-level QC
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
    target_demoData[, pData(target_demoData)$GeneDetectionRate >= .2]

dim(target_demoData)
    #102 samples using negative control 2 and gene detection 0.05
    #294 samples using negative control 1 and gene detection 0.05
    #226 samples using negative control 1 and gene detection 0.1
    #126 samples using negative control 1 and gene detection 0.15

#GENE DETECTION RATE across the study
    #creating a gene list (goi) to review

library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
    fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

#calculating percent of segments in which a certain gene is detected 

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

gene_detect_plot <- ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
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

ggsave(paste0(output_dir, "/gene_detect_plot.png"), gene_detect_plot)

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
#9567 features using NG2 cutoff
#13206 features using NG1 cutoff
#15364 features using NG1 and segment QC 0.1
#17145 features using NG1 and segment QC 0.15

####STOP: save some intermediate files
write.csv(exprs(target_demoData_detected_10), paste0(output_dir, "/raw_counts_geneQC_10.csv"))

#create counts file and feature annotation file with the negative probes separate
all_exprs <- exprs(target_demoData_detected_10)
feature_data <- pData(featureData(target_demoData_detected_10))

#remove the aggregated negative probe
features_noNeg <- subset(feature_data, feature_data$CodeClass == "Endogenous")
exprs_noNeg <- all_exprs[rownames(features_noNeg), ]

#add on the unaggregated negative probes
save_neg_counts$DetectedSegments <- c(rep("neg", nrow(save_neg_counts)))
save_neg_counts$DetectionRate <- c(rep("neg", nrow(save_neg_counts)))
save_neg_counts$RTS_ID <- NULL
save_neg_counts$ProbeID <- NULL
features_add_neg <- rbind(features_noNeg, save_neg_counts) #shared columns: 

neg_counts_sub <- neg_counts[, colnames(exprs_noNeg)]
exprs_add_neg <- rbind(exprs_noNeg, neg_counts_sub) #works 

#save the version with negative controls
write.csv(features_add_neg, paste0(output_dir, "/features_w_neg.csv"))
write.csv(exprs_add_neg, paste0(output_dir, "/counts_w_neg.csv"))

#save the probe-level information to plug into DEG analysis
write.csv(pData(featureData(target_demoData_detected_10)), paste0(output_dir, "/feature_data.csv"))

#done! ready for normalization, visualization, and DEG analysis
