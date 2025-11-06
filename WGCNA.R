#WGCNA.R: weighted gene correlation network analysis, testing association of disease/diagnosis with gene modules

#INPUT: 
  #counts file with counts data 
  #metadata file with metadata
  #output_dir you want results to save to
  #EdgeR_filt: whether you want to perform EdgeR filtering
  #immune_genes: whether you want to filter down to a set of relevant immune genes

#OUTPUT: 
  #gene_modules.txt: all the genes and their assigned module
  #sorted_gene_modules.txt: all the genes and their assigned module, sorted by module
  #gee_results.csv: p-values and beta for the association of each module with a particular diagnosis
  #dendrogram.png: tree showing relationships between each module of genes

#Methods Applied: 
  #EdgeR filtering
  #filtering down to immune genes
  #gene correlation analysis 
  #GEE for testing association between module eigengenes and traits 

#load all necessary packages
library("WGCNA")
library("tidyverse")
library("edgeR")
library("limma")
library("readxl")
library("data.table")
library("glmtoolbox")
library("geesmv")


countFile <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/unfiltered_analysis/batch_correction/batch_effect/quant_06_16_25/counts_w_neg.csv"
output_dir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/final/NG2_11_04_25/WGCNA_analysis"
metadataFile <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/metadata/skin_layer/all_samples.csv"

#specify whether I want to add voom normalization and/or quantile normalization
edgeR_filt <- TRUE
immune_genes <- TRUE
filter = "Keratinocyte"

#change to the output directory, create it if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)}
setwd(output_dir)

countFile <- read.csv(countFile) %>% as.data.frame()
metadata <- fread(metadataFile) %>% as.data.frame()

#clean the counts data up
counts <- countFile
rownames(counts) <- counts$X
counts$X <- NULL
dcc_kept <- sub("\\.dcc$", "", colnames(counts))
dcc_kept <- gsub("\\.", "-", dcc_kept)
colnames(counts) <- dcc_kept

#subset the metadata down to the ROIs kept after QC
metadata <- subset(metadata, metadata$Sample_ID %in% dcc_kept)
metadata <- metadata[match(dcc_kept, metadata$Sample_ID), ]
metadata$Sample_ID == dcc_kept
metadata$V1 <- NULL

#remove the samples that do not have custom diagnosis assigned 
metadata <- subset(metadata, !(is.na(Custom_Diagnosis)))
counts <- counts[, metadata$Sample_ID]

#remove negative controls from the counts data
counts <- counts[!grepl("^RTS", rownames(counts)), ]

probes_save <- rownames(counts)
samples_save <- colnames(counts)

counts_save <- counts
metadata_save <- metadata

#apply filtering by cell type
filter = "Keratinocyte"
metadata <- subset(metadata, Custom_CellType == filter)
counts <- counts[ ,metadata$Sample_ID]

#save the names of the probes and samples for later
probes <- rownames(counts)
samples <- colnames(counts)

#create design matrix based on diagnosis
design <- model.matrix(~0 + Custom_Diagnosis, data = metadata)

#performing EdgeR filtering
keep <- filterByExpr(counts, design, min.total.count = 15, min.count = 10)
counts <- counts[keep, ]

#convert counts to DGE List
dge <- DGEList(counts = counts)

#perform TMM normalization
dge_tmm <- calcNormFactors(dge, method = "TMM")

#perform voom normalization with TMM factors and log cpm scaling
v <- voom(dge_tmm, design, plot = FALSE)

#extract the normalized counts
new_counts <- v$E

#filter down to only immune-related genes or genes from the scRNAseq paper
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

    #which ones are in my list of genes 

    all_genes_immune <- c(symbols_immune, genes_imm, septic_genes)
    all_genes_immune <- toupper(all_genes_immune)
    all_genes_immune <- unique(all_genes_immune)

    my_genes_immune <- intersect(genes_to_annotate, all_genes_immune)
    length(my_genes_immune) #down to 2451 genes that are in the innate immunity, immunity, and septic shock databases 

    cell_marker_genes <- c("NRXN1", "MLANA", "LYVE1", "GJA4", "ACKR1", "KRT19", "KRT18", "DCD", "KRT14", "KRT5", "KRT10", "KRT2", "ACTA2", "DCN", "COL1A1", "CD163", "CLEC9A", "LAMP3", "ITGAX", "LILRA4", "CD79A", "MKI67", "FOXP3", "KLRD1", "KLRB1", "CCL5", "CD8A", "CXCL13", "IL22", "IL13", "CD4", "CD3D")
    #fibroblast markers: DCN, COL1A1
    #keratinocyte markers: KRT5, KRT10, KRT2, ACTA2

    DE_genes_article <- c("CD3D", "CD4", "CD8A", "IL13", "TXN", "CAPG", "LY6D", "CSF2", "JAML", "CD96", "CD9", "CD63", "ITGAE", "IL26", "EGLN3", "IL22", "LINC00892", "IL9R", "NCF4", "MIIP", "CCL4", "NKG7", "CCL5", "GZMH", "CCL4L2", "GZMA", "INFG", "CXCR4", "CLEC2B", "ARL4C", "THEMIS", "PLAAT4", "CXCL13", "CCR7", "SELL", "PGM2L1", "PASK", "TCF7", "KIR3DL2", "CD27", "TSHZ2", "SESN3", "TTC39C", "LTB", "COX5A", "CORO1B", "NDFIP1", "C1QBP", "PIM2", "ID3", "CTLA4", "LINC01480", "RCSD1")
    Tcell_genes <- c("CD4", "CD52", "IL7R", "CXCL13", "IL22", "IL4", "IL13", "TNFRSF4", "TNFRSF18", "KLRB1", "KLRD1", "CD8A", "INFG", "CCL5", "NKG7", "CCL4", "NKG7", "CCL4", "GZMA", "GZMK", "CTSW", "FOXP3", "CTLA4", "MKI67", "IL17A", "IL17F")
    FB_genes <- c("CCL19", "APOD", "APOE", "C3", "MFAP5", "CCN5", "FBLN1", "FBN1", "MMP2", "COL18A1", "APCDD1", "COMP", "ASPN", "COL11A", "TNN", "COCH", "MKI67", "ACTA2", "TAGLN", "RGS5", "CCL26", "COL6A6", "COL6A5", "POSTN", "TNC", "C1QTNF1", "TGFBI")
    KC_genes <- c("KRT6A", "KRT6B", "KRT6C", "KRT16", "KRT17", "S100A2", "S100A7", "S100A9", "SERPINB3", "SERPINB4", "FABP5", "SPRR1B", "IFI27", "KRT1", "KRT2", "KRT10", "DMKN", "KRTDAP", "CXCL9", "CXCL10", "CXCL11", "GBP4", "KRT5", "KRT14", "KRT15", "CHI3L1", "MGST1", "KRT23", "MKI67", "DCD", "PIP", "KRT18", "KRT19", "SNORC", "GBP1", "HLA-DRA", "CD74", "GLUL", "ASPN", "NMU")
    non_T_cell <- c("CD79A", "JCHAIN", "LILRA4", "CD207", "CD1A", "ITGAX", "IL1B", "CLEC10A", "LAMP3", "CCR7", "CCL22", "XCR1", "CLEC9A", "MRC1", "CD163", "CD68", "MAF")
    skin_top_clone <- c("CD4", "CD8A", "IL13", "AREG", "ANXA1", "IL22", "ETS2", "SOCS1", "AHR", "RAB11FIP1", "PLA2G16", "TGIF1", "CLU", "IL17A", "IL17F", "CCL20", "CTSH", "PTMS", "GPR15", "MGAT4A", "GBP5", "GYG1", "IL26", "F2R", "CEBPD", "CD7", "LAG3", "CLEC2B", "IFNG", "NKG7", "CCL4", "CCL5", "GZMH", "GZMB", "CCL4L2", "GZMA", "PRF1", "KLRK1", "SELL", "PGM2L1", "CCR7", "TCF7", "PASK", "ID3", "TSHZ2", "KIR3DL2", "LTB", "HACD1", "CD27", "TNFRSF4", "TIAM1", "IGFL2", "PIM2", "NINJ2", "NME1")

    all_together <- c(DE_genes_article, Tcell_genes, FB_genes, KC_genes, non_T_cell, skin_top_clone) %>% unique()

    immune_plus_paper <- c(all_together, my_genes_immune) %>% unique()

    #take only the counts, normalized counts, and weights for these designated genes
    new_counts <- new_counts[rownames(new_counts) %in% immune_plus_paper, ]
}

#make the normalized counts the official counts table to work with downstream
counts <- new_counts

######### CLEANING AND NORMALIZATION DONE #######################


######### NETWORK ANALYSIS #######################

#transpose the data frame
counts_t <- t(counts)

allowWGCNAThreads()

#testing different powers to see what will work best for these data
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(
  counts_t,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5, 
  corFnc = "bicor"
  )

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

#from here, view the plot, look for a power where we are above the line on the left plot, but still on the curve of the right plot

#my chosen power is 6 for this analysis
picked_power = 6
temp_cor <- cor
cor <- WGCNA::cor  # Force it to use WGCNA cor function (fix a namespace conflict issue)

#grouping genes into modules based on correlation
netwk <- blockwiseModules(counts_t,                

                          checkMissingData = TRUE,
                          corType = "bicor",         

                          # == Adjacency Function ==
                          power = picked_power,               
                          networkType = "signed",

                          # == Tree and Block Options ==
                          deepSplit =4,
                          pamRespectsDendro = F,
                          #detectCutHeight = 0.75,
                          minModuleSize = 5,
                          maxBlockSize = 4000,

                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.10,

                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",

                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor

#look at the modules:

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)

# Plot the dendrogram and the module colors underneath
png("dendrogram.png", width = 1200, height = 800, res = 150)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()
netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)

#make a dataframe of each gene and the module it belongs to:
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

#save this data frame to our output folder
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

######### NETWORK IS CREATED #######################


######### ANALYZING ASSOCIATIONS BETWEEN GENE MODULES AND TRAITS #######################

#First, we start with some visualizations

#isolate module eigengenes
MEs0 <- moduleEigengenes(counts_t, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add diagnoses of each sample as a column
MEs0$treatment = metadata$Custom_Diagnosis

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

#make a heatmap of the module eigengenes and association with diagnosis
module_vs_trait_heat <- mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

#save that plot to my output folder
ggsave(
  filename = "module_trait_relationships.png",
  plot = module_vs_trait_heat,
  width = 8,
  height = 6,
  dpi = 300
)

#Now that we've visualized differences in module-trait relationships, we will perform statistical testing

#extract the eigengenes 
module_eigengenes <- MEs0 

#create dataframe with diagnosis variables as binary 
CTCLvsOther <- ifelse(metadata$Custom_Diagnosis == "CTCL", 1, 0)

PSOvsOther <- ifelse(metadata$Custom_Diagnosis == "PSO", 1, 0)

ADvsOther <- ifelse(metadata$Custom_Diagnosis == "AD", 1, 0)

DRvsOther <- ifelse(metadata$Custom_Diagnosis == "DR", 1, 0)

all_groups <- cbind(CTCLvsOther, PSOvsOther, ADvsOther, DRvsOther) %>% as.data.frame()
rownames(all_groups) <- metadata$Sample_ID

#create a variable for the modules we want to test for association
modules_to_check <- names(module_eigengenes)[1: ncol(module_eigengenes) - 1]

# create a dataframe to store results
gee_results <- data.frame(
  Comparison = character(),
  Module = character(),
  Beta = numeric(),
  Pvalue = numeric(),
  stringsAsFactors = FALSE
)

#create a vector to store failed calculations - GEE can fail if there is not enough variation to calculate
failed_models <- c()

# Loop through each module eigengene
for (module in modules_to_check) {

  # For each contrast: 
  for (comparison in colnames(all_groups)){
        cond1 <- strsplit(comparison, "")[[1]][1:(nchar(comparison)-7)] %>% paste(collapse = "")

        #make a dataframe with the eigengens
        model_data <- data.frame(
            Eigengene = as.numeric(module_eigengenes[[module]]), #module eigengenes
            Trait = relevel(as.factor(all_groups[[comparison]]), ref = "1"), #diagnosis of each sample
            DP = as.factor(metadata$DP) #patient ID of each sample
            )

        #try to fit gee and perform differential expression testing
        tryCatch({
          gee_model <- glmgee(
          formula = Eigengene ~ Trait, #eigengene and response, trait as predictor
          id = DP, #Patient ID as grouping variable
          data = model_data,
          corstr = "Exchangeable",
          family = gaussian()) #gaussian distribution for the eigengenes

          estimate = coef(gee_model)["Trait0", ] #extract the coefficient

          #ok, now adjust the beta and variance with KC estimator
          model_kc = GEE.var.kc(formula = Eigengene ~ Trait, id = "DP", family=gaussian, data = model_data, corstr="exchangeable")

          #estimate of the variance for beta
          beta_var = model_kc$cov.beta[["Trait0"]] #new variance
          std_error = sqrt(beta_var) #calculate std error
          z = estimate/std_error #z-score
          p = 2 * pnorm(-abs(z)) #two-sided p-val

          #save the results
          results_to_keep <- data.frame(Comparison = comparison, Module = module, Beta = estimate, Pvalue = p)

          gee_results <- rbind(gee_results, results_to_keep) # nolint
      
          NULL  # return nothing on success
        }, error = function(e) {
        message("Error fitting model for ", module, comparison, ": ", e$message)
        failed_models <<- c(failed_models, paste(module, comparison))
        NULL
        })
  }
}

#round the p-values to 5 decimal points
gee_results$Pvalue <- round(gee_results$Pvalue, 5)

#save the results to a csv
write.csv(gee_results, "gee_results.csv", row.names = FALSE)

#do ADvsPSO - subset samples down to just AD and PSO, make into binary variable
md_ADvsPSO <- subset(metadata, Custom_Diagnosis %in% c("AD", "PSO"))

ADvsPSO_eigengenes <- module_eigengenes[md_ADvsPSO$Sample_ID, ]
ADvsPSO_counts <- counts_t[md_ADvsPSO$Sample_ID, ]

#make a vector of modules to test
modules_to_check <- colnames(ADvsPSO_eigengenes)[1:ncol(ADvsPSO_eigengenes) - 1]

#make a vector of PSO vs AD
trait_vector = relevel(as.factor(ifelse(md_ADvsPSO$Custom_Diagnosis == "PSO", 1, 0)), ref = "1")

#for each module 
for (module in modules_to_check){

    #the comparison is ADvsPSO
    comparison = "ADvsPSO"

    #cond1 is AD or Trait0 in our vector
    cond1 <- "Trait0"

        #creating our data frame for the model input
        model_data <- data.frame(
            Eigengene = as.numeric(ADvsPSO_eigengenes[[module]]),
            Trait = trait_vector,
            DP = as.factor(md_ADvsPSO$DP)
            )

        #try to fit the model and perform a differential expression test
        tryCatch({
          gee_model <- glmgee(
          formula = Eigengene ~ Trait,
          id = DP,
          data = model_data,
          corstr = "Exchangeable",
          family = gaussian())

          estimate = coef(gee_model)["Trait0", ]

          #ok, now adjust the beta and variance with KC estimator
          model_kc = GEE.var.kc(formula = Eigengene ~ Trait, id = "DP", family=gaussian, data = model_data, corstr="exchangeable")

          #estimate of the variance for beta
          beta_var = model_kc$cov.beta[["Trait0"]]
          std_error = sqrt(beta_var)
          z = estimate/std_error
          p = 2 * pnorm(-abs(z))
        
          results_to_keep <- data.frame(Comparison = comparison, Module = module, Beta = estimate, Pvalue = p)

          gee_results <- rbind(gee_results, results_to_keep) # nolint
      
          NULL  # return nothing on success
        }, error = function(e) { #if the model fails, we save this to the failed models
          message("Error fitting model for ", module, comparison, ": ", e$message)
          failed_models <<- c(failed_models, paste(module, comparison))
          NULL
        })}

#write these results to gee_results, should ahve all the other comparisons also already included
write.csv(gee_results, "gee_results.csv", row.names = FALSE)


#save all the genes and which modules they belong to, sorted: 
modules_sorted = module_df %>% arrange(colors)

write_delim(modules_sorted,
            file = "sorted_gene_modules.txt",
            delim = "\t")


############# YAY, ONE CELL TYPE DONE ############################################################

############# ############# ############# ############# ############# ############# ############# 

############# ############# ############# ############# ############# ############# ############# 


#if you are analyzing a different cell type, and filtering of genes still leaves a high number (>300 or so ), you can just to keep only the genes with the highest coefficients of variance

cv_values <- apply(as.matrix(counts), 1, function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
})

counts$CV <- cv_values

# Sort and subset top 200 most variable genes
top200 <- counts[order(counts$CV, decreasing = TRUE), ][1:200, ]
