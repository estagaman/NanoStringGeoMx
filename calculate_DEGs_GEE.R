#calculate_DEGs_GEE.R: calculate differentially expressed genes between diagnoses within one cell type

#INPUT: 
    #countFile: path to a csv file with genes as rows and samples as columns 
    #output_dir: directory you would like output to go to. The script will create an additional directory in this location for the specific cell type analyzed 
    #metadataFile: path to a csv file with samples under column Sample_ID, patient identifier in column DP, cell type in column Custom_CellType, and Diagnosis in column Custom_Diagnosis
    #filter: a string with the cell_type you would like to analyze- must match an entry in column "Custom_CellType" of metadata file
    #EdgeR_filt: boolean equals TRUE if you want EdgeR filtering
    #immune_genes: boolean equals TRUE if you want to filter down to specified immune genes 

#OUTPUT 
    #folder with same name as "filter" provided by input 
    #inside of the folder, you will find: 
        #one csv file per diagnosis contrast, including each gene, its beta value, P-value, adjusted P-value (q), and which diagnosis had higher expression

#Methods Applied: 
    #filtering with EdgeR
    #filtering to immune genes of interest
    #TMM and voom normalization with weight calculation
    #GEE with kauermann-carroll correction 
    #FDR correction by contrast

#load all necessary packages
library("tidyverse")
library("edgeR")
library("limma")
library("preprocessCore")
library("readxl")
library("data.table")
library("glmtoolbox")
library("geesmv")

countFile = "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/unfiltered_analysis/batch_correction/batch_effect/quant_06_16_25/counts_w_neg.csv"
output_dir <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/final/NG2_11_04_25/DEG_GEE_gaussian"
metadataFile <- "/Users/elise/Downloads/multi_omics_erythroderma/spatial.trans/results/metadata/skin_layer/all_samples.csv"

filter = "Keratinocyte" #cell type you would like to analyze for differences

#specify how I want to filter the genes down
edgeR_filt <- TRUE #I want to do EdgeR filtering - deduces which genes have any chance of being found as differentiallly expressed
immune_genes <- TRUE #I want to only test immune-related genes from the list of immune genes

countFile <- read.csv(countFile) %>% as.data.frame()
metadata <- fread(metadataFile) %>% as.data.frame()

#set working directory to the output directory, create it if it doesn't exist
if (!dir.exists(output_dir)) { 
  dir.create(output_dir)}
setwd(output_dir)

###### DATA CLEANING ###########################

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

#remove the samples that do not have custom diagnosis assigned (not important to this study)
metadata <- subset(metadata, !(is.na(Custom_Diagnosis)))
counts <- counts[, metadata$Sample_ID]

#remove negative controls from the counts data
counts <- counts[!grepl("^RTS", rownames(counts)), ]

#save the names of the probes and the samples
probes_save <- rownames(counts)
samples_save <- colnames(counts)

#save the full counts table and full metadata
counts_save <- counts
metadata_save <- metadata

########## DATA CLEANING DONE ##########

########## FILTERING STEP ##########

##ok, I just want to test keratinocyte ROIs, so I am filtering for that here
metadata <- subset(metadata, Custom_CellType == filter)
counts <- counts[ ,metadata$Sample_ID]

probes <- rownames(counts) #saving the probe names and sample names again
samples <- colnames(counts)


##Now, I am doing EdgeR filtering to leave out any genes with loo low of signal to test
design <- model.matrix(~0 + Custom_Diagnosis, data = metadata) #make a design matrix

if (edgeR_filt == TRUE){ #in this case EdgeR_filt is TRUE because I set it at the beginning of the script

    keep <- filterByExpr(counts, design, min.total.count = 15, min.count = 10) #this filters genes down by expression level

} else {
    keep <- c(rep(TRUE, nrow(counts)))
}

colnames(counts) <- samples #make sure rownames of counts are set correctly
rownames(counts) <- probes #same with sample names

#filter down to the genes kept by EdgeR
counts <- counts[keep, ]

##Now, I am going to voom-transform the counts to log counts per million
#Note: voom also computes "precision weights" which account for differences in variance across each sample
#this can help account for heteroscedasticity, as these weights can be plugged back into the model 

dge <- DGEList(counts = counts) #convert to DGE List

dge_tmm <- calcNormFactors(dge, method = "TMM") #perform TMM normalization - preferred by voom-limma

v <- voom(dge_tmm, design, plot = FALSE) #voom normalization, which incorporates TMM factors and calculates weights

weights <- v$weights #extract the weights
new_counts <- v$E #extract the normalized counts

#make sure weights axes are named for later
rownames(weights) <- rownames(counts)
colnames(weights) <- colnames(counts)

#filter down to only immune-related genes
if (immune_genes == TRUE){
    genes_to_annotate <- probes #get a list of all genes in our counts table

    #load in files of immune-related genes from literature
    genes_immune <- read_xls("/Users/elise/Downloads/innatedb_curated_genes.xls", col_names = TRUE) 
    genes_imm <- read_xlsx("/Users/elise/Downloads/InnateDB_genes.xlsx")
    septic_genes <- read_xlsx("/Users/elise/Downloads/septic_shock_genes.xlsx")

    #extract the names of the genes, format them correctly
    genes_imm <- genes_imm$name
    septic_genes <- septic_genes$name

    genes_immune <- genes_immune[, "Gene Symbol"]
    genes_immune <- genes_immune %>% as.data.frame() %>% unique()
    symbols_immune <- genes_immune[, "Gene Symbol"]

    #find which genes from these lists are in my counts table
    all_genes_immune <- c(symbols_immune, genes_imm, septic_genes)
    all_genes_immune <- toupper(all_genes_immune)
    all_genes_immune <- unique(all_genes_immune)

    #find the intersecting genes
    my_genes_immune <- intersect(genes_to_annotate, all_genes_immune)
    length(my_genes_immune) #down to 2451 genes that are in the innate immunity, immunity, and septic shock databases 

    #now load in the genes from the paper
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
    counts <- counts[rownames(counts) %in% immune_plus_paper, ]
    new_counts <- new_counts[rownames(new_counts) %in% immune_plus_paper, ]
    weights <- weights[rownames(weights) %in% immune_plus_paper, ]
}

#create a transposed version of the counts table, which I don't think we'll need but just in case
counts_t <- t(counts)

######## FILTERING DONE ########################

######## STATISTICAL TESTING ########################
#function make the metadata table with all possible contrasts I want to test
make_comparisons_df <- function(metadata){
    CTCLvsOther <- ifelse(metadata$Custom_Diagnosis == "CTCL", 1, 0) #CTCLvsOther

    PSOvsOther <- ifelse(metadata$Custom_Diagnosis == "PSO", 1, 0) #PSOvsOther

    ADvsOther <- ifelse(metadata$Custom_Diagnosis == "AD", 1, 0) #ADvsOther

    DRvsOther <- ifelse(metadata$Custom_Diagnosis == "DR", 1, 0) #DRvsOther

    all_groups <- cbind(CTCLvsOther, PSOvsOther, ADvsOther, DRvsOther) %>% as.data.frame() #bind all of these columns together into a dataframe

    rownames(all_groups) <- metadata$Sample_ID #name each row by the sample

    return(all_groups) #return this dataframe
}

#function to perform the GEE differential expression testing
perform_gee <- function(counts, new_counts, metadata, all_groups, filter, weights){

  # Create somewhere to save the GEE results
  hub_gene_results <- data.frame(
    Comparison = character(),
    Gene = character(),
    Beta = numeric(),
    Pvalue = numeric(),
    Higher = character(),
    stringsAsFactors = FALSE
  )

  # Create somewhere to save the failed results caught by the TryCatch()
  failed_results <- c()

  # For each contrast: ie CTCLvsOther, PSOvsOther, DRvsOther, etc.
  for (comparison in colnames(all_groups)) {

    cond1 = substr(comparison, start = 1, stop = nchar(comparison) - 7) #cond1 will isolate one side of the contrast ie. for CTCLvsOther, cond1 = "CTCL"

    # We want to test each gene individually
    for (gene in rownames(counts)) {

      # Create data for the model
      model_data <- data.frame(
        gene_test = as.numeric(new_counts[gene, ]), #take the normalized counts for that gene
        weights = as.numeric(weights[gene, ]), #take the weights for that gene
        Trait = relevel(as.factor(all_groups[[comparison]]), ref = "1"), #include the Diagnosis associated with each sample
        DP = as.factor(metadata$DP) #include the Patient/DP as a random effect
      )

      # Try to fit the model and compute statistics
      tryCatch({

        # Fit the GEE model
        gee_model <- glmgee(
          formula = gene_test ~ Trait, #test with gene expression as response and Diagnosis as the predictor
          id = DP, #patient ID as the clustering level
          data = model_data, #use the model_data we generated
          corstr = "Exchangeable", #assuming all ROIs from same patient are equally correlated (not necessarily true for spatial data)
          family = gaussian(), #use gaussian distribution because the counts are non-integers (normalized) and precision weights incorporated for heteroskedasticity
          weights = weights #incorporate weights
        )

        # Extract the estimate for Trait0 "Other"
        estimate <- coef(gee_model)["Trait0",] #extract the estimate from this model

        # Adjust the beta and variance with KC estimator for low sample size and low cluster size, parameters same as before
        model_kc <- GEE.var.kc(
          formula = gene_test ~ Trait,
          id = "DP",
          family = gaussian,
          data = model_data,
          corstr = "exchangeable"
        )

        beta_var <- model_kc$cov.beta[["Trait0"]] #extract the new variance of beta 
        std_error <- sqrt(beta_var) #get the standard error
        z <- estimate / std_error #calculate z-score
        p <- 2 * pnorm(-abs(z)) #calculate two-sided p-value

        #calculate which group in our test is higher
        higher <- ifelse(estimate > 0, "Other", cond1)

        # Save successful results
        results_to_keep <- data.frame(
          Comparison = comparison,
          Gene = gene,
          Beta = estimate,
          Pvalue = p, 
          Higher = higher
        )

        #bind this to our results for all genes
        hub_gene_results <- rbind(hub_gene_results, results_to_keep)

      }, error = function(e) {
        # If there was an error, log it and save "failed" results
        message("Error fitting model for ", gene, " ", comparison, ": ", e$message)
        failed_results <<- c(failed_results, paste(gene, comparison))

        results_to_keep <- data.frame(
          Comparison = comparison,
          Gene = gene,
          Beta = NA,
          Pvalue = NA,
          Higher = NA
        )
        hub_gene_results <- rbind(hub_gene_results, results_to_keep)

      })
    }
  }

  # When done, return all the results at once
  return(hub_gene_results)
}

#this function is for if you want to compare between 2 specified groups, instead of one group vs all other
perform_gee_2groups <- function(counts, new_counts, metadata, all_groups, filter, weights, group1, group2){

  #filter down to just those two groups
  metadata_filt <- subset(metadata, metadata$Custom_Diagnosis == group1 | metadata$Custom_Diagnosis == group2)

  #filter down counts and weights to match
  new_counts_filt <- new_counts[ ,metadata_filt$Sample_ID]
  weights_filt <- weights[ ,metadata_filt$Sample_ID]

  #data is already voom normalized, weights already calculated, already filtered down to immune genes 

  # Create somewhere to save the GEE results
  hub_gene_results <- data.frame(
    Comparison = character(),
    Gene = character(),
    Beta = numeric(),
    Pvalue = numeric(),
    Higher = character(),
    stringsAsFactors = FALSE
  )

  # Create somewhere to save the failed results caught by the TryCatch()
  failed_results <- c()

  #create new vector for Diagnosis. Only testing one comparison 
  diag_vector <- ifelse(metadata_filt$Custom_Diagnosis == group1, 1, 0) %>% as.factor()

  # We want to test each gene individually
  for (gene in rownames(counts)) {

    # Create data for the model
    model_data <- data.frame(
      gene_test = as.numeric(new_counts_filt[gene, ]), #take the normalized counts for that gene
      weights = as.numeric(weights_filt[gene, ]), #take the weights for that gene
      Trait = relevel(diag_vector, ref = "1"), #include the Diagnosis associated with each sample
      DP = as.factor(metadata_filt$DP) #include the Patient/DP as a random effect
    )

      # Try to fit the model and compute statistics
    tryCatch({

      # Briefly use glm.nb to find the dispersion parameter theta
      #fitting_for_theta <- glm.nb(formula = gene_test ~ Trait, data = model_data)
      #theta_fitted <- fitting_for_theta$theta
      #^ I would only do this if I were doing negative binomial. In this case, my counts are transformed and use precision weights, so I am using gaussian

      # Fit the GEE model
      gee_model <- glmgee(
        formula = gene_test ~ Trait, #test with gene expression as response and Diagnosis as the predictor
        id = DP, #patient ID as the clustering level
        data = model_data, #use the model_data we generated
        corstr = "Exchangeable", #all ROIs from same patient are equally correlated
        family = gaussian(), #use gaussian distribution because the counts are non-integers (normalized) and precision weights incorporated for heteroskedasticity
        weights = weights #incorporate weights
      )

        # Extract the estimate for Trait0
      estimate <- coef(gee_model)["Trait0",] #extract the estimate from this model

      # Adjust the beta and variance with KC estimator, parameters same as before
      model_kc <- GEE.var.kc(
        formula = gene_test ~ Trait,
        id = "DP",
        family = gaussian,
        data = model_data,
        corstr = "exchangeable"
      )

      beta_var <- model_kc$cov.beta[["Trait0"]] #extract the new variance of beta 
      std_error <- sqrt(beta_var) #get the standard error
      z <- estimate / std_error #calculate z-score
      p <- 2 * pnorm(-abs(z)) #calculate two-sided p-value

      #calculate which group had higher expression of the gene
      higher = ifelse(estimate > 0, group2, group1)

      # Save successful results
      results_to_keep <- data.frame(
        Comparison = paste0(group1, "vs", group2),
        Gene = gene,
        Beta = estimate,
        Pvalue = p, 
        Higher = higher
      )
      hub_gene_results <- rbind(hub_gene_results, results_to_keep)

    }, error = function(e) {
      # If there was an error, log it and save "failed"
      message("Error fitting model for ", gene, ": ", e$message)
      failed_results <<- c(failed_results, gene)

      results_to_keep <- data.frame(
        Comparison = paste(group1, "vs", group2),
        Gene = gene,
        Beta = NA,
        Pvalue = NA,
        Higher = NA
      )

      hub_gene_results <- rbind(hub_gene_results, results_to_keep)
    })
  }

  # When done, return all the results at once
  return(hub_gene_results)

}

#perform FDR correction on a dataframe of genes and p-values, containing multiple contrasts

FDR_correction <- function(all_results){

  #split the results by contrast
  results_split <- split(all_results, all_results$Comparison)

  #create a list for the new adjusted results
  adjusted_results = list()

  #for each contrast
  for (i in 1:length(results_split)){

    #arrange results in order of p_value 
    table = results_split[i] %>% unname() %>% as.data.frame() %>% arrange(Pvalue, desc = FALSE)

    #adjust the pvalues in new column FDR_adj_P
    table$FDR_adj_P = p.adjust(table$Pvalue, method = "BH")

    #create a column noting whether expression is differentially significant
    table$DE = ifelse(table$FDR_adj_P < 0.05, "DE", "not DE")

    #save this new table to the same index in our list of tables
    adjusted_results <- c(adjusted_results, list(table))
  }

  #return the adjusted results
  return(adjusted_results)

}

all_groups_df <- make_comparisons_df(metadata)
hub_gene_results <- perform_gee(counts, new_counts, metadata, all_groups_df, "Keratinocyte", weights)
hub_gene_results_ADvsPSO <- perform_gee_2groups(counts, new_counts, metadata, all_groups, filter, weights, group1 = "AD", group2 = "PSO")

#combine all the results together as needed
all_results <- rbind(hub_gene_results, hub_gene_results_ADvsPSO)

#do FDR correction
all_results_FDR_corrected <- FDR_correction(all_results)

#save all of these results to the DEG folder
if (!dir.exists(filter)) { 
  dir.create(filter)}

for (i in 1:length(all_results_FDR_corrected)){
  compare_name = all_results_FDR_corrected[[i]]$Comparison[1]
  write.csv(all_results_FDR_corrected[[i]], paste0(filter, "/DE_genes_", compare_name, ".csv"))
}

############# DONE ##################################################
