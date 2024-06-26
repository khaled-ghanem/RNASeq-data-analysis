#1) Load the required libraries----
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(apeglm)
library(NOISeq)
library(tximport)
#2) Read meta data----
# Read the ReadSetTable
read_set_table <- read_excel("meta/ReadSetTable_Sample2_01182024.xlsx")

# Read the TC_Table
tc_table <- read_excel("meta/TC_table_sample2_01182024.xlsx")

# Select the 'ReadSetID' and ' SampleID' columns from the 'read_set_table' dataframe
study_design <- dplyr::select(read_set_table, "ReadSetID", "SampleID")

study_design <- as.data.frame(study_design)
rownames(study_design) <- study_design[,1]

#3) extract the expression data----
# Set the source and destination folder paths
source_folder <- "data2/"
destination_folder <- "extracted_data/"

# Create the destination folder if it doesn't exist
if (!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

# Iterate through each file and copy it to the destination folder
for (set_id in study_design$ReadSetID) {
  # Construct the source and destination file paths
  source_file <- paste0(source_folder, set_id, "_ST/t_data.ctab.gz")
  destination_file <- paste0(destination_folder, set_id, "_t_data.ctab")
  
  # Copy the file to the destination folder
  #file.copy(source_file, destination_file)
  
  # If the file is compressed, extract it using gzfile
  if (file.exists(source_file)) {
    con <- gzfile(source_file, "rb")
    content <- readLines(con)
    close(con)
    
    # Save the content to the same file with .tsv extension
    writeLines(content, destination_file)
    
  } else {
    cat("Error: File", source_file, "not found.\n")
  }
}
# Create a vector of file paths
files <- paste0(destination_folder, study_design$ReadSetID, "_t_data.ctab")


#4) Apply DE using DESeq2----

# Function to check if a file name contains any word from the list
contains_word <- function(file_name, words) {
  any(sapply(words, function(word) grepl(word, file_name, ignore.case = TRUE)))
}

# Create the result folder if it doesn't exist
if (!dir.exists("results/")) {
  dir.create("results/")
}

# Create the normalized_data folder if it doesn't exist
if (!dir.exists("normalized_data/")) {
  dir.create("normalized_data/")
}

i <- 1
for (TC_pair in tc_table$`TC sets`) {
  cat("iteration number", i ,"\n")
  i <- i+1
  TC_pair_vector <- strsplit(TC_pair, "_")[[1]]
  TC_pair_study_design <- study_design[study_design$SampleID %in% TC_pair_vector, ]
  TC_pair_study_design$label <- as.factor(ifelse(TC_pair_study_design$SampleID == TC_pair_vector[1], 'T' , 'C'))
  
  # Create a pattern by combining all ReadSetID elements with "|"
  pattern <- paste(TC_pair_study_design$ReadSetID, collapse = "|")
  # Use grep to select files based on the condition
  TC_pair_files <- files[grep(pattern, files)]
  
  #print(TC_pair_files)
  #print(TC_pair_study_design$ReadSetID)
  
  tmp <- read_tsv(files[1])
  tx2gene <- tmp[, c("t_name", "gene_id")]
  # if you want to run DTE just uncomment line 93
  txi <- tximport(TC_pair_files,
                  type = "stringtie", 
                  tx2gene = tx2gene , 
                  #txOut=FALSE
                  )
  colnames(txi$abundance) <- TC_pair_study_design$ReadSetID
  colnames(txi$length) <- TC_pair_study_design$ReadSetID
  colnames(txi$counts) <- TC_pair_study_design$ReadSetID
  # Filter out rows where 80% or more of the columns have a value equal to 0 in that row
  selected_rows <- rowMeans(txi$counts == 0) < 0.8
  txi$counts <- txi$counts[selected_rows, ]
  txi$abundance <- txi$abundance[selected_rows, ]
  txi$length <- txi$length[selected_rows, ]
  
  if(any(duplicated(TC_pair_study_design$SampleID)))
  {
    dds <- DESeqDataSetFromTximport(txi,
                                  colData = TC_pair_study_design,
                                  design = ~ label)
    
  
    
    # Perform size factor estimation and normalization
    dds <- DESeq(dds)
    # Extract normalized counts
    normalized_counts <- counts(dds, normalized = TRUE)
    
    # Perform differential expression analysis
    results <- results(dds)
    #summary(results)
    
    results <- lfcShrink(dds, coef = "label_T_vs_C" , type = "apeglm")
    
    # Extract significant differentially expressed genes
    DE_genes <- subset(results, padj < 0.05 & abs(log2FoldChange) > 1)
    
  }
  else
  {
    # create NOISeq object
    mydata <- readData(data = txi$counts, 
                       length = txi$length[,1],  
                       factors = TC_pair_study_design
                       )
    # apply tmm normalization
    normalized_counts = tmm(assayData(mydata)$exprs , 
                long = txi$length[,1] , 
                lc = 0.5 
                )
    
    # Detect differentially expressed features
    mynoiseq = noiseq(mydata,norm = "tmm", factor = "label")
    #head(mynoiseq@results[[1]])
    
    # Select differentially expressed
    DE_genes = degenes(mynoiseq, q = 0.8, M = NULL)   
    
    ## Explanation of NOISeq Output
    # M: Log2 fold change between treatment and control groups (T_mean - C_mean).
    # D: The NOISeq dispersion measure.
    # prob: The probability of differential expression.
    # ranking: The ranking of the gene based on the probability of differential expression.
    
  }
  
  
  DE_genes.dataframe <- as.data.frame(DE_genes)
  #View(DE_genes.dataframe)
  write.csv(DE_genes.dataframe, paste('results/',TC_pair,'_Pair_DE_Genes.xlsx'))
  
  # save the preprocessed counts 
  write.csv(normalized_counts, paste('normalized_data/',TC_pair,'_Pair_Normalized_Expression.xlsx'))
  
}






