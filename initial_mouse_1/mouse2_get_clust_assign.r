library(chai)
library(scSHC)
library(RaceID)
library(SC3)
library(SingleCellExperiment)
#> Warning: package 'GenomeInfoDb' was built under R version 4.3.3
library(CHOIR)
setwd("/Users/lodimk2/Documents/CMSC691_CHAI_GNN/initial_mouse_2")


# Load data
data("baron_mouse_1")

# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(baron_mouse_1)))
# Add logcounts 
sce <- scuttle::logNormCounts(sce)

# Run All Clustering Algorithms
sce <- get_clust_assignments(sce, n_cores = 12, svm_max = 200, max_k = 15)

create_matrix <- function(assignment_df) {

    clusters <- unique(assignment_df)
    # Create an empty matrix
    similarity_matrix <- matrix(0, nrow = nrow(assignment_df), ncol = nrow(assignment_df))

    # Assign similarity scores
    for (i in 1:length(clusters)) {
    idx <- assignment_df$clust_assign == clusters[i]
    similarity_matrix[idx, idx] <- 1
    }

   diag(similarity_matrix) <- 1

   return(similarity_matrix)
}

create_matrix_list <- function(sce) {
    dataframe <- as.data.frame(colData(sce))
    print(dataframe)
    sim_matrix_list <- list()
    
    for (col_name in names(dataframe)) {
        if (grepl("_assign", col_name)) {
            clusters <- unique(dataframe[[col_name]])
            similarity_matrix <- matrix(0, nrow = nrow(dataframe), ncol = nrow(dataframe))
            for (i in 1:length(clusters)) {
            idx <- dataframe[[col_name]] == clusters[i]
            similarity_matrix[idx, idx] <- 1
            }
            diag(similarity_matrix) <- 1
            #sim_mat <- create_matrix(dataframe[[col_name]])
            write.csv(similarity_matrix, paste0("alg_assignments/",col_name,"_matrix.csv"))
            sim_matrix_list[[col_name]] <- similarity_matrix
        }
    }
    
    print("Similarity Matrix List Created")
    return(sim_matrix_list)
}


# Specify the directory path
directory_path <- "alg_assignments"

# Check if the directory exists
if (!file.exists(directory_path)) {
  # If the directory doesn't exist, create it
  dir.create(directory_path)
  print(paste("Directory", directory_path, "created successfully."))
} else {
  print(paste("Directory", directory_path, "already exists."))
}

similarity_matrix_list = create_matrix_list(sce)

best_k <- 10


avgsim_time <- system.time(sce <- CHAI_AvgSim(sce,best_k,eval = FALSE))

avgsim_time_sce <- CHAI_SNF(sce,best_k,eval = FALSE)

write.csv(data.frame(colData(sce))$CHAI_AvgSim_assign, "chai_avgsim_labels.csv")

write.csv(data.frame(colData(sce))$chai_snf_assign, "chai_snf_labels.csv")
