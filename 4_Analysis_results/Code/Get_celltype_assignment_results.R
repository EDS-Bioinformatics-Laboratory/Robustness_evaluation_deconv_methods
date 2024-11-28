################################################################################
# Where do the deconvolution method assigns the missing cell types proportion?
# 
# Celltype assignment plot
# We calculate the difference between the predictions of a deconvolution method
# with a removal of one or more celltype/s from reference dataset and
# the predictions at baseline scenario in mentioned order
# 
# The increase in a particular celltype proportion infers that the celltype was
# assigned to the cells when the actual one is missing from reference data
# 
# The decrease in a particular celltype proportion infers that the deconvolution
# method failed to identify the cells belonging to a cell type present in
# reference data compared to how it performed at baseline
# 
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################

# Initialize environment
source("Init_env.R")

args <- commandArgs(trailingOnly = TRUE)
st_num <- args[1]
rm_scenarios <- strsplit(args[2], ",")[[1]]

# all deconvolution methods are considered if command line arg is missing
if (length(args)<=2) {
  methods_ <- c("cell2location",
                "RCTD",
                "CARD",
                "SCDC",
                "MuSiC",
                "Stereoscope",
                "Seurat",
                "SPOTlight")
} else {
  methods_ <- strsplit(args[3], ",")[[1]]
}


# path to deconvolution methods results for baseline scenarios
b_Method_Res <- "../../3_ST_methods/Results/rm0/"
decon_results <- "../../3_ST_methods/Results/"


# function for calculating cell type assignment
get_celltype_reassignment <- function(base_, pred_) {
  
  missing.ct <- setdiff(colnames(base_), colnames(pred_))
  common.ct <- intersect(colnames(base_), colnames(pred_))
  
  for (q in 1:dim(base_)[1]) {
    for (s in colnames(base_)) {
      if (s %in% common.ct) {
        pred_[q, s] <- pred_[q, s] - base_[q, s]
      }
    }
  }
  
  d <- data.frame(colSums(pred_)/sum(colSums(base_)[missing.ct]))
  colnames(d) <- c("value")
  d <- tibble::rownames_to_column(d, "celltype")
  rownames(d) <- d$celltype
  
  # assign the proportion total from baseline results to the calculated colSums
  d[missing.ct, 1] <- missing.ct
  
  # ordering row names
  d <- d[order(row.names(d)), ]
  
  # returns a row in matrix stored in dflist list
  return(d$value)
}


# reading the deconvolution method results
get_celltype_results <- function(z, y, Method_Res, dflist_, num_ct) {
  # z: index of ST dataset
  # y: index of sc reference dataset for the removal scenario
  # Method_Res: path to deconvolution methods results for removal scenarios
  # num_ct: number of celltypes in sc-reference dataset for removal scenario
  # y value is fixed for baseline results as theres only one reference dataset
  
  counter_ <- 1 # to store results in the array
  
  # dflist[ , , ] --> removed celltype, all celltypes, method
  
  if ("cell2location" %in% methods_) {
    results.Cell2Loc <- 
      utils::read.csv(paste0(Method_Res, "Cell2Location/Cell2Location.", z, "-", y, ".csv"), row.names = 1)
    colnames(results.Cell2Loc) <- sub("q05cell_abundance_w_sf_", "", colnames(results.Cell2Loc))
    row_sum <- rowSums(results.Cell2Loc)
    results.Cell2Loc <- results.Cell2Loc %>% mutate_all(~ . /row_sum)
    
    b.results.Cell2Loc <- 
      utils::read.csv(paste0(b_Method_Res, "Cell2Location/Cell2Location.", z, "-", "1.csv"), row.names = 1)
    colnames(b.results.Cell2Loc) <- sub("q05cell_abundance_w_sf_", "", colnames(b.results.Cell2Loc))
    row_sum <- rowSums(b.results.Cell2Loc)
    b.results.Cell2Loc <- b.results.Cell2Loc %>% mutate_all(~ . /row_sum)
    
    dflist_[y, , counter_] <- get_celltype_reassignment(b.results.Cell2Loc, results.Cell2Loc)
    counter_ <- counter_ + 1
  }
  
  if ("RCTD" %in% methods_) {
    results.RCTD <- read.csv(paste0(Method_Res, "RCTD/results.RCTD.", z, "-", y, ".csv"),
                             row.names = 1)
    b.results.RCTD <- read.csv(paste0(b_Method_Res, "RCTD/results.RCTD.", z, "-", "1.csv"),
                               row.names = 1)
    
    dflist_[y, , counter_] <- get_celltype_reassignment(b.results.RCTD, results.RCTD)
    counter_ <- counter_ + 1
  }
  
  if ("CARD" %in% methods_) {
    if (file.exists(paste0(Method_Res, "CARD/results.CARD.", z, "-", y, ".csv")) != T) {
      results.CARD <- matrix(NA, nrow = 1600, ncol = num_ct) %>% data.frame()
      colnames(results.CARD) <- colnames(results.RCTD)
      rownames(results.CARD) <- sprintf("Spot_%s",seq(1:1600))
      message("CARD ", z, y)
    } else {
      results.CARD <- read.csv(paste0(Method_Res, "CARD/results.CARD.", z, "-", y, ".csv"), row.names = 1)
    }
    
    if (file.exists(paste0(b_Method_Res, "CARD/results.CARD.", z, "-", "1.csv")) != T) {
      b.results.CARD <- matrix(NA, nrow = 1600, ncol = 13) %>% data.frame()
      colnames(b.results.CARD) <- colnames(b.results.RCTD)
      rownames(b.results.CARD) <- sprintf("Spot_%s",seq(1:1600))
      message("CARD ", z, y)
    } else {
      b.results.CARD <- read.csv(paste0(b_Method_Res, "CARD/results.CARD.", z, "-", "1.csv"), row.names = 1)
    }
    
    dflist_[y, , counter_] <- get_celltype_reassignment(b.results.CARD, results.CARD)
    counter_ <- counter_ + 1
  }
  
  if ("SCDC" %in% methods_) {
    results.SCDC <- read.csv(paste0(Method_Res, "SCDC/results.SCDC.", z, "-", y, ".csv"), row.names = 1)
    b.results.SCDC <- read.csv(paste0(b_Method_Res, "SCDC/results.SCDC.", z, "-", "1.csv"), row.names = 1)
    
    dflist_[y, , counter_] <- get_celltype_reassignment(b.results.SCDC, results.SCDC)
    counter_ <- counter_ + 1
  }
  
  if ("MuSiC" %in% methods_) {
    results.MuSiC <- read.csv(paste0(Method_Res, "MuSiC/results.MuSiC.", z, "-", y, ".csv"), row.names = 1)
    b.results.MuSiC <- read.csv(paste0(b_Method_Res, "MuSiC/results.MuSiC.", z, "-", "1.csv"), row.names = 1)
    
    dflist_[y, , counter_] <- get_celltype_reassignment(b.results.MuSiC, results.MuSiC)
    counter_ <- counter_ + 1
  }
  
  if ("Stereoscope" %in% methods_) {
    results.Stereoscope <- 
      utils::read.csv(paste0(Method_Res, "Stereoscope/Stereoscope_result.", z, "-", y, ".csv"), row.names = 1)
    b.results.Stereoscope <- 
      utils::read.csv(paste0(b_Method_Res, "Stereoscope/Stereoscope_result.", z, "-", "1.csv"), row.names = 1)
    
    dflist_[y, , counter_] <- get_celltype_reassignment(b.results.Stereoscope, results.Stereoscope)
    counter_ <- counter_ + 1
  }
  
  if ("Seurat" %in% methods_) {
    if (file.exists(paste0(Method_Res, "Seurat/results.Seurat.", z, "-", y, ".csv")) != T) {
      results.Seurat <- matrix(NA, nrow = 1600, ncol = num_ct) %>% data.frame()
      colnames(results.Seurat) <- colnames(b.results.RCTD)
      rownames(results.Seurat) <- sprintf("Spot_%s",seq(1:1600))
      message("Seurat ", z, y)
    } else {
      results.Seurat <- read.csv(paste0(Method_Res, "Seurat/results.Seurat.", z, "-", y, ".csv"), row.names = 1)
      colnames(results.Seurat) <- gsub("\\.", "_", colnames(results.Seurat))
    }
    
    if (file.exists(paste0(b_Method_Res, "Seurat/results.Seurat.", z, "-", "1.csv")) != T) {
      b.results.Seurat <- matrix(NA, nrow = 1600, ncol = 13) %>% data.frame()
      colnames(b.results.Seurat) <- colnames(b.results.RCTD)
      rownames(b.results.Seurat) <- sprintf("Spot_%s",seq(1:1600))
      message("Seurat ", z, y)
    } else {
      b.results.Seurat <- read.csv(paste0(b_Method_Res, "Seurat/results.Seurat.", z, "-", "1.csv"), row.names = 1)
      colnames(b.results.Seurat) <- gsub("\\.", "_", colnames(b.results.Seurat))
    }
    
    dflist_[y, , counter_] <- get_celltype_reassignment(b.results.Seurat, results.Seurat)
    counter_ <- counter_ + 1
  }
  
  if ("SPOTlight" %in% methods_) {
    results.SPOTlight <- read.csv(paste0(Method_Res, "SPOTlight/results.SPOTlight.", z, "-", y, ".csv"),
                                  row.names = 1)
    b.results.SPOTlight <- read.csv(paste0(b_Method_Res, "SPOTlight/results.SPOTlight.", z, "-", "1.csv"),
                                    row.names = 1)
    
    dflist_[y, , counter_] <- get_celltype_reassignment(b.results.SPOTlight, results.SPOTlight)
    counter_ <- counter_ + 1
  }
  
  return(dflist_)
}


if ("rm1" %in% rm_scenarios) {
  # removal of one celltype
  Method_Res <- paste0(decon_results, "rm1/")
  
  # param: total # sc reference datasets, total # celltypes, # methods
  dflist_ <- array(-99, dim = c(13, 13, length(methods_)))
  
  for (y in 1:13) {
    dflist_ <- get_celltype_results(st_num, y, Method_Res, dflist_, 12)
  }
  saveRDS(dflist_, file = paste0(Results, "CT_assign_rm1.", st_num,"_ST.rds"))
}


if ("rm2" %in% rm_scenarios) {
  # removal of two celltypes
  Method_Res <- paste0(decon_results, "rm2/")
  
  # param: total # sc reference datasets, total # celltypes, # methods
  dflist_ <- array(-99, dim = c(5, 13, length(methods_)))
  
  for (y in 1:5) {
    dflist_ <- get_celltype_results(st_num, y, Method_Res, dflist_, 11)
  }
  saveRDS(dflist_, file = paste0(Results, "CT_assign_rm2.", st_num,"_ST.rds"))
}


if ("rm3" %in% rm_scenarios) {
  # removal of three celltypes
  Method_Res <- paste0(decon_results, "rm3/")
  
  # param: total # sc reference datasets, total # celltypes, # methods
  dflist_ <- array(-99, dim = c(5, 13, length(methods_)))
  
  for (y in 1:5) {
    dflist_ <- get_celltype_results(st_num, y, Method_Res, dflist_, 10)
  }
  saveRDS(dflist_, file = paste0(Results, "CT_assign_rm3.", st_num,"_ST.rds"))
}


# if ("rm5" %in% rm_scenarios) {
#   # removal of five celltypes
#   Method_Res <- paste0(decon_results, "rm5/")
#   
#   # param: total # sc reference datasets, total # celltypes, # methods
#   dflist_ <- array(-99, dim = c(1, 13, length(methods_)))
#   
#   for (y in 1:1) {
#     dflist_ <- get_celltype_results(st_num, y, Method_Res, dflist_, 8)
#   }
#   saveRDS(dflist_, file = paste0(Results, "CT_assign_rm5.RDS"))
# }
# 
# if ("rm10" %in% rm_scenarios) {
#   # removal of ten celltypes
#   Method_Res <- paste0(decon_results, "rm10/")
#   
#   # param: total # sc reference datasets, total # celltypes, # methods
#   dflist_ <- array(-99, dim = c(5, 13, length(methods_)))
#   
#   for (y in 1:5) {
#     dflist_ <- get_celltype_results(st_num, y, Method_Res, dflist_, 3)
#   }
#   saveRDS(dflist_, file = paste0(Results, "CT_assign_rm10.RDS"))
# }
# 
# if ("rm11" %in% rm_scenarios) {
#   # removal of ten celltypes
#   Method_Res <- paste0(decon_results, "rm11/")
#   
#   # param: total # sc reference datasets, total # celltypes, # methods
#   dflist_ <- array(-99, dim = c(5, 13, length(methods_)))
#   
#   for (y in 1:5) {
#     dflist_ <- get_celltype_results(st_num, y, Method_Res, dflist_, 2)
#   }
#   saveRDS(dflist_, file = paste0(Results, "CT_assign_rm11.RDS"))
# }