################################################################################
# This script calculates the difference in the RMSE values for all the removal
# scenarios (wtih cell type mismatch) compared to the baseline scenario (with
# no cell type mismatch)
# 
# We used a vector (length = number of spots in ST data) of calculated RMSE values
# Incase of multiple RMSE vectors (depending on number of single cell reference
# datasets available for a particular removal scenario) for a single method and
# scenario, we considered the mean vector for calculation. 
# 
# 
# RMSE.mean.1_ST.rds, RMSE.mean.2_ST.rds, RMSE.mean.3_ST.rds
# Each one of the above object refers to one ST dataset
# Each object has 7 lists, one for each scenario (baseline + 6 removal scenario)
# Each sub-list is a matrix with columns as methods and rows as spots
# [i,j] is the mean RMSE value of spot i for method j over multiple single cell
# reference datasets
# 
# 
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################

# Initialize environment
source("Init_env.R")

args <- commandArgs(trailingOnly = TRUE)
st_num <- args[1]

# all removal scenarios are considered if command line argument is missing
if (length(args)<=1) {
  rm_scenarios <- c("rm1", "rm2", "rm3", "rm5", "rm10", "rm11")
} else {
  rm_scenarios <- strsplit(args[2], ",")[[1]]
}

# all deconvolution methods are considered if command line argument is missing
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

st_data_path <- "../../2_Simulating_ST_data/Results/Spatial_Data/"
decon_results <- "../../3_ST_methods/Results/"

# reading the deconvolution method results
read_results <- function(z, y, Method_Res) {
  # z refers to ST dataset number
  # y refers to sc reference dataset number for the removal scenario
  # Method_Res refers to path of deconvolution methods results
  
  list_ <- c()
  
  spots.metadata <- readRDS(paste0(st_data_path, "metadata.", z, ".rds"))
  
  if ("cell2location" %in% methods_) {
    results.Cell2Loc <- 
      utils::read.csv(paste0(Method_Res, "Cell2Location/Cell2Location.", z, "-", y, ".csv"), row.names = 1)
    colnames(results.Cell2Loc) <- sub("q05cell_abundance_w_sf_", "", colnames(results.Cell2Loc))
    row_sum <- rowSums(results.Cell2Loc)
    results.Cell2Loc <- results.Cell2Loc %>% mutate_all(~ . /row_sum)
    
    # common celltypes from the ground truth results and scenario results
    missing.ct <- setdiff(colnames(spots.metadata), colnames(results.Cell2Loc))
    common.ct <- intersect(colnames(spots.metadata), colnames(results.Cell2Loc))
    
    results.Cell2Loc <- results.Cell2Loc[, common.ct]
    if (length(missing.ct)!=0) {
      for (c in 1:length(missing.ct)){
        results.Cell2Loc <- cbind(results.Cell2Loc, "_" = rep(0, dim(spots.metadata)[1]))
      }
    }
    colnames(results.Cell2Loc) <- c(common.ct, missing.ct)
    results.Cell2Loc <- results.Cell2Loc[, order(colnames(results.Cell2Loc))]
    
    list_ <- c(list_, list(results.Cell2Loc))
  }
  
  if ("RCTD" %in% methods_) {
    results.RCTD <- read.csv(paste0(Method_Res, "RCTD/results.RCTD.", z, "-", y, ".csv"),
                             row.names = 1)
    
    # common celltypes from the ground truth results and scenario results
    missing.ct <- setdiff(colnames(spots.metadata), colnames(results.RCTD))
    common.ct <- intersect(colnames(spots.metadata), colnames(results.RCTD))
    
    results.RCTD <- results.RCTD[, common.ct]
    if (length(missing.ct)!=0) {
      for (c in 1:length(missing.ct)){
        results.RCTD <- cbind(results.RCTD, "_" = rep(0, dim(spots.metadata)[1]))
      }
    }
    colnames(results.RCTD) <- c(common.ct, missing.ct)
    results.RCTD <- results.RCTD[, order(colnames(results.RCTD))]
    
    list_ <- c(list_, list(results.RCTD))
  }
  
  if ("CARD" %in% methods_) {
    if (file.exists(paste0(Method_Res, "CARD/results.CARD.", z, "-", y, ".csv")) != T) {
      results.CARD <- matrix(NA, nrow = 1600, ncol = 13) %>% data.frame()
      colnames(results.CARD) <- colnames(results.RCTD)
      rownames(results.CARD) <- sprintf("Spot_%s",seq(1:1600))
      message("CARD ", z, y)
    } else {
      results.CARD <- read.csv(paste0(Method_Res, "CARD/results.CARD.", z, "-", y, ".csv"), row.names = 1)
    }
    
    # common celltypes from the ground truth results and scenario results
    missing.ct <- setdiff(colnames(spots.metadata), colnames(results.CARD))
    common.ct <- intersect(colnames(spots.metadata), colnames(results.CARD))
    
    results.CARD <- results.CARD[, common.ct]
    if (length(missing.ct)!=0) {
      for (c in 1:length(missing.ct)){
        results.CARD <- cbind(results.CARD, "_" = rep(0, dim(spots.metadata)[1]))
      }
    }
    colnames(results.CARD) <- c(common.ct, missing.ct)
    results.CARD <- results.CARD[, order(colnames(results.CARD))]
    
    list_ <- c(list_, list(results.CARD))
  }
  
  if ("SCDC" %in% methods_) {
    results.SCDC <- read.csv(paste0(Method_Res, "SCDC/results.SCDC.", z, "-", y, ".csv"), row.names = 1)
    
    # common celltypes from the ground truth results and scenario results
    missing.ct <- setdiff(colnames(spots.metadata), colnames(results.SCDC))
    common.ct <- intersect(colnames(spots.metadata), colnames(results.SCDC))
    
    results.SCDC <- results.SCDC[, common.ct]
    if (length(missing.ct)!=0) {
      for (c in 1:length(missing.ct)){
        results.SCDC <- cbind(results.SCDC, "_" = rep(0, dim(spots.metadata)[1]))
      }
    }
    colnames(results.SCDC) <- c(common.ct, missing.ct)
    results.SCDC <- results.SCDC[, order(colnames(results.SCDC))]
    
    list_ <- c(list_, list(results.SCDC))
  }
  
  if ("MuSiC" %in% methods_) {
    results.MuSiC <- read.csv(paste0(Method_Res, "MuSiC/results.MuSiC.", z, "-", y, ".csv"), row.names = 1)
    
    # common celltypes from the ground truth results and scenario results
    missing.ct <- setdiff(colnames(spots.metadata), colnames(results.MuSiC))
    common.ct <- intersect(colnames(spots.metadata), colnames(results.MuSiC))
    
    results.MuSiC <- results.MuSiC[, common.ct]
    if (length(missing.ct)!=0) {
      for (c in 1:length(missing.ct)){
        results.MuSiC <- cbind(results.MuSiC, "_" = rep(0, dim(spots.metadata)[1]))
      }
    }
    colnames(results.MuSiC) <- c(common.ct, missing.ct)
    results.MuSiC <- results.MuSiC[, order(colnames(results.MuSiC))]
    
    list_ <- c(list_, list(results.MuSiC))
  }
  
  if ("Stereoscope" %in% methods_) {
    results.Stereoscope <- 
      utils::read.csv(paste0(Method_Res, "Stereoscope/Stereoscope_result.", z, "-", y, ".csv"), row.names = 1)
    
    # common celltypes from the ground truth results and scenario results
    missing.ct <- setdiff(colnames(spots.metadata), colnames(results.Stereoscope))
    common.ct <- intersect(colnames(spots.metadata), colnames(results.Stereoscope))
    
    results.Stereoscope <- results.Stereoscope[, common.ct]
    if (length(missing.ct)!=0) {
      for (c in 1:length(missing.ct)){
        results.Stereoscope <- cbind(results.Stereoscope, "_" = rep(0, dim(spots.metadata)[1]))
      }
    }
    colnames(results.Stereoscope) <- c(common.ct, missing.ct)
    results.Stereoscope <- results.Stereoscope[, order(colnames(results.Stereoscope))]
    
    list_ <- c(list_, list(results.Stereoscope))
  }
  
  if ("Seurat" %in% methods_) {
    if (file.exists(paste0(Method_Res, "Seurat/results.Seurat.", z, "-", y, ".csv")) != T) {
      results.Seurat <- matrix(NA, nrow = 1600, ncol = 13) %>% data.frame()
      colnames(results.Seurat) <- colnames(results.RCTD)
      rownames(results.Seurat) <- sprintf("Spot_%s",seq(1:1600))
      message("Seurat ", z, y)
    } else {
      results.Seurat <- read.csv(paste0(Method_Res, "Seurat/results.Seurat.", z, "-", y, ".csv"), row.names = 1)
      colnames(results.Seurat) <- gsub("\\.", "_", colnames(results.Seurat))
    }
    
    # common celltypes from the ground truth results and scenario results
    missing.ct <- setdiff(colnames(spots.metadata), colnames(results.Seurat))
    common.ct <- intersect(colnames(spots.metadata), colnames(results.Seurat))
    
    results.Seurat <- results.Seurat[, common.ct]
    if (length(missing.ct)!=0) {
      for (c in 1:length(missing.ct)){
        results.Seurat <- cbind(results.Seurat, "_" = rep(0, dim(spots.metadata)[1]))
      }
    }
    colnames(results.Seurat) <- c(common.ct, missing.ct)
    results.Seurat <- results.Seurat[, order(colnames(results.Seurat))]
    
    list_ <- c(list_, list(results.Seurat))
  }
  
  if ("SPOTlight" %in% methods_) {
    results.SPOTlight <- read.csv(paste0(Method_Res, "SPOTlight/results.SPOTlight.", z, "-", y, ".csv"),
                                  row.names = 1)
    
    # common celltypes from the ground truth results and scenario results
    missing.ct <- setdiff(colnames(spots.metadata), colnames(results.SPOTlight))
    common.ct <- intersect(colnames(spots.metadata), colnames(results.SPOTlight))
    
    # creating a vector of zeroes in deconvolution results for missing celltype/s
    results.SPOTlight <- results.SPOTlight[, common.ct]
    if (length(missing.ct)!=0) {
      for (c in 1:length(missing.ct)){
        results.SPOTlight <- cbind(results.SPOTlight, "_" = rep(0, dim(spots.metadata)[1]))
      }
    }
    colnames(results.SPOTlight) <- c(common.ct, missing.ct)
    results.SPOTlight <- results.SPOTlight[, order(colnames(results.SPOTlight))]
    
    list_ <- c(list_, list(results.SPOTlight))
  }
  
  return(c(list_, list(spots.metadata)))
}

# calculate RMSE values among the spots
calculate_RMSE <- function(sc_num, st_num, results_path) {
  # sc_num: total number of sc reference dataset for the scenario
  # st_num: index of ST dataset
  # results_path: path to deconvolution method results for a scenario
  
  # reading deconvolution results
  Method_Res <- results_path
  
  # concatenated list of results from all the deconvolution methods
  if("cell2location" %in% methods_) {
    cell2loc_ <- list()
  }
  if("RCTD" %in% methods_) {
    rctd_ <- list()
  }
  if("CARD" %in% methods_) {
    card_ <- list()
  }
  if("SCDC" %in% methods_) {
    scdc_ <- list()
  }
  if("MuSiC" %in% methods_) {
    music_ <- list()
  }
  if("Stereoscope" %in% methods_) {
    stereoscope_ <- list()
  }
  if("Seurat" %in% methods_) {
    seurat_ <- list()
  }
  if("SPOTlight" %in% methods_) {
    spotlight_ <- list() 
  }
  gt_ <- list()
  
  
  aa <- 1
  # reading results for baseline scenario
  for (y in 1:sc_num) { # number of sc-ref data
    for (z in st_num:st_num) { # number of st-data
      
      r <- read_results(z, y, Method_Res)
      
      counter_ <- 1
      
      if("cell2location" %in% methods_) {
        cell2loc_[[aa]] <- r[[counter_]]
        counter_ <- counter_ + 1
      }
      if("RCTD" %in% methods_) {
        rctd_[[aa]] <- r[[counter_]]
        counter_ <- counter_ + 1
      }
      if("CARD" %in% methods_) {
        card_[[aa]] <- r[[counter_]]
        counter_ <- counter_ + 1
      }
      if("SCDC" %in% methods_) {
        scdc_[[aa]] <- r[[counter_]]
        counter_ <- counter_ + 1
      }
      if("MuSiC" %in% methods_) {
        music_[[aa]] <- r[[counter_]]
        counter_ <- counter_ + 1
      }
      if("Stereoscope" %in% methods_) {
        stereoscope_[[aa]] <- r[[counter_]]
        counter_ <- counter_ + 1
      }
      if("Seurat" %in% methods_) {
        seurat_[[aa]] <- r[[counter_]]
        counter_ <- counter_ + 1
      }
      if("SPOTlight" %in% methods_) {
        spotlight_[[aa]] <- r[[counter_]]
        counter_ <- counter_ + 1
      }
      gt_[[aa]] <- r[[counter_]]
      
      aa <- aa + 1
    }
  }
  
  
  # calculating JSD for each simulated data
  methodsResult <- list()
  
  if("cell2location" %in% methods_) {
    methodsResult <- c(methodsResult, list(cell2loc_))
  }
  if("RCTD" %in% methods_) {
    methodsResult <- c(methodsResult, list(rctd_))
  }
  if("CARD" %in% methods_) {
    methodsResult <- c(methodsResult, list(card_))
  }
  if("SCDC" %in% methods_) {
    methodsResult <- c(methodsResult, list(scdc_))
  }
  if("MuSiC" %in% methods_) {
    methodsResult <- c(methodsResult, list(music_))
  }
  if("Stereoscope" %in% methods_) {
    methodsResult <- c(methodsResult, list(stereoscope_))
  }
  if("Seurat" %in% methods_) {
    methodsResult <- c(methodsResult, list(seurat_))
  }
  if("SPOTlight" %in% methods_) {
    methodsResult <- c(methodsResult, list(spotlight_))
  }
  
  
  
  # calculating JSD between a spot composition in ground truth and prediction made
  # by a method
  rmseMat <- matrix(nrow = nrow(gt_[[1]]), ncol = length(methods_)) %>% data.frame()
  colnames(rmseMat) <- methods_
  
  rmseList <- list()
  
  # Remove notation in the entire R session using options() function
  # options(scipen = 100)
  for (m in 1:length(gt_)) {
    for (p in 1:ncol(rmseMat)) {
      for (q in 1:nrow(rmseMat)) {
        act_ <- as.numeric(gt_[[m]][q, ])
        pred_ <- as.numeric(methodsResult[[p]][[m]][q, ])
        
        if (any(is.na(pred_))) {
          rmseMat[q, p] <- NA
        } else {
          x <- rbind(act_, pred_)
          
          rmseMat[q, p] <- Metrics::rmse(act_, pred_)
        }
      }
    }
    rmseList[[m]] <- rmseMat
  }
  
  return(rmseList)
}


# list of JSD vectors calculated using all the reference dataset for a 
# particular scenario
r.list <- list()
counter_ <- 1

## baseline 
# reading baseline results
Method_Res <- paste0(decon_results, "rm0/")

rmseList <- calculate_RMSE(1, st_num, Method_Res)

# mean JSD value for a method across all the spots
r.list[[counter_]] <- rmseList[[1]]
counter_ <- counter_ + 1


## removal of 1 cell type
if ("rm1" %in% rm_scenarios) {
  Method_Res <- paste0(decon_results, "rm1/")
  
  rmseList <- calculate_RMSE(13, st_num, Method_Res)
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], rmseList[[6]],
                      rmseList[[7]], rmseList[[8]], rmseList[[9]],
                      rmseList[[10]], rmseList[[11]], rmseList[[12]],
                      rmseList[[13]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[counter_]] <- temp_array2
  counter_ <- counter_ + 1
}



## removal of 2 cell types
if ("rm2" %in% rm_scenarios) {
  Method_Res <- paste0(decon_results, "rm2/")
  
  rmseList <- calculate_RMSE(5, st_num, Method_Res)
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[counter_]] <- temp_array2
  counter_ <- counter_ + 1
}

## removal of 3 cell types
if ("rm3" %in% rm_scenarios) {
  Method_Res <- paste0(decon_results, "rm3/")
  
  rmseList <- calculate_RMSE(5, st_num, Method_Res)
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[counter_]] <- temp_array2
  counter_ <- counter_ + 1
}


## removal of 5 cell types
if ("rm5" %in% rm_scenarios) {
  Method_Res <- paste0(decon_results, "rm5/")
  
  rmseList <- calculate_RMSE(1, st_num, Method_Res)
  
  r.list[[counter_]] <- rmseList[[1]]
  counter_ <- counter_ + 1
}


## removal of 10 cell types
if ("rm10" %in% rm_scenarios) {
  Method_Res <- paste0(decon_results, "rm10/")
  
  rmseList <- calculate_RMSE(5, st_num, Method_Res)
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[counter_]] <- temp_array2
  counter_ <- counter_ + 1
}


## removal of 11 cell types
if ("rm11" %in% rm_scenarios) {
  Method_Res <- paste0(decon_results, "rm11/")
  
  rmseList <- calculate_RMSE(5, st_num, Method_Res)
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[counter_]] <- temp_array2
  counter_ <- counter_ + 1
}


saveRDS(r.list, paste0(Results, "RMSE.means.", st_num,"_ST.rds"))
