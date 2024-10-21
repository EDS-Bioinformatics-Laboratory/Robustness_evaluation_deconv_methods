################################################################################
# This script calculates the difference in the RMSE values for all the removal
# scenarios (wtih cell type mismatch) compared to the baseline scenario (with
# no cell type mismatch)
# 
# 
# We used a vector (length = number of spots in ST data) of calculated RMSE values
# Incase of multiple RMSE vectors (depending on number of single cell reference
# datasets available for a particular removal scenario) for a single method and
# scenario, we considered the mean vector for calculation. 
# 
# 
# RMSE.mean.list.4ST.RDS (also for median, ranking, entire RMSE values as a list)
# Refers to 4 list, one for each ST dataset
# Each list has 7 sub-lists, one for each scenario (baseline + 6 removal scenario)
# Each sub-list is a matrix with columns as methods and rows as spots
# [i,j] is the mean JSD value of spot i for method j over multiple single cell
# reference dataset
# 
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################



# Initialize environment
source("Init_env.R")

methods_ <- c("Cell2Location",
              "RCTD",
              "CARD",
              "SCDC",
              "MuSiC",
              "Stereoscope",
              "Seurat",
              "SPOTlight")

st_data_path <- "../../2_Simulating_ST_data/Results/Spatial_Data/"

# reading the deconvolution method results
read_results <- function(z, y, Method_Res) {
  # z refers to ST dataset number
  # y refers to sc reference dataset number for the removal scenario
  # Method_Res refers to path of deconvolution methods results
  
  spots.metadata <- readRDS(paste0(st_data_path, "metadata.", z, ".rds"))
  
  results.RCTD <- read.csv(paste0(Method_Res, "RCTD/results.RCTD.", z, "-", y, ".csv"),
                           row.names = 1)
  
  results.SPOTlight <- read.csv(paste0(Method_Res, "SPOTlight/results.SPOTlight.", z, "-", y, ".csv"),
                                row.names = 1)
  
  if (file.exists(paste0(Method_Res, "Seurat/results.Seurat.", z, "-", y, ".csv")) != T) {
    results.Seurat <- matrix(NA, nrow = dim(results.SPOTlight)[1], ncol = dim(results.SPOTlight)[2]) %>% data.frame()
    colnames(results.Seurat) <- colnames(results.SPOTlight)
    rownames(results.Seurat) <- rownames(results.SPOTlight)
    message("Seurat ",z,y)
  } else {
    results.Seurat <- read.csv(paste0(Method_Res, "Seurat/results.Seurat.", z, "-", y, ".csv"), row.names = 1)
    colnames(results.Seurat) <- colnames(results.SPOTlight)
  }
  
  results.Stereoscope <- 
    utils::read.csv(paste0(Method_Res, "Stereoscope/Stereoscope_result.", z, "-", y, ".csv"), row.names = 1)
  
  results.Cell2Loc <- 
    utils::read.csv(paste0(Method_Res, "Cell2Location/Cell2Location.", z, "-", y, ".csv"), row.names = 1)
  colnames(results.Cell2Loc) <- colnames(results.Stereoscope)
  row_sum <- rowSums(results.Cell2Loc)
  results.Cell2Loc <- results.Cell2Loc %>% mutate_all(~ . /row_sum)
  
  if (file.exists(paste0(Method_Res, "CARD/results.CARD.", z, "-", y, ".csv")) != T) {
    results.CARD <- matrix(NA, nrow = dim(results.SPOTlight)[1], ncol = dim(results.SPOTlight)[2]) %>% data.frame()
    colnames(results.CARD) <- colnames(results.SPOTlight)
    rownames(results.CARD) <- rownames(results.SPOTlight)
    message("CARD ",z,y)
  } else {
    results.CARD <- read.csv(paste0(Method_Res, "CARD/results.CARD.", z, "-", y, ".csv"), row.names = 1)
  }
  
  results.SCDC <- read.csv(paste0(Method_Res, "SCDC/results.SCDC.", z, "-", y, ".csv"), row.names = 1)
  
  results.MuSiC <- read.csv(paste0(Method_Res, "MuSiC/results.MuSiC.", z, "-", y, ".csv"), row.names = 1)
  
  
  ## adjusting the scenario results as per the baseline results 
  # common celltypes from the ground truth results and scenario results
  missing.ct <- setdiff(colnames(spots.metadata), colnames(results.RCTD))
  common.ct <- intersect(colnames(spots.metadata), colnames(results.RCTD))
  
  # creating a zero vector in method prediction results for missing celltype/s
  # to compare with the baseline predictions
  results.RCTD <- results.RCTD[, common.ct]
  if (length(missing.ct!=0)) {
    for (c in 1:length(missing.ct)){
      results.RCTD <- cbind(results.RCTD, "_" = rep(0, dim(spots.metadata)[1]))
    }
  }
  colnames(results.RCTD) <- c(common.ct, missing.ct)
  results.RCTD <- results.RCTD[, order(colnames(results.RCTD))]
  
  results.SPOTlight <- results.SPOTlight[, common.ct]
  if (length(missing.ct!=0)) {
    for (c in 1:length(missing.ct)){
      results.SPOTlight <- cbind(results.SPOTlight, "_" = rep(0, dim(spots.metadata)[1]))
    }
  }
  colnames(results.SPOTlight) <- c(common.ct, missing.ct)
  results.SPOTlight <- results.SPOTlight[, order(colnames(results.SPOTlight))]
  
  results.Seurat <- results.Seurat[, common.ct]
  if (length(missing.ct!=0)) {
    for (c in 1:length(missing.ct)){
      results.Seurat <- cbind(results.Seurat, "_" = rep(0, dim(spots.metadata)[1]))
    }
  }
  colnames(results.Seurat) <- c(common.ct, missing.ct)
  results.Seurat <- results.Seurat[, order(colnames(results.Seurat))]
  
  results.Stereoscope <- results.Stereoscope[, common.ct]
  if (length(missing.ct!=0)) {
    for (c in 1:length(missing.ct)){
      results.Stereoscope <- cbind(results.Stereoscope, "_" = rep(0, dim(spots.metadata)[1]))
    }
  }
  colnames(results.Stereoscope) <- c(common.ct, missing.ct)
  results.Stereoscope <- results.Stereoscope[, order(colnames(results.Stereoscope))]
  
  results.Cell2Loc <- results.Cell2Loc[, common.ct]
  if (length(missing.ct!=0)) {
    for (c in 1:length(missing.ct)){
      results.Cell2Loc <- cbind(results.Cell2Loc, "_" = rep(0, dim(spots.metadata)[1]))
    }
  }
  colnames(results.Cell2Loc) <- c(common.ct, missing.ct)
  results.Cell2Loc <- results.Cell2Loc[, order(colnames(results.Cell2Loc))]
  
  results.CARD <- results.CARD[, common.ct]
  if (length(missing.ct!=0)) {
    for (c in 1:length(missing.ct)){
      results.CARD <- cbind(results.CARD, "_" = rep(0, dim(spots.metadata)[1]))
    }
  }
  colnames(results.CARD) <- c(common.ct, missing.ct)
  results.CARD <- results.CARD[, order(colnames(results.CARD))]
  
  results.SCDC <- results.SCDC[, common.ct]
  if (length(missing.ct!=0)) {
    for (c in 1:length(missing.ct)){
      results.SCDC <- cbind(results.SCDC, "_" = rep(0, dim(spots.metadata)[1]))
    }
  }
  colnames(results.SCDC) <- c(common.ct, missing.ct)
  results.SCDC <- results.SCDC[, order(colnames(results.SCDC))]
  
  results.MuSiC <- results.MuSiC[, common.ct]
  if (length(missing.ct!=0)) {
    for (c in 1:length(missing.ct)){
      results.MuSiC <- cbind(results.MuSiC, "_" = rep(0, dim(spots.metadata)[1]))
    }
  }
  colnames(results.MuSiC) <- c(common.ct, missing.ct)
  results.MuSiC <- results.MuSiC[, order(colnames(results.MuSiC))]
  
  return(list(results.RCTD, results.SPOTlight, results.Seurat, results.CARD,
              results.SCDC, results.MuSiC, results.Stereoscope, results.Cell2Loc,
              spots.metadata))
}

# calculate RMSE values among the spots
calculate_RMSE <- function(sc_num, st_num, results_path) {
  # sc_num: number of sc reference dataset
  # st_num: index of ST dataset
  # results_path: path to deconvolution method results for a scenario
  
  # reading deconvolution results
  Method_Res <- results_path
  
  # concatenated list of results from all the deconvolution methods
  rctd_ <- list()
  spotlight_ <- list() 
  seurat_ <- list()
  stereoscope_ <- list()
  cell2loc_ <- list()
  card_ <- list()
  scdc_ <- list()
  music_ <- list()
  gt_ <- list()
  
  aa <- 1
  # reading results for baseline scenario
  for (y in 1:sc_num) { # number of sc-ref data
    for (z in st_num:st_num) { # number of st-data
      
      r <- read_results(z, y, Method_Res)
      
      rctd_[[aa]] <- r[[1]]
      spotlight_[[aa]] <- r[[2]]
      seurat_[[aa]] <- r[[3]]
      card_[[aa]] <- r[[4]]
      scdc_[[aa]] <- r[[5]]
      music_[[aa]] <- r[[6]]
      stereoscope_[[aa]] <- r[[7]]
      cell2loc_[[aa]] <- r[[8]]
      gt_[[aa]] <- r[[9]]
      
      aa <- aa + 1
    }
  }
  
  
  # calculating RMSE for each simulated data
  methodsResult <- list(cell2loc_,
                        rctd_,
                        card_,
                        scdc_,
                        music_,
                        stereoscope_,
                        seurat_,
                        spotlight_)
  
  # calculating RMSE between a spot composition in ground truth and prediction made
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


mean.RMSE.list <- list() # list of mean RMSE values for all methods (1 value for all spots)
med.RMSE.list <- list() # list of mean RMSE values for all methods
RMSE.rank.list <- list() # list of ranks based on RMSE values for all methods
rmse.lists <- list() # list of RMSE values for all methods (1 value for each spots)


for (st in 1:4) {
  # data frame of RMSE means for all the methods (mean is of RMSE vectors 
  # calculated for reference dataset for a particular)
  mean.RMSE <- matrix(nrow = 7, ncol = length(methods_)) %>% data.frame()
  colnames(mean.RMSE) <- methods_
  rownames(mean.RMSE) <- c("Baseline: 13 CT present", "Scenario 1: 12 CT present",
                           "Scenario 2: 11 CT present", "Scenario 3: 10 CT present",
                           "Scenario 4: 8 CT present", "Scenario 5: 3 CT present",
                           "Scenario 6: 2 CT present")
  
  # data frame of RMSE medians for all the methods (median is of RMSE vectors 
  # calculated for reference dataset for a particular)
  med.RMSE <- matrix(nrow = 7, ncol = length(methods_)) %>% data.frame()
  colnames(med.RMSE) <- methods_
  rownames(med.RMSE) <- c("Baseline: 13 CT present", "Scenario 1: 12 CT present",
                          "Scenario 2: 11 CT present", "Scenario 3: 10 CT present",
                          "Scenario 4: 8 CT present", "Scenario 5: 3 CT present",
                          "Scenario 6: 2 CT present")
  
  # dataframe to rank the median RMSE value for a scenario (median is of RMSE
  # vectors calculated for reference dataset for a particular)
  RMSE.rank <- matrix(nrow = 7, ncol = length(methods_)) %>% data.frame()
  colnames(RMSE.rank) <- methods_
  rownames(RMSE.rank) <- c("Baseline: 13 CT present", "Scenario 1: 12 CT present",
                           "Scenario 2: 11 CT present", "Scenario 3: 10 CT present",
                           "Scenario 4: 8 CT present", "Scenario 5: 3 CT present",
                           "Scenario 6: 2 CT present")
  
  # list of RMSE vectors calculated using all the reference dataset for a 
  # particular scenario
  r.list <- list()
  
  
  ## baseline 
  # reading baseline results
  Method_Res <- "../../3_ST_methods/Results/rm0/"
  
  rmseList <- calculate_RMSE(1, st, Method_Res)
  
  # mean RMSE value for a method across all the spots
  mean.RMSE[1, ] <- colMeans(rmseList[[1]])
  med.RMSE[1, ] <- apply(rmseList[[1]], 2, median)
  RMSE.rank[1, ] <- rank(-med.RMSE[1, ])
  r.list[[1]] <- rmseList[[1]]
  
  
  ## removal of 1 cell type
  # reading results for removal of one cell type from single cell ref data
  Method_Res <- "../../3_ST_methods/Results/rm1/"
  
  rmseList <- calculate_RMSE(13, st, Method_Res)
  
  # getting the column means for each scenario first in a dataframe, followed by
  # column mean to get the mean RMSE value for a method across all scenarios
  rm.mean <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  rm.median <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  for (i in 1:length(rmseList)) {
    t <- colMeans(rmseList[[i]])
    t2 <- apply(rmseList[[1]], 2, median)
    rm.mean[i, ] <- t
    rm.median[i, ] <- t2
  }
  mean.RMSE[2, ] <- colMeans(rm.mean)
  med.RMSE[2, ] <- colMeans(rm.median)
  RMSE.rank[2, ] <- rank(-med.RMSE[2, ])
  
  temp_array <- abind::abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], rmseList[[6]],
                      rmseList[[7]], rmseList[[8]], rmseList[[9]],
                      rmseList[[10]], rmseList[[11]], rmseList[[12]],
                      rmseList[[13]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[2]] <- temp_array2
  
  
  ## removal of 2 cell types
  # reading results for removal of two cell types from single cell ref data
  Method_Res <- "../../3_ST_methods/Results/rm2/"
  
  rmseList <- calculate_RMSE(5, st, Method_Res)
  
  # getting the column means for each scenario first in a dataframe, followed by
  # column mean to get the mean RMSE value for a method across all scenarios
  rm.mean <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  rm.median <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  for (i in 1:length(rmseList)) {
    t <- colMeans(rmseList[[i]])
    t2 <- apply(rmseList[[1]], 2, median)
    rm.median[i, ] <- t2
    rm.mean[i, ] <- t
  }
  mean.RMSE[3, ] <- colMeans(rm.mean)
  med.RMSE[3, ] <- colMeans(rm.median)
  RMSE.rank[3, ] <- rank(-med.RMSE[3, ])
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[3]] <- temp_array2
  
  
  ## removal of 3 cell types
  # reading results for removal of three cell types from single cell ref data
  Method_Res <- "../../3_ST_methods/Results/rm3/"
  
  rmseList <- calculate_RMSE(5, st, Method_Res)
  
  # getting the column means for each scenario first in a dataframe, followed by
  # column mean to get the mean RMSE value for a method across all scenarios
  rm.mean <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  rm.median <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  for (i in 1:length(rmseList)) {
    t <- colMeans(rmseList[[i]])
    rm.mean[i, ] <- t
    t2 <- apply(rmseList[[1]], 2, median)
    rm.median[i, ] <- t2
  }
  med.RMSE[4, ] <- colMeans(rm.median)
  mean.RMSE[4, ] <- colMeans(rm.mean)
  RMSE.rank[4, ] <- rank(-med.RMSE[4, ])
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[4]] <- temp_array2
  
  
  ## removal of 5 cell types
  # reading results for removal of five cell types from single cell ref data
  Method_Res <- "../../3_ST_methods/Results/rm5/"
  
  rmseList <- calculate_RMSE(1, st, Method_Res)
  
  # getting the column means for each scenario first in a dataframe, followed by
  # column mean to get the mean RMSE value for a method across all scenarios
  rm.mean <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  rm.median <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  for (i in 1:length(rmseList)) {
    t <- colMeans(rmseList[[i]])
    rm.mean[i, ] <- t
    t2 <- apply(rmseList[[1]], 2, median)
    rm.median[i, ] <- t2
  }
  med.RMSE[5, ] <- colMeans(rm.median)
  mean.RMSE[5, ] <- colMeans(rm.mean)
  RMSE.rank[5, ] <- rank(-med.RMSE[5, ])
  r.list[[5]] <- rmseList[[1]]
  
  
  ## removal of 10 cell types
  # reading results for removal of ten cell types from single cell ref data
  Method_Res <- "../../3_ST_methods/Results/rm10/"
  
  rmseList <- calculate_RMSE(5, st, Method_Res)
  
  # getting the column means for each scenario first in a dataframe, followed by
  # column mean to get the mean RMSE value for a method across all scenarios
  rm.mean <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  rm.median <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  for (i in 1:length(rmseList)) {
    t <- colMeans(rmseList[[i]])
    rm.mean[i, ] <- t
    t2 <- apply(rmseList[[1]], 2, median)
    rm.median[i, ] <- t2
  }
  med.RMSE[6, ] <- colMeans(rm.median)
  mean.RMSE[6, ] <- colMeans(rm.mean, na.rm = T)
  RMSE.rank[6, ] <- rank(-med.RMSE[6, ])
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[6]] <- temp_array2
  
  
  ## removal of 11 cell types
  # reading results for removal of eleven cell types from single cell ref data
  Method_Res <- "../../3_ST_methods/Results/rm11/"
  
  rmseList <- calculate_RMSE(5, st, Method_Res)
  
  # getting the column means for each scenario first in a dataframe, followed by
  # column mean to get the mean RMSE value for a method across all scenarios
  rm.mean <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  rm.median <- matrix(nrow = length(rmseList), ncol = length(methods_)) %>% data.frame()
  for (i in 1:length(rmseList)) {
    t <- colMeans(rmseList[[i]])
    # message(i)
    # print(t)
    rm.mean[i, ] <- t
    t2 <- apply(rmseList[[1]], 2, median)
    rm.median[i, ] <- t2
  }
  med.RMSE[7, ] <- colMeans(rm.median)
  mean.RMSE[7, ] <- colMeans(rm.mean, na.rm = T)
  ff <- as.vector(med.RMSE[7, ])
  ff[is.na(ff)] <- 1
  ff
  RMSE.rank[7, ] <- rank(-ff)
  
  temp_array <- abind(rmseList[[1]], rmseList[[2]], rmseList[[3]],
                      rmseList[[4]], rmseList[[5]], along=3)
  temp_array2 <- as.data.frame(apply(temp_array, 1:2, mean))
  r.list[[7]] <- temp_array2
  
  
  # storing lists
  mean.RMSE.list[[st]] <- mean.RMSE
  med.RMSE.list[[st]] <- med.RMSE
  RMSE.rank.list[[st]] <- RMSE.rank
  rmse.lists[[st]] <- r.list
}


saveRDS(mean.RMSE.list, paste0(Results, "RMSE.mean.4ST.RDS"))
saveRDS(RMSE.rank.list, paste0(Results, "RMSE.rank.4ST.RDS"))
saveRDS(med.RMSE.list, paste0(Results, "RMSE.median.4ST.RDS"))
saveRDS(rmse.lists, paste0(Results, "RMSE.mean.list.4ST.RDS"))


