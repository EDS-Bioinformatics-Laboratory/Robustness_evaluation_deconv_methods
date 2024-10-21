################################################################################
# Funky plot for ranking performance of deconvolution methods
# cell type assignment for removal of 1 cell type from reference data
# JSD values for all the removal of cell type scenarios
# RMSE values for all the removal of cell type scenarios
# 
# https://funkyheatmap.dynverse.org
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################


#### Initialise environment ####
source("Init_env.R")

methods_ <- c("CARD",
              "Cell2Location",
              "MuSiC",
              "RCTD",
              "SCDC",
              "Seurat",
              "SPOTlight",
              "Stereoscope")



# normalize data with min-max scaling using functionalities in caret package
norm_min_max <- function(x) {
  # The "range" transformation scales the data to be within rangeBounds.
  process <- caret::preProcess(as.data.frame(x), method = c("range"))
  norm_scale <- predict(process, as.data.frame(x))
  return (norm_scale)
}

# JSD results: list of JSD mean vectors for each ST dataset
jsd.mean.list <- readRDS(paste0("../../4_Analysis_results/Results/",
                                "JSD.mean.4ST.RDS"))
jsd.mean.list <- jsd.mean.list[-3]
jsd.rank.list <- readRDS(paste0("../../4_Analysis_results/Results/",
                                "JSD.rank.4ST.RDS"))
jsd.rank.list <- jsd.rank.list[-3]


# RMSE results: list of RMSE mean vectors for each ST dataset
rmse.mean.list <- readRDS(paste0("../../4_Analysis_results/Results/",
                                 "RMSE.mean.4ST.RDS"))
rmse.mean.list <- rmse.mean.list[-3]
rmse.rank.list <- readRDS(paste0("../../4_Analysis_results/Results/",
                                 "RMSE.rank.4ST.RDS"))
rmse.rank.list <- rmse.rank.list[-3]


#### Funky plots ####

jmat.df <- data.frame()
for (q in 1:3) {
  jmat <- jsd.mean.list[[q]]
  jmat.df <- rbind(jmat.df, jmat)
}

# # average ranking based jsd values for all 3 ST
# jmat.df <- rbind(jsdRank = (colSums(jsd.rank.list[[1]]) + colSums(jsd.rank.list[[2]]) + colSums(jsd.rank.list[[3]]))/3,
#                  jmat.df)

rmat.df <- data.frame()
for (q in 1:3) {
  rmat <- rmse.mean.list[[q]]
  rmat.df <- rbind(rmat.df, rmat)
}

# # average ranking based rmse values for all 3 ST
# rmat.df <- rbind(rmseRank = (colSums(rmse.rank.list[[1]]) + colSums(rmse.rank.list[[2]]) + colSums(rmse.rank.list[[3]]))/3,
#                  rmat.df)

jr.mat.df <- rbind(jmat.df, rmat.df)

# lower jsd and rmse values refers to higher ranking and thus subtracted from 1
jr.mat.df <- 1 - jr.mat.df


# average ranking based jsd and rmse values for all 3 ST
# jr.mat.df <- rbind(jr.mat.df, jsdRank = (colSums(jsd.rank.list[[1]]) + colSums(jsd.rank.list[[2]]) + colSums(jsd.rank.list[[3]]))/3)
# jr.mat.df <- rbind(jr.mat.df, rmseRank = (colSums(rmse.rank.list[[1]]) + colSums(rmse.rank.list[[2]]) + colSums(rmse.rank.list[[3]]))/3)

jr.mat.df <- rbind(jr.mat.df, jsdRank1 = colMeans(jsd.rank.list[[1]]))
jr.mat.df <- rbind(jr.mat.df, rmseRank1= colMeans(rmse.rank.list[[1]]))

jr.mat.df <- rbind(jr.mat.df, jsdRank2 = colMeans(jsd.rank.list[[2]]))
jr.mat.df <- rbind(jr.mat.df, rmseRank2= colMeans(rmse.rank.list[[2]]))

jr.mat.df <- rbind(jr.mat.df, jsdRank3 = colMeans(jsd.rank.list[[3]]))
jr.mat.df <- rbind(jr.mat.df, rmseRank3= colMeans(rmse.rank.list[[3]]))


# adding the overall ranking of the methods.
# ranking for jsd and rmse has been aggregated to represent one value
jr.mat.df <- rbind(jr.mat.df, overallRank1 = (jr.mat.df["rmseRank1",] + jr.mat.df["jsdRank1",]) / 2)
jr.mat.df <- rbind(jr.mat.df, overallRank2 = (jr.mat.df["rmseRank2",] + jr.mat.df["jsdRank2",]) / 2)
jr.mat.df <- rbind(jr.mat.df, overallRank3 = (jr.mat.df["rmseRank3",] + jr.mat.df["jsdRank3",]) / 2)


jr.mat.df2 <- data.frame(t(jr.mat.df))
colnames(jr.mat.df2) <- rownames(jr.mat.df)
jr.mat.df2 <- jr.mat.df2 %>% rownames_to_column("id")

# column info for ST dataset 1
column_info1 <- tribble(
  ~id,	~group,	~name,	~geom,	~palette,	~options,
  "id",	"Methods",	"",	"text",	"palette",	list(hjust = 0, width = 6),
  
  "overallRank1", "Overall1", "", "bar", "palette0", lst(width = 4),
  
  "Baseline: 13 CT present",	"Accuracy1",	"J",	"bar",	"palette1",	lst(width = 2.5),
  "Baseline: 13 CT present3",	"Accuracy1",	"R",	"bar",	"palette1",	lst(width = 2.5),
  
  "Scenario 1: 12 CT present",	"RobustnessDataset11",	"J",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 1: 12 CT present3",	"RobustnessDataset11",	"R",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 2: 11 CT present",	"RobustnessDataset12",	"J",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 2: 11 CT present3",	"RobustnessDataset12",	"R",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 3: 10 CT present",	"RobustnessDataset13",	"J",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 3: 10 CT present3",	"RobustnessDataset13",	"R",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 4: 8 CT present",	"RobustnessDataset14",	"J",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 4: 8 CT present3",	"RobustnessDataset14",	"R",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 5: 3 CT present",	"RobustnessDataset15",	"J",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 5: 3 CT present3",	"RobustnessDataset15",	"R",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 6: 2 CT present",	"RobustnessDataset16",	"J",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 6: 2 CT present3",	"RobustnessDataset16",	"R",	"funkyrect",	"palette2",	lst(width = 1)
)

# column info for ST dataset 2
column_info2 <- tribble(
  ~id,	~group,	~name,	~geom,	~palette,	~options,
  "id",	"Methods",	"",	"text",	"palette",	list(hjust = 0, width = 6),
  
  "overallRank2", "Overall2", "", "bar", "palette0", lst(width = 4),
  
  "Baseline: 13 CT present1",	"Accuracy2",	"JSD",	"bar",	"palette1",	lst(width = 2.5),
  "Baseline: 13 CT present11",	"Accuracy2",	"RMSE",	"bar",	"palette1",	lst(width = 2.5),
  
  "Scenario 1: 12 CT present1",	"RobustnessDataset21",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 1: 12 CT present11",	"RobustnessDataset21",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 2: 11 CT present1",	"RobustnessDataset22",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 2: 11 CT present11",	"RobustnessDataset22",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 3: 10 CT present1",	"RobustnessDataset23",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 3: 10 CT present11",	"RobustnessDataset23",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 4: 8 CT present1",	"RobustnessDataset24",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 4: 8 CT present11",	"RobustnessDataset24",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 5: 3 CT present1",	"RobustnessDataset25",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 5: 3 CT present11",	"RobustnessDataset25",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 6: 2 CT present1",	"RobustnessDataset26",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 6: 2 CT present11",	"RobustnessDataset26",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1)
)

# column info for ST dataset 3
column_info3 <- tribble(
  ~id,	~group,	~name,	~geom,	~palette,	~options,
  "id",	"Methods",	"",	"text",	"palette",	list(hjust = 0, width = 6),
  
  "overallRank3", "Overall3", "", "bar", "palette0", lst(width = 4),
  
  "Baseline: 13 CT present2",	"Accuracy3",	"JSD",	"bar",	"palette1",	lst(width = 2.5),
  "Baseline: 13 CT present21",	"Accuracy3",	"RMSE",	"bar",	"palette1",	lst(width = 2.5),
  
  "Scenario 1: 12 CT present2",	"RobustnessDataset31",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 1: 12 CT present21",	"RobustnessDataset31",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 2: 11 CT present2",	"RobustnessDataset32",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 2: 11 CT present21",	"RobustnessDataset32",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 3: 10 CT present2",	"RobustnessDataset33",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 3: 10 CT present21",	"RobustnessDataset33",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 4: 8 CT present2",	"RobustnessDataset34",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 4: 8 CT present21",	"RobustnessDataset34",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 5: 3 CT present2",	"RobustnessDataset35",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 5: 3 CT present21",	"RobustnessDataset35",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 6: 2 CT present2",	"RobustnessDataset36",	"JSD",	"funkyrect",	"palette2",	lst(width = 1),
  "Scenario 6: 2 CT present21",	"RobustnessDataset36",	"RMSE",	"funkyrect",	"palette2",	lst(width = 1)
)

# column groups for ST dataset 1
col_groups1 <- tibble::tibble(col1 = c("Method",
                                       "Overall",
                                       "Accuracy",
                                       "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch"), #"Robustness for ST1",
                             Category = c("",
                                          "",
                                          "0 CT",
                                          "1 CT", "2 CT", "3 CT", "5 CT", "10 CT", "11 CT"),
                             group = c("Methods",
                                       "Overall1",
                                       "Accuracy1",
                                       "RobustnessDataset11", "RobustnessDataset12", "RobustnessDataset13", "RobustnessDataset14", "RobustnessDataset15", "RobustnessDataset16"),
                             palette = c("palette",
                                         "palette0",
                                         "palette1",
                                         "palette2", "palette2", "palette2", "palette2", "palette2", "palette2"),
)

# column groups for ST dataset 2
col_groups2 <- tibble::tibble(col1 = c("Method",
                                       "Overall",
                                       "Accuracy",
                                       "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch"), #"Robustness for ST1",
                              Category = c("",
                                           "",
                                           "0 CT",
                                           "1 CT", "2 CT", "3 CT", "5 CT", "10 CT", "11 CT"),
                              group = c("Methods",
                                        "Overall2",
                                        "Accuracy2",
                                        "RobustnessDataset21", "RobustnessDataset22", "RobustnessDataset23", "RobustnessDataset24", "RobustnessDataset25", "RobustnessDataset26"),
                              palette = c("palette",
                                          "palette0",
                                          "palette1",
                                          "palette2", "palette2", "palette2", "palette2", "palette2", "palette2"),
)

# column groups for ST dataset 3
col_groups3 <- tibble::tibble(col1 = c("Method",
                                       "Overall",
                                       "Accuracy",
                                       "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch", "Robustness during cell type mismatch"), #"Robustness for ST1",
                              Category = c("",
                                           "",
                                           "0 CT",
                                           "1 CT", "2 CT", "3 CT", "5 CT", "10 CT", "11 CT"),
                              group = c("Methods",
                                        "Overall3",
                                        "Accuracy3",
                                        "RobustnessDataset31", "RobustnessDataset32", "RobustnessDataset33", "RobustnessDataset34", "RobustnessDataset35", "RobustnessDataset36"),
                              palette = c("palette",
                                          "palette0",
                                          "palette1",
                                          "palette2", "palette2", "palette2", "palette2", "palette2", "palette2"),
)

row_info <- tribble(
  ~group, ~id,
  "ST methods", "Cell2Location",
  "ST methods", "RCTD",
  "ST methods", "CARD",
  "ST methods", "Stereoscope",
  "ST methods", "Seurat",
  "ST methods", "SPOTlight",
  "Bulk methods", "SCDC",
  "Bulk methods", "MuSiC"
)

funky1 <- funkyheatmap::funky_heatmap(norm_min_max(jr.mat.df2),
                        column_info = column_info1,
                        column_groups = col_groups1,
                        row_info = row_info,
                        scale_column = F,
                        expand = c(xmin = 2, xmax = 2, ymin = 10, ymax = 10),
                        # col_annot_angle = 10,
                        col_annot_offset = 2)
funky1

funky2 <- funkyheatmap::funky_heatmap(norm_min_max(jr.mat.df2),
                        column_info = column_info2,
                        column_groups = col_groups2,
                        row_info = row_info,
                        scale_column = F,
                        expand = c(xmin = 2, xmax = 2, ymin = 10, ymax = 10),
                        # col_annot_angle = 80,
                        col_annot_offset = 3)
funky2

funky3 <- funkyheatmap::funky_heatmap(norm_min_max(jr.mat.df2),
                        column_info = column_info3,
                        column_groups = col_groups3,
                        row_info = row_info,
                        scale_column = F,
                        expand = c(xmin = 2, xmax = 2, ymin = 10, ymax = 10),
                        # col_annot_angle = 80,
                        col_annot_offset = 3)
funky3


png(file = paste0(Results, "Fig_funky_rank_1.1.png"), 
    res = 450, width = funky1$width, height = funky1$height, units = "in")
funky1
dev.off()

png(file = paste0(Results, "Fig_funky_rank_2.png"), 
    res = 450, width = funky1$width, height = funky1$height, units = "in")
funky2
dev.off()

png(file = paste0(Results, "Fig_funky_rank_3.png"), 
    res = 450, width = funky1$width, height = funky1$height, units = "in")
funky3
dev.off()