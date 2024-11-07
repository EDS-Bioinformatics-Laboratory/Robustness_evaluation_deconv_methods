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

args <- commandArgs(trailingOnly = TRUE)
st_num <- args[1]

methods_ <- c("Cell2Location",
              "RCTD",
              "CARD",
              "SCDC",
              "MuSiC",
              "Stereoscope",
              "Seurat",
              "SPOTlight")



# normalize data with min-max scaling using functionalities in caret package
norm_min_max <- function(x) {
  # The "range" transformation scales the data to be within rangeBounds.
  process <- caret::preProcess(as.data.frame(x), method = c("range"))
  norm_scale <- predict(process, as.data.frame(x))
  return (norm_scale)
}

# JSD/RMSE results: list of JSD?RMSE mean vectors for the specified ST dataset

jsd.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                            "JSD.means.", st_num, "_ST.rds"))
rmse.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                             "RMSE.means.", st_num, "_ST.rds"))


# we calculate mean JSD for a removal scenario over all the spots and later
# substract it from 1, since lower the jsd value higher the rank
jsd.mat <- as.data.frame(1- do.call(rbind, lapply(jsd.lists, colMeans)))
rownames(jsd.mat) <- c("0ctj", "1ctj", "2ctj", "3ctj", "5ctj", "10ctj", "11ctj")

# replace NA
jsd.mat[is.na(jsd.mat)] <- 0

# we ranks the methods for each removal scenario based on jsd values and then
# take the mean over all scenarios to get a final rank. Later this ranking is
# normalised higher the rank value, better the methods performed overall
jsd.rank <- t(as.data.frame(((apply(jsd.mat, 1, rank)))))
# jsd.rank <- as.data.frame((rowMeans(apply(jsd.mat, 1, rank))))
# colnames(jsd.rank) <- "ranks.jsd"

# jsd.mat <- rbind(jsd.mat, t(jsd.rank))


# we calculate mean RMSE for a removal scenario over all the spots and later
# substract it from 1, since lower the jsd value higher the rank
rmse.mat <- as.data.frame(1- do.call(rbind, lapply(rmse.lists, colMeans)))
rownames(rmse.mat) <- c("0ctr", "1ctr", "2ctr", "3ctr", "5ctr", "10ctr", "11ctr")

# replace Na
rmse.mat[is.na(rmse.mat)] <- 0

# we ranks the methods for each removal scenario based on rmse values and then
# take the mean over all scenarios to get a final rank. Later this ranking is
# normalised higher the rank value, better the methods performed overall
rmse.rank <- t(as.data.frame(((apply(rmse.mat, 1, rank)))))
# rmse.rank <- as.data.frame((rowMeans(apply(rmse.mat, 1, rank))))
# colnames(rmse.rank) <- "ranks.rmse"

# rmse.mat <- rbind(rmse.mat, t(rmse.rank))

tmp <- as.data.frame(t(rbind(jsd.rank, rmse.rank)))

overall.rank <- as.data.frame((colMeans(jsd.rank) + colMeans(rmse.rank)) / 2)
colnames(overall.rank) <- "overall.rank"

tmp <- cbind(tmp, overall.rank) %>%
  rownames_to_column("id")


# column info for ST dataset 1
column_info <- tidyr::tribble(
  ~id,	~group,	~name,	~geom,	~options,
  "id",	"Methods",	"",	"text",	list(palette = "methods", hjust = 0, width = 6, fontface = "italic", size=5),
  
  "overall.rank", "Overall1", "", "bar", list(palette = "ranks", width = 4.5),
  
  "0ctj",	"Accuracy1",	"JSD",	"bar",	list(palette = "accuracy", width = 3),
  "0ctr",	"Accuracy1",	"RMSE",	"bar",	list(palette = "accuracy", width = 3),
  
  "1ctj",	"RobustnessDataset11",	"JSD",	"funkyrect",	list(palette = "rm", width = 1.2),
  "1ctr",	"RobustnessDataset11",	"RMSE",	"funkyrect",	list(palette = "rm", width = 1.2),
  "2ctj",	"RobustnessDataset12",	"JSD",	"funkyrect",	list(palette = "rm", width = 1.2),
  "2ctr",	"RobustnessDataset12",	"RMSE",	"funkyrect",	list(palette = "rm", width = 1.2),
  "3ctj",	"RobustnessDataset13",	"JSD",	"funkyrect",	list(palette = "rm", width = 1.2),
  "3ctr",	"RobustnessDataset13",	"RMSE",	"funkyrect",	list(palette = "rm", width = 1.2),
  "5ctj",	"RobustnessDataset14",	"JSD",	"funkyrect",	list(palette = "rm", width = 1.2),
  "5ctr",	"RobustnessDataset14",	"RMSE",	"funkyrect",	list(palette = "rm", width = 1.2),
  "10ctj",	"RobustnessDataset15",	"JSD",	"funkyrect",	list(palette = "rm", width = 1.2),
  "10ctr",	"RobustnessDataset15",	"RMSE",	"funkyrect",	list(palette = "rm", width = 1.2),
  "11ctj",	"RobustnessDataset16",	"JSD",	"funkyrect",	list(palette = "rm", width = 1.2),
  "11ctr",	"RobustnessDataset16",	"RMSE",	"funkyrect",	list(palette = "rm", width = 1.2)
)


# column groups for ST dataset
col_groups <- tibble::tibble(col1 = c("Methods",
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
                             palette = c("methods",
                                         "ranks",
                                         "accuracy",
                                         "rm", "rm", "rm", "rm", "rm", "rm"),
)


row_info <- tidyr::tribble(
  ~group, ~id,
  "ST-based", "Cell2Location",
  "ST-based", "RCTD",
  "ST-based", "CARD",
  "ST-based", "Stereoscope",
  "ST-based", "Seurat",
  "ST-based", "SPOTlight",
  "Bulk rna-based", "SCDC",
  "Bulk rna-based", "MuSiC"
)

palettes <- list(
  methods = "Greys",
  ranks = (RColorBrewer::brewer.pal(9, "Blues")),
  accuracy = (RColorBrewer::brewer.pal(9, "Greens")),
  rm = (RColorBrewer::brewer.pal(9, "Reds"))
)

legends <- list(
  list(
    title = "Overall ranking",
    palette = "ranks",
    geom = "funkyrect",
    labels = c(" 8", "", "", "", "", "", "", "1"),
    size = c(.25, .375, .5, .6, .7, .8, .9, 1)
  ),
  list(
    title = "Baseline ranking",
    palette = "accuracy",
    geom = "funkyrect",
    labels = c(" 8", "", "", "", "", "", "", "1"),
    size = c(.25, .375, .5, .6, .7, .8, .9, 1)
  ),
  list(
    title = "Removal ranking",
    palette = "rm",
    geom = "funkyrect",
    labels = c(" 8", "", "", "", "", "", "", "1"),
    size = c(.25, .375, .5, .6, .7, .8, .9, 1)
  )
)


funky <- funkyheatmap::funky_heatmap(as.data.frame(norm_min_max(tmp)),
                                     column_info = column_info,
                                     column_groups = col_groups,
                                     row_info = row_info,
                                     palettes = palettes,
                                     legends = legends,
                                     scale_column = F,
                                     position_args = funkyheatmap::position_arguments(
                                       row_height = 1.2,
                                       row_space = 0.3,
                                       row_bigspace = 1.5,
                                       col_width = 1.2,
                                       col_space = 0.1,
                                       col_bigspace = 1,
                                       expand_xmin = 0,
                                       expand_xmax = 2,
                                       expand_ymin = 0,
                                       expand_ymax = 0,
                                       col_annot_angle = 30,
                                       col_annot_offset = 2
                                     )
)
funky


png(file = paste0(Results, "Fig_funky_rank_", st_num, "st.png"), 
    res = 450, width = funky$width + 1, height = funky$height, units = "in")
print(funky)
dev.off()
