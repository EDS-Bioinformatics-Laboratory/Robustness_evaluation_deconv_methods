################################################################################
# Part 1
# The script generates a bar plot with number of cells per celltype/cluster
# which could be used to cross-verify the number of cells per cluster in original
# analysis and reproduction efforts for the integrated single-cell reference data
#
# Part 2
# Create heatmap plots for mean JSD/RMSE values for all scenarios in each ST dataset
# 
# Part 3
# Create heatmap plots for median JSD/RMSE values for all scenarios in each ST dataset
# 
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################

#### Initialize environment ####
source("Init_env.R")

args <- commandArgs(trailingOnly = TRUE)
st_num <- args[1]

#### Part 1 ####
sc.ref.int.data <- readRDS(paste0("../../1_Generate_sc_ref_data/Results/sc.combined.nor.rds"))


df <- as.data.frame(table(sc.ref.int.data$blue.main))
colnames(df) <- c("Cell_types", "Freq")

png(paste0(Results, "Single_cell_integrated_data.png"), height = 5.5,
    width = 13.5, units = "in", res = 450)

g <- ggplot(data = df, aes(y = Cell_types, x = Freq, fill = Cell_types)) +
  geom_bar(stat = "identity") +
  geom_text(
    aes(label = Freq),
    position = position_dodge(width = 1),
    vjust = 0.5,
    hjust = 0
  ) +
  theme(
    plot.title = element_text(size = 14, vjust = 0.5, color = "#05445E", face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, vjust = 0.5, color = "dodgerblue4"),
    axis.title.y = element_text(size = 14, vjust = 0.5, color = "dodgerblue4"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey95", linewidth = 0.1)
  ) + 
  NoLegend() + 
  ggtitle("Cell type distribution") + ylab("Cell types") + xlab("Number of cells")

print(g)

dev.off()



#### Part 2 ####

for (s in st_num:st_num) {
  
  jsd.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                              "JSD.means.", s, "_ST.rds"))
  rmse.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                               "RMSE.means.", s, "_ST.rds"))
  
  
  df <- do.call(rbind, lapply(lapply(jsd.lists, colMeans), function(x) as.data.frame(t(x))))
  df_long <- df %>%
    mutate(Row = row_number()) %>%
    pivot_longer(-Row, names_to = "Column", values_to = "Value")
  
  df2 <- do.call(rbind, lapply(lapply(rmse.lists, colMeans), function(x) as.data.frame(t(x))))
  df_long2 <- df2 %>%
    mutate(Row = row_number()) %>%
    pivot_longer(-Row, names_to = "Column", values_to = "Value")
  
  
  png(paste0(Results, "Mean_JSD_RMSE_values_", s, "st.png"), height = 5.5,
      width = 13.5, units = "in", res = 450)
  
  g1 <- ggplot(df_long, aes(y = Column, x = factor(Row), fill = Value)) +
    geom_tile(aes(fill = Value), color = "white", lwd = 1.5) +
    geom_text(aes(label = round(Value, 4)), color = "black", size = 4) +
    scale_fill_gradient2(high = "dodgerblue4", na.value = "grey40") +
    scale_x_discrete(labels = c("0", "1", "2", "3", "5", "10", "11")) +
    xlab("# cell types removed from reference data") +
    ylab("Deconvolution methods") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, color = "brown4", hjust = 0, face="bold"),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 12, face="bold"),
      axis.title.y = element_text(size = 12, face="bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
    ) +
    NoLegend() +
    ggtitle("Mean JSD values across removal scenarios")
  
  g2 <- ggplot(df_long2, aes(y = Column, x = factor(Row), fill = Value)) +
    geom_tile(aes(fill = Value), color = "white", lwd = 1.5) +
    geom_text(aes(label = round(Value, 4)), color = "black", size = 4) +
    scale_fill_gradient2(high = "dodgerblue4", na.value = "grey40") +
    scale_x_discrete(labels = c("0", "1", "2", "3", "5", "10", "11")) +
    xlab("# cell types removed from reference data") +
    ylab("Deconvolution methods") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, color = "brown4", hjust = 0, face="bold"),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 12, face="bold"),
      axis.title.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
    ) +
    NoLegend() +
    ggtitle("Mean RMSE values across removal scenarios")
  
  print(g1+g2)
  
  dev.off()
}


#### Part 3 ####

for (s in st_num:st_num) {
  
  jsd.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                              "JSD.means.", s, "_ST.rds"))
  rmse.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                               "RMSE.means.", s, "_ST.rds"))
  

  df <- do.call(rbind, lapply(jsd.lists, function(mat) {
    # Calculate column medians using apply
    medians <- apply(mat, 2, median, na.rm = TRUE)
    # Convert to a dataframe and transpose
    as.data.frame(t(medians))
  }))
  
  df_long <- df %>%
    mutate(Row = row_number()) %>%
    pivot_longer(-Row, names_to = "Column", values_to = "Value")
  
  df2 <- do.call(rbind, lapply(rmse.lists, function(mat) {
    # Calculate column medians using apply
    medians <- apply(mat, 2, median, na.rm = TRUE)
    # Convert to a dataframe and transpose
    as.data.frame(t(medians))
  }))
  df_long2 <- df2 %>%
    mutate(Row = row_number()) %>%
    pivot_longer(-Row, names_to = "Column", values_to = "Value")
  
  
  png(paste0(Results, "Median_JSD_RMSE_values_", s, "st.png"), height = 5.5,
      width = 13.5, units = "in", res = 450)
  
  g1 <- ggplot(df_long, aes(y = Column, x = factor(Row), fill = Value)) +
    geom_tile(aes(fill = Value), color = "white", lwd = 1.5) +
    geom_text(aes(label = round(Value, 4)), color = "black", size = 4) +
    scale_fill_gradient2(high = "dodgerblue4", na.value = "grey40") +
    scale_x_discrete(labels = c("0", "1", "2", "3", "5", "10", "11")) +
    xlab("# cell types removed from reference data") +
    ylab("Deconvolution methods") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, color = "brown4", hjust = 0, face="bold"),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 12, face="bold"),
      axis.title.y = element_text(size = 12, face="bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
    ) +
    NoLegend() +
    ggtitle("Median JSD values across removal scenarios")
  
  g2 <- ggplot(df_long2, aes(y = Column, x = factor(Row), fill = Value)) +
    geom_tile(aes(fill = Value), color = "white", lwd = 1.5) +
    geom_text(aes(label = round(Value, 4)), color = "black", size = 4) +
    scale_fill_gradient2(high = "dodgerblue4", na.value = "grey40") +
    scale_x_discrete(labels = c("0", "1", "2", "3", "5", "10", "11")) +
    xlab("# cell types removed from reference data") +
    ylab("Deconvolution methods") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, color = "brown4", hjust = 0, face="bold"),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 12, face="bold"),
      axis.title.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
    ) +
    NoLegend() +
    ggtitle("Median RMSE values across removal scenarios")
  
  print(g1+g2)
  
  dev.off()
}

