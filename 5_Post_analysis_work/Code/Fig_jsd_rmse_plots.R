################################################################################
# Generates the baseline JSD/RMSE & delta JSD/RMSE plots from the results 
# calculated earlier in Get_JSD_results.R and Get_RMSE_results.R files under
# '4_Analysis_results' section
# 
# JSD.mean.1_ST.rds, JSD.mean.2_ST.rds, JSD.mean.3_ST.rds or
# RMSE.mean.1_ST.rds, RMSE.mean.2_ST.rds, RMSE.mean.3_ST.rds
# 
# Each one of the above object refers to one ST dataset
# Each object has 7 lists, one for each scenario (baseline + 6 removal scenario)
# Each sub-list is a matrix with columns as methods and rows as spots
# [i,j] is the mean JSD/RMSE value of spot i for method j over multiple single cell
# reference datasets
# 
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################


#### Initialize environment ####
source("Init_env.R")

args <- commandArgs(trailingOnly = TRUE)
st_num <- args[1]


# read the mean JSD/RMSE results for reproducing plots in the manuscript
jsd.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                            "JSD.means.", st_num, "_ST.rds"))
rmse.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                             "RMSE.means.", st_num, "_ST.rds"))

# The methods are ordered based on the baseline JSD result (even for rmse plots)
jsd.list <- jsd.lists[[1]]
new_col_order <- names(sort(apply(jsd.list, 2, median)))

# reordered JSD and RMSE results
jsd.lists <- lapply(jsd.lists, function(df) df[, new_col_order])
rmse.lists <- lapply(rmse.lists, function(df) df[, new_col_order])


# baseline scenarios
jsd.list <- jsd.lists[[1]]
jsdPlot1 <- reshape2::melt(jsd.list)

rmse.list <- rmse.lists[[1]]
rmsePlot1 <- reshape2::melt(rmse.list)


jsd.base <- ggplot(jsdPlot1, aes(x = variable, y = value)) +
  geom_violin(trim = T) +
  stat_summary(fun = median, geom = "point", shape = 4, color = "red", size = 3) +
  guides(fill = guide_legend("", nrow = 1)) +
  ylab("JSD") +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.key.size = unit(.8, "cm"),
    legend.box.margin = margin(t = -10, r = -5, b = -10, l = -20),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(size = 12, color = "#05445E", hjust = 0),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, angle = 90, vjust = 0.5, color = "dodgerblue4"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey95", linewidth = 0.1)
  ) + NoLegend() +
  ggtitle("Baseline scenario")

rmse.base <- ggplot(rmsePlot1, aes(x = variable, y = value)) +
  geom_violin(trim = T) +
  stat_summary(fun = median, geom = "point", shape = 4, color = "red", size = 3) +
  guides(fill = guide_legend("", nrow = 1)) +
  ylab("RMSE") +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2, 0.3), limits = c(0, 0.35)) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.key.size = unit(.8, "cm"),
    legend.box.margin = margin(t = -10, r = -5, b = -10, l = -20),
    panel.background = element_rect(fill = "white"),
    plot.title = element_text(size = 12, color = "#05445E", hjust = 0),
    axis.text.x = element_text(size = 12, color = "dodgerblue4", vjust = 1),
    axis.text.y = element_text(size = 13, color = "black", hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, angle = 90, vjust = 0.5, color = "dodgerblue4"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey95", linewidth = 0.1)
  ) + NoLegend() +
  ggtitle("Baseline scenario")

png(file = paste0(Results, "Fig_baseline_st_", st_num, ".png"),
    res = 450, width = 9.6, height = 6, units = "in")
print(plot_grid(jsd.base, rmse.base,  ncol = 1, rel_heights = c(.7, .8, .15)))
dev.off()


jsd.list <- jsd.lists
# removal scenarios
jsd.other <- lapply(1:1, function(m) {
  mat.plot <- list()
  # first row for baseline
  for (mp in 1:(length(jsd.lists)-1)) { 
    mat.plot[[mp]] <- jsd.list[[mp+1]] - jsd.list[[1]]
    mat.plot[[mp]] <- mat.plot[[mp]] %>% data.frame() %>% dplyr::mutate(scene = rep(mp, ))
  }
  
  mean.JSD.plot <- reshape2::melt(do.call(rbind, mat.plot), id = "scene")
  
  mean.JSD.plot[is.na(mean.JSD.plot)] <- -0.5
  
  ggplot(mean.JSD.plot, aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = as.factor(scene)), outlier.shape = NA, na.rm = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = .2) +
    ylab(expression(Delta*"JSD")) +
    # scale_y_continuous(breaks = c(-1, -0.50, 0.00, 0.50, 1.00), limits = c(-.75, 1)) +
    theme_classic() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, color = "#05445E", hjust = 0),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 13, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14, angle = 90, vjust = 0.5, color = "dodgerblue4"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "grey95", linewidth = 0.1)
    ) + NoLegend() +
    ggtitle("Cell type mismatch scenario")
})

rmse.list <- rmse.lists
rmse.other <- lapply(1:1, function(m) {
  mat.plot <- list()
  for (mp in 1:(length(rmse.list)-1)) { 
    mat.plot[[mp]] <- rmse.list[[mp+1]] - rmse.list[[1]]
    mat.plot[[mp]] <- mat.plot[[mp]] %>% data.frame() %>% dplyr::mutate(scene = rep(mp, ))
  }
  
  mean.RMSE.plot <- reshape2::melt(do.call(rbind, mat.plot), id = "scene")
  
  mean.RMSE.plot[is.na(mean.RMSE.plot)] <- -0.25
  
  ggplot(mean.RMSE.plot, aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = as.factor(scene)), outlier.shape = NA, na.rm = F)
  
  ggplot(mean.RMSE.plot, aes(x = variable, y = value)) +
    geom_boxplot(aes(fill = as.factor(scene)), outlier.shape = NA, na.rm = F) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = .2) + 
    ylab(expression(Delta*"RMSE")) +
    # ylim (-.3, .35) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.5, "cm"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, color = "#05445E", hjust = 0),
      axis.text.x = element_text(size = 12, color = "dodgerblue4", vjust = 1),
      axis.text.y = element_text(size = 13, color = "black", hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 14, angle = 90, vjust = 0.5, color = "dodgerblue4"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(colour = "grey95", linewidth = 0.1)
    ) + 
    guides(fill = guide_legend(nrow = 1)) +
    scale_fill_discrete(name = "# missing cell types     ",
                        labels = c("1", "2", "3", "5", "10", "11")) + 
    ggtitle("Cell type mismatch scenario")
})


png(file = paste0(Results, "Fig_removal_scenario_st_", st_num, ".png"),
    res = 450, width = 9.6, height = 6, units = "in")
print(plot_grid(jsd.other[[1]], rmse.other[[1]], ncol = 1, rel_heights = c(.7, .9)))
dev.off()



# Calculate median values for the delta JSD per combination of 
# mismatch scenario and deconvolution method
jsd.other <- lapply(1:1, function(m) {
  mat.plot <- list()
  # first row for baseline
  for (mp in 1:(length(jsd.lists)-1)) { 
    mat.plot[[mp]] <- jsd.list[[mp+1]] - jsd.list[[1]]
    mat.plot[[mp]] <- mat.plot[[mp]] %>% data.frame() %>% dplyr::mutate(scene = rep(mp, ))
  }
  
  mean.JSD.plot <- reshape2::melt(do.call(rbind, mat.plot), id = "scene")
  
  mean.JSD.plot[is.na(mean.JSD.plot)] <- -0.5
  mean.JSD.plot
})
delta.jsd.median <- aggregate(value ~ scene + variable, data = jsd.other[[1]], median)
# Put missing values back
delta.jsd.median %>% 
  mutate(value=ifelse(variable == "Seurat" & scene %in% c(5,6), NA, value)) %>%
  mutate(value=ifelse(variable == "CARD" & scene ==6, NA, value)) ->
  delta.jsd.median

# Plot the delta JSD median values and add a regression line
  gg1 <- ggplot(delta.jsd.median, aes(scene, value, colour=variable)) +
    geom_point()  + 
  geom_smooth(aes(scene, value, colour=variable), method=lm, se=FALSE) +
  facet_wrap(~ variable, nrow = 1, scales="free_x") +
  labs(x = "", y = expression(Delta*"JSD (median)")) + 
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5, 6),  # Customize breaks
    labels = c("1", "2", "3", "5", "10", "11")  # Custom labels
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, vjust = 1)) +
  NoLegend()

# Calculate median values for the delta RMSE per combination of 
# mismatch scenario and deconvolution method
rmse.other <- lapply(1:1, function(m) {
  mat.plot <- list()
  for (mp in 1:(length(rmse.list)-1)) { 
    mat.plot[[mp]] <- rmse.list[[mp+1]] - rmse.list[[1]]
    mat.plot[[mp]] <- mat.plot[[mp]] %>% data.frame() %>% dplyr::mutate(scene = rep(mp, ))
  }
  
  mean.RMSE.plot <- reshape2::melt(do.call(rbind, mat.plot), id = "scene")
  
  mean.RMSE.plot[is.na(mean.RMSE.plot)] <- -0.25
  mean.RMSE.plot
})
delta.rmse.median <- aggregate(value ~ scene + variable, data = rmse.other[[1]], median)
# Put missing values back
delta.rmse.median %>% 
  mutate(value=ifelse(variable == "Seurat" & scene %in% c(5,6), NA, value)) %>%
  mutate(value=ifelse(variable == "CARD" & scene ==6, NA, value)) ->
  delta.rmse.median

# Plot the delta RMSE median values and add a regression line
gg2 <- ggplot(delta.rmse.median, aes(scene, value, colour=variable)) +
  geom_point()  + 
  geom_smooth(aes(scene, value, colour=variable), method=lm, se=FALSE) +
  facet_wrap(~ variable, nrow = 1, scales="free_x") +
  labs(x = "# missing cell types", y = expression(Delta*"RMSE (median)")) + 
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5, 6),  # Customize breaks
    labels = c("1", "2", "3", "5", "10", "11")  # Custom labels
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8, vjust = 1),
        axis.title.x = element_text(size = 10, vjust = 1)) +
  NoLegend()

png(file = paste0(Results, "delta_median_st_", st_num, ".png"),
    res = 450, width = 9.6, height = 6, units = "in")
print(gg1+gg2+patchwork::plot_layout(nrow=2))
dev.off()