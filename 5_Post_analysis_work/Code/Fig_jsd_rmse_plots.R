################################################################################
# Generates the baseline JSD/RMSE & delta JSD/RMSE plots from the results 
# calculated earlier in delta_JSD_plots.R & delta_RMSE_plots.R files
# 
# JSD.mean.list.4ST.RDS / RMSE.mean.list.4ST.RDS
# Refers to 4 list, one for each ST dataset
# Each list has 7 sub-lists, one for each scenario (baseline + 6 removal scenario)
# Each sub-list is a matrix with columns as methods and rows as spots
# [i,j] is the mean JSD value of spot i for method j over multiple single cell
# reference dataset
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################


#### Initialize environment ####
source("Init_env.R")

# read the mean JSD/RMSE results for reproducing plots in the manuscript
jsd.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                            "JSD.mean.list.4ST.RDS"))
rmse.lists <- readRDS(paste0("../../4_Analysis_results/Results/",
                             "RMSE.mean.list.4ST.RDS"))

# apply((jsd.lists[[1]][[1]]), 2, median)
# apply((rmse.lists[[1]][[1]]), 2, median)

# change the value of 'st' variable based on which ST datasets to choose
# st=1 for 4-8|10-15; st=2 for 1-5|10-15, st=4 for 1-5|3-7 ST dataset config

for (st in 1:4) {
  if (st != 3) {
    
    jsdPlot1 <- reshape2::melt(jsd.lists[[st]])
    
    jsd.base <- ggplot(jsdPlot1, aes(x = variable, y = value)) +
      geom_violin(trim = T) +
      stat_summary(fun = median, geom = "point", shape = 4, color = "red", size = 3) +
      guides(fill = guide_legend("", nrow = 1)) +
      ylab("JSD") +
      scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
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
    
    jsd.other <- lapply(1:1, function(m) {
      mat.plot <- list()
      # first row for baseline
      for (mp in 1:6) { 
        mat.plot[[mp]] <- jsd.lists[[m]][[mp+1]] - jsd.lists[[m]][[1]]
        mat.plot[[mp]] <- mat.plot[[mp]] %>% dplyr::mutate(scene = rep(mp, ))
      }
      
      # mean.JSD.plot <- reshape2::melt(mat.plot[[m]])
      mean.JSD.plot <- reshape2::melt(rbind(mat.plot[[1]], mat.plot[[2]], mat.plot[[3]], mat.plot[[4]], 
                                            mat.plot[[5]], mat.plot[[6]]), id = "scene")
      mean.JSD.plot[is.na(mean.JSD.plot)] <- -0.5
      
      ggplot(mean.JSD.plot, aes(x = variable, y = value)) +
        geom_boxplot(aes(fill = as.factor(scene)), outlier.shape = NA, na.rm = F) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = .2) +
        ylab(expression(Delta*"JSD")) +
        scale_y_continuous(breaks = c(-1.00, -0.50, 0.00, 0.50, 1.00)) +
        theme_classic() +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12, color = "#05445E", hjust = 0),
          axis.text.x = element_text(size = 13, color = "dodgerblue4", vjust = 1),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, angle = 90, vjust = 0.5, color = "dodgerblue4"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95", size = 0.1)
        ) + NoLegend() +
        ggtitle("Removal scenario")
    })
    
    jsd.leg <- get_legend(jsd.other[[1]]
                          + theme(legend.position = "bottom",
                                  legend.title = element_text(size = 14, color = "black"),
                                  legend.text = element_text(size = 13),
                                  legend.key.size = unit(.7, "cm"))
                          + guides(fill = guide_legend(nrow = 1))
                          + scale_fill_discrete(name = "# missing cell types",
                                                labels = c(
                                                  "1",
                                                  "2",
                                                  "3",
                                                  "5",
                                                  "10",
                                                  "11")))
    
    
    png(file = paste0(Results, "Fig_jsd_st_", st, ".png"),
        res = 450, width = 9.6, height = 6, units = "in")
    print(plot_grid(jsd.base, jsd.other[[1]], jsd.leg, ncol = 1, labels = c("e", ""), rel_heights = c(.7, 1.4, .15)))
    dev.off()
    
    # tiff(file = paste0(Results, "Fig_jsd_st_",st,".tif"), res = 600, width = 9.6, height = 6, units = "in")
    # plot_grid(jsd.base, jsd.other[[1]], jsd.leg, ncol = 1, labels = c("e", ""), rel_heights = c(.7, 1.4, .15))
    # dev.off()
    
    
    
    rmsePlot1 <- reshape2::melt(rmse.lists[[st]])
    
    rmse.base <- ggplot(rmsePlot1, aes(x = variable, y = value)) +
      geom_violin(trim = T) +
      stat_summary(fun.y = median, geom = "point", shape = 4, color = "red", size = 3) +
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
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13, color = "black", hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, angle = 90, vjust = 0.5, color = "dodgerblue4"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey95", size = 0.1)
      ) + NoLegend() +
      ggtitle("Baseline scenario")
    
    rmse.other <- lapply(1:1, function(m) {
      mat.plot <- list()
      for (mp in 1:6) {
        mat.plot[[mp]] <- rmse.lists[[m]][[mp+1]] - rmse.lists[[m]][[1]]
        mat.plot[[mp]] <- mat.plot[[mp]] %>% dplyr::mutate(scene = rep(mp, ))
      }
      
      # mean.JSD.plot <- reshape2::melt(mat.plot[[m]])
      mean.RMSE.plot <- reshape2::melt(rbind(mat.plot[[1]], mat.plot[[2]],
                                             mat.plot[[3]], mat.plot[[4]],
                                             mat.plot[[5]], mat.plot[[6]]),
                                       id = "scene")
      
      mean.RMSE.plot[is.na(mean.RMSE.plot)] <- -0.25
      
      ggplot(mean.RMSE.plot, aes(x = variable, y = value)) +
        geom_boxplot(aes(fill = as.factor(scene)), outlier.shape = NA, na.rm = F) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = .3) + 
        ylab(expression(Delta*"RMSE")) +
        theme_classic() +
        theme(
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12, color = "#05445E", hjust = 0),
          axis.text.x = element_text(size = 13, color = "dodgerblue4", vjust = 1),
          axis.text.y = element_text(size = 13, color = "black", hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14, angle = 90, vjust = 0.5, color = "dodgerblue4"),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "grey95", size = 0.1)
        ) + NoLegend() +
        ggtitle("Removal scenario")
    })
    
    
    rmse.leg <- get_legend(rmse.other[[1]]
                           + theme(legend.position = "bottom",
                                   legend.title = element_text(size = 14, color = "black"),
                                   legend.text = element_text(size = 13),
                                   legend.key.size = unit(.7, "cm"))
                           + guides(fill = guide_legend(nrow = 1))
                           + scale_fill_discrete(name = "# missing cell types",
                                                 labels = c(
                                                   "1",
                                                   "2",
                                                   "3",
                                                   "5",
                                                   "10",
                                                   "11")))
    
    
    png(file = paste0(Results, "Fig_rmse_st_",st,".png"), res = 450, width = 9.6, height = 6, units = "in")
    print(plot_grid(rmse.base, rmse.other[[1]], rmse.leg, labels = c("f", ""), ncol = 1, rel_heights = c(.7, 1.4, .15)))
    dev.off()
    
    # tiff(file = paste0(Results, "Fig_rmse_st_",st,".tif"), res = 600, width = 9.6, height = 6, units = "in")
    # plot_grid(rmse.base, rmse.other[[1]], rmse.leg, labels = c("f", ""), ncol = 1, rel_heights = c(.7, 1.4, .15))
    # dev.off()
    
  }
}




