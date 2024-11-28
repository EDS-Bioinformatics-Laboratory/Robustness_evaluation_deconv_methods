################################################################################
# 
# This script plots the celltype assignment plots 
# 
# The celltype reassignment calculations for ST 1, 2, 3 and removal scenarios
# rm1, rm2, rm3 are all independent R objects
# e.g., for ST 1 and removal scenario rm1 -> CT_assign_rm1_1_ST.rds
# 
# Every object comprises a list of 8 matrices, one each for a deconvolution method
# Each column in the matrix refers to celltypes and rows refers to one or more
# celltypes removed from the reference single cell reference data
#
# Implemented by Utkarsh Mahamune
# Bioinformatics Laboratory | Amsterdam UMC (Location AMC)
################################################################################

#### Initialize environment ####
source("Init_env.R")

args <- commandArgs(trailingOnly = TRUE)
st_num <- args[1]

methods_ <- c("cell2location",
              "RCTD",
              "CARD",
              "SCDC",
              "MuSiC",
              "Stereoscope",
              "Seurat",
              "SPOTlight")


sc.data.ref <- readRDS(paste0("../../1_Generate_sc_ref_data/Results/",
                              "sc.ref.data.rds"))
celltypes <- sort(unique(sc.data.ref$blue.main))
rm(sc.data.ref)


#### Celltype reassignment plots for rm1 scenario ####
dflist <- readRDS(paste0("../../4_Analysis_results/Results/", "CT_assign_rm1.", st_num,"_ST.rds"))

png(file = paste0(Results, "Fig_Celltype_assignment_rm1_", st_num, "st.png"),
    res = 450, width = 10.8, height = 4.5, units = "in")

ct.assignment.plots <- lapply(1:dim(dflist)[3], function(m) {
  dp <- dflist[,,m]
  colnames(dp) <- celltypes
  rownames(dp) <- celltypes
  
  dp <- reshape2::melt(dp)
  
  if (m %in% c(1, 5)) {
    ggplot(dp, aes(y = forcats::fct_rev(Var1), x = Var2)) +
      geom_tile(aes(fill = value), color = "white", lwd = .5) +
      scale_fill_gradient2(low = "brown3", high = "deepskyblue4", mid = "white",
                           midpoint = 0, na.value = "grey40") +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        legend.box.margin = margin(t = -10, r = 20, b = -20, l = -10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"),
        plot.title = element_text(size = 8, color = "dodgerblue4", hjust = 0),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.2, hjust = .95),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
      ) +
      ggtitle(methods_[m]) +
      scale_x_discrete(labels = c("Adipocytes" = "Adipo",
                                  "B_cells" = "B",
                                  "CD4_T_cells" = "CD4 T",
                                  "CD8_T_cells" = "CD8 T",
                                  "Endothelial_cells" = "Endo",
                                  "Fibroblasts" = "Fibro",
                                  "HSC" = "HSC",
                                  "Macrophages" = "Macro",
                                  "Monocytes" = "Mono",
                                  "Myocytes" = "Myo",
                                  "Neutrophils" = "Neutro",
                                  "NK_cells" = "NK",
                                  "Skeletal_muscle" = "SkeMus")) +
      scale_y_discrete(labels = c("Adipocytes" = "Adipo",
                                  "B_cells" = "B",
                                  "CD4_T_cells" = "CD4 T",
                                  "CD8_T_cells" = "CD8 T",
                                  "Endothelial_cells" = "Endo",
                                  "Fibroblasts" = "Fibro",
                                  "HSC" = "HSC",
                                  "Macrophages" = "Macro",
                                  "Monocytes" = "Mono",
                                  "Myocytes" = "Myo",
                                  "Neutrophils" = "Neutro",
                                  "NK_cells" = "NK",
                                  "Skeletal_muscle" = "SkeMus"))
  } else {
    ggplot(dp, aes(y = forcats::fct_rev(Var1), x = Var2)) +
      geom_tile(aes(fill = value), color = "white", lwd = .5) +
      scale_fill_gradient2(low = "brown3", high = "deepskyblue4", mid = "white",
                           midpoint = 0, na.value = "grey40") +
      xlab(methods_[m]) +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        legend.box.margin = margin(t = -10, r = 20, b = -20, l = -10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"),
        plot.title = element_text(size = 8, color = "dodgerblue4", hjust = 0),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.2, hjust = .95),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
      ) +
      ggtitle(methods_[m]) +
      scale_x_discrete(labels = c("Adipocytes" = "Adipo",
                                  "B_cells" = "B",
                                  "CD4_T_cells" = "CD4 T",
                                  "CD8_T_cells" = "CD8 T",
                                  "Endothelial_cells" = "Endo",
                                  "Fibroblasts" = "Fibro",
                                  "HSC" = "HSC",
                                  "Macrophages" = "Macro",
                                  "Monocytes" = "Mono",
                                  "Myocytes" = "Myo",
                                  "Neutrophils" = "Neutro",
                                  "NK_cells" = "NK",
                                  "Skeletal_muscle" = "SkeMus"))
  }
})

# using cowplot (plot_grid function)
l1 <- cowplot::plot_grid(ct.assignment.plots[[1]], NULL, ct.assignment.plots[[2]], NULL,
                         ct.assignment.plots[[3]], NULL, ct.assignment.plots[[4]], NULL,
                         nrow = 1, align = "hv", rel_widths = c(1, -.12, 1, -.12, 1, -.12, 1, -.12))

l2 <- cowplot::plot_grid(ct.assignment.plots[[5]], NULL, ct.assignment.plots[[6]], NULL,
                         ct.assignment.plots[[7]], NULL, ct.assignment.plots[[8]], NULL,
                         nrow = 1, align = "hv", rel_widths = c(1, -.12, 1, -.12, 1, -.12, 1, -.12))
prow_ct <- cowplot::plot_grid(l1, NULL, l2,
                              align = 'hv', ncol = 1, rel_heights = c(1, 0, 1))

y.grob <- grid::textGrob("Removed cell type (in grey)", 
                         gp=gpar(fontface = "bold", col = "dodgerblue4", fontsize = 10), rot = 90)
print(gridExtra::grid.arrange(arrangeGrob(prow_ct, left = y.grob)))

dev.off()




#### Celltype reassignment plots for rm2 scenario ####
dflist <- readRDS(paste0("../../4_Analysis_results/Results/", "CT_assign_rm2.", st_num,"_ST.rds"))

png(file = paste0(Results, "Fig_Celltype_assignment_rm2_", st_num, "st.png"),
    res = 450, width = 10.8, height = 3.2, units = "in")

ct.assignment.plots2 <- lapply(1:dim(dflist)[3], function(m) {
  dp <- dflist[,,m]
  colnames(dp) <- celltypes
  rownames(dp) <- c("CD4 T & CD8 T", "Myo & SkeMus", "HSC & NK", "Macro & Mono",
                    "Adipo & Fibro")
  
  dp <- reshape2::melt(dp)
  
  if (m %in% c(1, 5)) {
    ggplot(dp, aes(y = forcats::fct_rev(Var1), x = Var2)) +
      geom_tile(aes(fill = value), color = "white", lwd = .5) +
      scale_fill_gradient2(low = "brown3", high = "deepskyblue4", mid = "white",
                           midpoint = 0, na.value = "grey40") +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        legend.box.margin = margin(t = -10, r = 20, b = -20, l = -10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"),
        plot.title = element_text(size = 8, color = "dodgerblue4", hjust = 0),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.2, hjust = .95),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
      ) +
      ggtitle(methods_[m]) +
      scale_x_discrete(labels = c("Adipocytes" = "Adipo",
                                  "B_cells" = "B",
                                  "CD4_T_cells" = "CD4 T",
                                  "CD8_T_cells" = "CD8 T",
                                  "Endothelial_cells" = "Endo",
                                  "Fibroblasts" = "Fibro",
                                  "HSC" = "HSC",
                                  "Macrophages" = "Macro",
                                  "Monocytes" = "Mono",
                                  "Myocytes" = "Myo",
                                  "Neutrophils" = "Neutro",
                                  "NK_cells" = "NK",
                                  "Skeletal_muscle" = "SkeMus"))
  } else {
    ggplot(dp, aes(y = forcats::fct_rev(Var1), x = Var2)) +
      geom_tile(aes(fill = value), color = "white", lwd = .5) +
      scale_fill_gradient2(low = "brown3", high = "deepskyblue4", mid = "white",
                           midpoint = 0, na.value = "grey40") +
      xlab(methods_[m]) +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        legend.box.margin = margin(t = -10, r = 20, b = -20, l = -10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"),
        plot.title = element_text(size = 8, color = "dodgerblue4", hjust = 0),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.2, hjust = .95),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
      ) +
      ggtitle(methods_[m]) +
      scale_x_discrete(labels = c("Adipocytes" = "Adipo",
                                  "B_cells" = "B",
                                  "CD4_T_cells" = "CD4 T",
                                  "CD8_T_cells" = "CD8 T",
                                  "Endothelial_cells" = "Endo",
                                  "Fibroblasts" = "Fibro",
                                  "HSC" = "HSC",
                                  "Macrophages" = "Macro",
                                  "Monocytes" = "Mono",
                                  "Myocytes" = "Myo",
                                  "Neutrophils" = "Neutro",
                                  "NK_cells" = "NK",
                                  "Skeletal_muscle" = "SkeMus"))
  }
  
})


# using cowplot (plot_grid function)
l12 <- cowplot::plot_grid(ct.assignment.plots2[[1]], NULL, ct.assignment.plots2[[2]], NULL,
                          ct.assignment.plots2[[3]], NULL, ct.assignment.plots2[[4]], NULL,
                          nrow = 1, align = "hv", rel_widths = c(1, -.12, 1, -.12, 1, -.12, 1, -.1))

l22 <- cowplot::plot_grid(ct.assignment.plots2[[5]], NULL, ct.assignment.plots2[[6]], NULL,
                          ct.assignment.plots2[[7]], NULL, ct.assignment.plots2[[8]], NULL,
                          nrow = 1, align = "hv", rel_widths = c(1, -.12, 1, -.12, 1, -.12, 1, -.1))
prow_ct2 <- cowplot::plot_grid(l12, NULL, l22,
                               align = 'hv', ncol = 1, rel_heights = c(1, 0, 1))

y.grob2 <- grid::textGrob("Removed cell types (in grey)",
                          gp=gpar(fontface = "bold", col = "dodgerblue4", fontsize = 10), rot = 90)
print(gridExtra::grid.arrange(arrangeGrob(prow_ct2, left = y.grob2)))


dev.off()


#### Celltype reassignment plots for rm2 scenario ####
dflist <- readRDS(paste0("../../4_Analysis_results/Results/", "CT_assign_rm3.", st_num,"_ST.rds"))

png(file = paste0(Results, "Fig_Celltype_assignment_rm3_", st_num, "st.png"),
    res = 450, width = 10.8, height = 3.8, units = "in")

ct.assignment.plots3 <- lapply(1:dim(dflist)[3], function(m) {
  dp <- dflist[,,m]
  colnames(dp) <- celltypes
  rownames(dp) <- c("NK, CD4 T\n& CD8 T", "Myo, SkeMus\n& Fibro", "B, CD4 T\n& CD8 T",
                    "NK, Macro \n& Mono", "Adipo, Fibro\n Endo")
  
  dp <- reshape2::melt(dp)
  
  if (m %in% c(1, 5)) {
    ggplot(dp, aes(y = forcats::fct_rev(Var1), x = Var2)) +
      geom_tile(aes(fill = value), color = "white", lwd = .5) +
      scale_fill_gradient2(low = "brown3", high = "deepskyblue4", mid = "white",
                           midpoint = 0, na.value = "grey40") +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        legend.box.margin = margin(t = -10, r = 20, b = -20, l = -10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"),
        plot.title = element_text(size = 8, color = "dodgerblue4", hjust = 0),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.2, hjust = .95),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey70", size = 0.1)
      ) +
      ggtitle(methods_[m]) +
      scale_x_discrete(labels = c("Adipocytes" = "Adipo",
                                  "B_cells" = "B",
                                  "CD4_T_cells" = "CD4 T",
                                  "CD8_T_cells" = "CD8 T",
                                  "Endothelial_cells" = "Endo",
                                  "Fibroblasts" = "Fibro",
                                  "HSC" = "HSC",
                                  "Macrophages" = "Macro",
                                  "Monocytes" = "Mono",
                                  "Myocytes" = "Myo",
                                  "Neutrophils" = "Neutro",
                                  "NK_cells" = "NK",
                                  "Skeletal_muscle" = "SkeMus"))
  } else {
    ggplot(dp, aes(y = forcats::fct_rev(Var1), x = Var2)) +
      geom_tile(aes(fill = value), color = "white", lwd = .5) +
      scale_fill_gradient2(low = "brown3", high = "deepskyblue4", mid = "white",
                           midpoint = 0, na.value = "grey40") +
      xlab(methods_[m]) +
      theme_classic() +
      theme(
        legend.title = element_blank(),
        legend.box.margin = margin(t = -10, r = 20, b = -20, l = -10),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.4, "cm"),
        plot.title = element_text(size = 8, color = "dodgerblue4", hjust = 0),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.2, hjust = .95),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey70", linewidth = 0.1)
      ) +
      ggtitle(methods_[m]) +
      scale_x_discrete(labels = c("Adipocytes" = "Adipo",
                                  "B_cells" = "B",
                                  "CD4_T_cells" = "CD4 T",
                                  "CD8_T_cells" = "CD8 T",
                                  "Endothelial_cells" = "Endo",
                                  "Fibroblasts" = "Fibro",
                                  "HSC" = "HSC",
                                  "Macrophages" = "Macro",
                                  "Monocytes" = "Mono",
                                  "Myocytes" = "Myo",
                                  "Neutrophils" = "Neutro",
                                  "NK_cells" = "NK",
                                  "Skeletal_muscle" = "SkeMus"))
  }
  
})


# using cowplot (plot_grid function)
l13 <- cowplot::plot_grid(ct.assignment.plots3[[1]], NULL, ct.assignment.plots3[[2]], NULL,
                          ct.assignment.plots3[[3]], NULL, ct.assignment.plots3[[4]], NULL,
                          nrow = 1, align = "hv", rel_widths = c(1, -.12, 1, -.12, 1, -.12, 1, -.1))

l23 <- cowplot::plot_grid(ct.assignment.plots3[[5]], NULL, ct.assignment.plots3[[6]], NULL,
                          ct.assignment.plots3[[7]], NULL, ct.assignment.plots3[[8]], NULL,
                          nrow = 1, align = "hv", rel_widths = c(1, -.12, 1, -.12, 1, -.12, 1, -.1))
prow_ct3 <- cowplot::plot_grid(l13, NULL, l23,
                               align = 'hv', ncol = 1, rel_heights = c(1, 0, 1))

y.grob2 <- grid::textGrob("Removed cell types (in grey)",
                          gp=gpar(fontface = "bold", col = "dodgerblue4", fontsize = 10), rot = 90)
print(gridExtra::grid.arrange(arrangeGrob(prow_ct3, left = y.grob2)))

dev.off()