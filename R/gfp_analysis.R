
# Integrative analysis
DefaultAssay(object = gfp_combined) <- "integrated"

gfp_combined <- subset(gfp_combined, cells = setdiff(names(gfp_combined@active.ident), c(doublets[[1]],doublets[[3]])))

gfp_combined <- ScaleData(object = gfp_combined,
                          vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"),
                          verbose = T)

gfp_combined <- analyze_merged(gfp_combined, group.levels = stages,
                               verbose = T, npcs = 50, dims = 1:20, nnei = 70, k.param = 70, min.dist = 0.8, spread = 1.5, resolution = 0.2)

# Visualization

plot_merge(gfp_combined)

plot_cluster(gfp_combined, label = F)

plot_split(gfp_combined, colors = get_colors(seq(levels(gfp_combined$seurat_clusters)), pal = "Set3"))

plot_features(gfp_combined, features = c("Krt14","Csn3","Lyz2","Ptn","S100a9","Col1a1","Cd3d","Il17ra"), ncol = 4)

# Identify cell markers

gfp_markers <- FindAllMarkers(gfp_combined, only.pos = T, logfc.threshold = 0.1)

# Heatmap

plot_heatmap(gfp_combined, gfp_markers, 6, cluster_pal = c("Set3"))

# Relabeling

gfp_labels <- c("Basal 1", 
                "Luminal 1", 
                "Macrophage", 
                "Luminal 2",
                "Neutrophil",
                "NS Immune",
                "Basal 2", 
                "Fibroblast")

gfp_levels <- gfp_labels[c(1,7,2,4,3,5,6,8)]

gfp_combined <- rename_cluster(gfp_combined, gfp_labels)

# New visualization

gfp_index <- c(5:8,4,2,10,12)
plot_cluster(gfp_combined, label = F, levels = gfp_levels, self_set_color = T, self_colors = get_colors(gfp_index))

# Statistics

plot_stat(gfp_combined, "group_count", group_levels = stages, cluster_levels = gfp_levels, plot_ratio = 3)

plot_stat(gfp_combined, "cluster_count", group_levels = stages, cluster_levels = gfp_levels,
          self_set_color = T, self_colors = get_colors(gfp_index))

plot_stat(gfp_combined, "prop_fill", group_levels = stages, cluster_levels = gfp_levels, plot_ratio = 3,
          self_set_color = T, self_colors = get_colors(gfp_index))

plot_stat(gfp_combined, "prop_diverge", group_levels = stages, cluster_levels = gfp_levels, plot_ratio = 0.8)

plot_stat(gfp_combined, "prop_multi", group_levels = stages, cluster_levels = gfp_levels, plot_ratio = 3)


# DE analysis

gfp_diff <- find_diff_genes(dataset = gfp_combined, clusters = gfp_levels, groups = stages, logfc = 0)

gfp_hallmark_gsea <- test_GSEA(gfp_diff, gfp_levels, pathways.hallmark)

plot_GSEA(gfp_hallmark_gsea, p_cutoff = 0.1, levels = gfp_levels)



# Test

hist(gfp_combined$nCount_RNA, breaks = 100)







