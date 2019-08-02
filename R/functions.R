
library(tidyverse)

folder_names <- c("data", "refs", "R", "analysis", "figures", "man")
map(folder_names, dir.create) #purrr-fect way

# Load packages and refs for analysis

library(Seurat)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(colormap)
library(ggrepel)
library(reshape2)
library(formattable)
library(magrittr)
library(sctransform)
library(ggpubr)
library(scales)

library(clusterProfiler)
library(org.Mm.eg.db)
library(fgsea)
library(biomaRt)

library(Matrix)
library(fields)
library(KernSmooth)
library(modes)
library(ROCR)
library(DoubletFinder)

library(foreach)
library(doParallel)

mm_hs <- read_tsv("~/Documents/r_projects/HER2_scRNAseq/refs/mm_hs.txt", col_names = T)
pathways.hallmark <- gmtPathways("~/Documents/r_projects/HER2_scRNAseq/refs/h.all.v6.2.symbols.gmt")

# Helper functions

#version 1
vlookup <- function(list, table, col_1, col_2) {
        new_list <- table[[col_2]][match(list, table[[col_1]])]
        return(new_list)
}

load_10X_mito_cc <- function(dir = NULL, raw_data = NULL, gcol = 2, proj_name, min_cells = 5, org = "human", cc = T) {
        
        if (!is.null(dir))  raw_data <- Read10X(data.dir = dir, gene.column = gcol)
        
        data <- CreateSeuratObject(counts = raw_data, project = proj_name, min.cells = min_cells)
        
        mito_pattern <- ifelse(org == 'human', "^MT-", "^mt-")
        
        if (org == "mouse") {
                
                s.genes <- vlookup(cc.genes$s.genes, mm_hs, 2, 1)
                g2m.genes <- vlookup(cc.genes$g2m.genes, mm_hs, 2, 1)
                
                s.genes <- s.genes[!is.na(s.genes)]
                g2m.genes <- g2m.genes[!is.na(g2m.genes)]
                
        } else {
                
                s.genes <- cc.genes$s.genes
                g2m.genes <- cc.genes$g2m.genes
        }
        
        data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = mito_pattern)
        
        if (cc) data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
        
        data$group <- proj_name
        
        return(data)
}

plot_qc <- function(data_list, metrics) {
        
        qc <- list()
        
        for (i in seq(length(data_list))) {
                qc[[i]] <- tibble(value = data_list[[i]][[metrics]][[1]],
                                  sample = names(data_list)[i])
        }
        qc <- do.call(rbind, qc)
        ggplot(qc) + 
                geom_boxplot(aes(x = sample, y = value, fill = sample)) +
                scale_fill_brewer(palette = "Set3") +
                labs(y = metrics)
}

filter_norm_10X <- function(dataset, nfeature = 500, mito = 10, nfeatures = 2000) {
        
        expr1 <- FetchData(dataset, vars = "nFeature_RNA")
        expr2 <- FetchData(dataset, vars = "percent.mt")
        
        dataset <- dataset[, which(x = expr1 > nfeature & expr2 < mito)]
        dataset %<>% 
                NormalizeData(verbose = FALSE) %>%
                FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures)
        
}

filter_sctrans_10X <- function(dataset, nfeature = 500, mito = 10, 
                               vars = c("percent.mt","nCount_RNA","S.Score","G2M.Score")) {
        
        expr1 <- FetchData(dataset, vars = "nFeature_RNA")
        expr2 <- FetchData(dataset, vars = "percent.mt")
        
        dataset <- dataset[, which(x = expr1 > nfeature & expr2 < mito)]
        dataset %<>% 
                SCTransform(vars.to.regress = vars, verbose = T)
        
}

find_doublets <- function(dataset, dims = 1:20, ratio = 0.05, resolution = 0.4, txt) {
        
        dataset %<>%
                ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>%
                RunPCA(features = VariableFeatures(dataset)) %>%
                RunUMAP(dims = dims) %>%
                FindNeighbors(dims = dims) %>%
                FindClusters(resolution = resolution)
        
        ## pK Identification
        sweep.res <- paramSweep_v3(dataset, PCs = dims)
        sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        
        ## Homotypic Doublet Proportion Estimate
        homotypic.prop <- modelHomotypic(dataset$seurat_clusters)
        nExp_poi <- round(ratio*length(Idents(dataset)))
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        dataset <- doubletFinder_v3(dataset, PCs = dims, pN = 0.25, pK = 0.1, nExp = nExp_poi.adj, reuse.pANN = F)
        
        barcodes <- names(dataset@active.ident)[dataset[[paste("DF.classifications_0.25_0.1_", as.character(nExp_poi.adj), sep = "")]] == "Doublet"]
        
        paste(barcodes, txt, sep = "")
        
}

analyze_merged <- function(dataset, group.levels,
                           verbose = T, npcs = 50,
                           reduction = "umap", dims = 1:20, nnei = 30, min.dist = 0.3, spread = 1, n.epochs = 500, 
                           k.param = 20,
                           resolution = 0.8) {
        
        dataset$group <- factor(dataset$group, levels = group.levels)
        dataset %<>% 
                RunPCA(npcs = npcs, verbose = verbose) %>%
                RunUMAP(reduction = "pca", dims = dims, n.neighbors = nnei, min.dist = min.dist, spread = spread, n.epochs = n.epochs) %>%
                FindNeighbors(reduction = reduction, dims = 1:2, k.param = k.param) %>%
                FindClusters(resolution = resolution, algorithm = 3)
        
}

get_palette <- function(ncol, palette = c("Paired", "Set2")) {
        
        n1 <- brewer.pal.info[palette[1],][[1]]
        n2 <- brewer.pal.info[palette[2],][[1]]
        ful_pal <- c(brewer.pal(n = n1, name = palette[1]), brewer.pal(n = n2, name = palette[2]))
        pal <- ful_pal[1:ncol]
        return(pal)
        
}

get_colors <- function(v, pal = "Paired") {
        
        ncol <- brewer.pal.info[pal,][[1]]
        
        if (sum(!(v %in% 1:ncol)) > 0) {
                stop("Please input a valid numeric vector")
        }
        return(brewer.pal(n = ncol, name = pal)[v])
        
}

plot_merge <- function(dataset, reduction = "umap", group.by = "group", 
                       colors = c('#7bccc4','#f03b20'), legend.title = "Group", labels = levels(dataset$group)) {
        
        p <- DimPlot(object = dataset, reduction = reduction, group.by = group.by)
        
        p  + scale_color_manual(values = colors, name = legend.title, labels = labels) +
                theme(panel.border = element_rect(colour = "black", fill=NA, size=1, linetype = 1),
                      axis.line=element_blank(),
                      aspect.ratio = 1
                )
}

plot_split <- function(dataset, reduction = "umap", split.by = "group", 
                       colors = c('#7bccc4','#f03b20'), legend.title = "Cluster", labels = levels(dataset$seurat_clusters)) {
        
        p <- DimPlot(object = dataset, reduction = reduction, split.by = split.by)
        
        p  + scale_color_manual(values = colors, name = legend.title, labels = labels) +
                theme(panel.border = element_rect(colour = "black", fill = NA, size = 1, linetype = 1),
                      strip.text.x = element_text(face = "plain", vjust = 1),
                      axis.line=element_blank(),
                      aspect.ratio = 1
                )
}

plot_cluster <- function(dataset, reduction = "umap", label = T, levels = NULL,
                         self_set_color = F,
                         self_colors,
                         palette = c("Set2", "Paired")) {
        
        ncol <- length(levels(Idents(dataset)))
        colors <- if (self_set_color) self_colors else (get_palette(ncol, palette))
        
        tmp <- dataset
        
        if (is.null(levels) == F) {
                Idents(tmp) <- factor(Idents(tmp), levels = levels)
        }
        
        p <- DimPlot(object = tmp, reduction = reduction, label = label)
        p + scale_color_manual(values = colors, name = "Clusters") +
                theme(panel.border = element_rect(colour = "black", fill=NA, size=1, linetype = 1),
                      axis.line=element_blank(),
                      strip.text = element_blank(),
                      aspect.ratio = 1
                )
}

plot_features <- function(dataset, features, ncol) {
        
        DefaultAssay(dataset) <- "RNA"
        p_gene <- FeaturePlot(object = dataset, 
                              features = features, 
                              min.cutoff = "q9",
                              cols = rev(colormap(colormap = colormaps$density, nshades = 72, format = "hex",
                                                  alpha = 1, reverse = FALSE)), combine = F)
        
        p_gene <- lapply(X = p_gene, 
                         FUN = function(x) 
                                 x + theme(plot.title = element_text(face = 'plain'),
                                           panel.border = element_rect(colour = "black", fill=NA, size=1, linetype = 1),
                                           axis.line=element_blank(),
                                           axis.title.x=element_blank(),
                                           axis.title.y=element_blank(),
                                           legend.position = 'none',
                                           aspect.ratio = 1)
        )
        CombinePlots(plots = p_gene, ncol = ncol)
}

find_markers <- function(dataset, use.par = F, ncores = 4, min.diff.pct = -Inf) {
        
        DefaultAssay(object = dataset) <- "RNA"
        
        if (use.par == F) {
                
                markers <- list()
                
                for (i in 0:(length(levels(Idents(dataset)))-1)) {
                        
                        m <- FindConservedMarkers(object = dataset, ident.1 = i, grouping.var = "group", verbose = T, min.diff.pct = min.diff.pct)
                        m %<>%
                                add_column(feature = rownames(m), .before = 1) %>%
                                add_column(cluster = i, .after = 1)
                        
                        markers[[i+1]] <- m
                }
                
        } else {
                
                registerDoParallel(cores = ncores)
                
                markers <- foreach(i = 0:(length(levels(Idents(dataset)))-1)) %dopar% {
                        
                        m <- FindConservedMarkers(object = dataset, ident.1 = i, grouping.var = "group", verbose = T, min.diff.pct = min.diff.pct)
                        m %<>%
                                add_column(feature = rownames(m), .before = 1) %>%
                                add_column(cluster = i, .after = 1)
                        m
                }
        }
        markers <- do.call("rbind", markers)
        
        rownames(markers) <- NULL
        markers <- as_tibble(markers)
        
        return(markers)
}

get_top_genes <- function(dataset, markers, n) {
        
        int_features <- rownames(dataset@assays$integrated@scale.data)
        
        a <- str_subset(colnames(markers), "logFC")
        
        df <- markers %>%
                filter(feature %in% int_features) %>%
                mutate(logFC = (!!sym(a[1])) + (!!sym(a[2]))) %>%
                arrange(desc(logFC)) %>%
                group_by(cluster) %>%
                filter(row_number() <= n) %>%
                arrange(cluster)
        
        return(df$feature)
}

plot_heatmap <- function(dataset, markers, nfeatures) {
        
        df <- as_tibble(cbind(colnames(dataset), dataset$seurat_clusters, dataset$group))
        colnames(df) <- c("barcode","cluster","group")
        df$cluster <- as.numeric(df$cluster)
        df %<>%
                arrange(cluster, group)
        
        genes <- get_top_genes(dataset, markers, nfeatures)
        
        p_heat <- DoHeatmap(object = dataset, assay = 'integrated', features = genes, 
                            group.bar = F, cells = df$barcode, raster = F, draw.lines = F)
        
        p_pos_y <- ggplot_build(plot = p_heat)$layout$panel_params[[1]]$y.range
        
        ncol <- length(levels(Idents(dataset)))
        
        pal1 <- get_palette(ncol = ncol)
        col1 <- pal1[as.numeric(df$cluster)]
        
        pal2 <- c('#7bccc4','#f03b20')
        col2 <- pal2[as.numeric(factor(df$group))]
        
        p_heat + 
                annotation_raster(t(col2), -Inf, Inf, max(p_pos_y)+0.5, max(p_pos_y)+1.5) +
                annotation_raster(t(col1), -Inf, Inf, max(p_pos_y)+2, max(p_pos_y)+3) +
                coord_cartesian(ylim = c(0, max(p_pos_y)+4), clip = 'off') +
                scale_fill_gradient2(low = '#377eb8', high = '#e41a1c', mid = 'white', midpoint = 0) +
                guides(colour="none")
        
}

rename_cluster <- function(dataset, labels) {
        
        if (length(labels) != length(levels(Idents(dataset)))) {
                stop("Length of new names must be the same with old names.")
        } else {
                current.name <- levels(Idents(dataset))
                new.name <- labels
                
                Idents(dataset) <- plyr::mapvalues(x = Idents(dataset), from = current.name, to = new.name)
                return(dataset)
        }
}

plot_stat <- function(dataset, plot_type, 
                      group_levels, cluster_levels,
                      self_set_color = F,
                      self_colors,
                      group_colors = c('#7bccc4','#f03b20'),
                      palette = c("Paired", "Set2"),
                      plot_ratio = 1) {
        
        stat <- as_tibble(cbind(group = as.character(dataset$group), cluster = as.character(Idents(dataset))))
        stat %<>%
                mutate(group = factor(group, levels = group_levels),
                       cluster = factor(cluster, levels = cluster_levels)) %>%
                group_by(group, cluster) %>%
                summarise(n = n()) %>%
                mutate(freq = n / sum(n))
        
        cluster_colors <- if(self_set_color) self_colors else (get_palette(length(levels(Idents(dataset))), palette = palette))
        
        thm <- theme(aspect.ratio = plot_ratio,
                     axis.title.x = element_blank())
        
        switch(plot_type,
               group_count = stat %>%
                       group_by(group) %>%
                       summarise(sum(n)) %>%
                       ggplot() +
                       geom_col(aes(x = group, y = `sum(n)`, fill = group)) +
                       geom_text(aes(x = group, y = `sum(n)`, label = `sum(n)`), 
                                 vjust = -1) +
                       scale_fill_manual(values = c('#7bccc4','#f03b20'), name = "Group") + 
                       labs(y = "Counts") + thm,
               
               cluster_count = stat %>%
                       group_by(cluster) %>%
                       summarise(sum(n)) %>%
                       ggplot() +
                       geom_col(aes(x = cluster, y = `sum(n)`, fill = cluster)) +
                       geom_text(aes(x = cluster, y = `sum(n)`, label = `sum(n)`), 
                                 vjust = -1) +
                       scale_fill_manual(values = cluster_colors, name = "Clusters") + 
                       labs(y = "Counts") + 
                       theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
                       thm,
               
               prop_fill = ggplot(stat) + 
                       geom_bar(aes(x = group, y = freq, fill = cluster), position = "fill", stat = "identity") +
                       scale_y_continuous(labels = scales::percent) +
                       scale_fill_manual(values = cluster_colors, name = "Clusters") +
                       labs(y = "Proportion") + thm,
               
               prop_diverge = stat %>%
                       mutate(cluster = fct_rev(cluster)) %>%
                       mutate(freq = ifelse(group == group_levels[1], -freq, freq)) %>%
                       mutate(freq = round(freq, 3)) %>%
                       ggplot() + 
                       geom_bar(aes(x=cluster, y = freq, fill = group), stat = "identity") +
                       geom_text(aes(x=cluster, y = freq + 0.03 * sign(freq), label = percent(abs(freq), digits = 1))) +
                       coord_flip() +
                       scale_fill_manual(values = group_colors, name = "Stages") +
                       scale_y_continuous(breaks = pretty(c(stat$freq, -stat$freq)),
                                          labels = scales::percent(abs(pretty(c(stat$freq, -stat$freq))))) +
                       labs(x = NULL, y = "Proportion"),
               
               prop_multi = stat %>%
                       mutate(freq = round(freq, 3)) %>%
                       ggplot() + 
                       geom_bar(aes(x = group, y = freq, fill = group), stat = "identity") +
                       geom_text(aes(x = group, y = freq, label = scales::percent(freq)), vjust = -0.5) +
                       scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), labels = scales::percent_format(accuracy = 0.1)) +
                       facet_wrap(~cluster, ncol = 4, scales = "free") +
                       scale_fill_manual(values = group_colors, name = "Stages") +
                       labs(x = NULL, y = "Proportion"),
               
               stop("Unknown plot type")
        )
        
        
}

find_diff_genes <- function(dataset, clusters, groups, logfc = 0.25) {
        
        dataset$celltype.group <- paste(Idents(object = dataset), dataset$group, sep = "_")
        Idents(object = dataset) <- "celltype.group"
        
        
        de <- list()
        
        for (i in seq(length(clusters))) {
                
                d <- FindMarkers(dataset, 
                                  ident.1 = paste(clusters[i], groups[2], sep = "_"),
                                  ident.2 = paste(clusters[i], groups[1], sep = "_"),
                                  logfc.threshold = logfc,
                                  assay = "RNA")
                d %<>%
                        add_column(feature = rownames(d), .before = 1) %>%
                        add_column(cluster = clusters[i], .after = 1)
                
                de[[i]] <- as_tibble(d)
        }
        
        de <- do.call("rbind", de)
        
        de

}

test_GSEA <- function(diff, clusters, pathway) {
        
        gsea_res <- list()
        
        for (i in seq(length(clusters))) {
                
                data <- diff %>%
                        filter(cluster == clusters[i]) %>%
                        right_join(mm_hs, by = c("feature" = "mouse")) %>%
                        replace_na(list(avg_logFC = 0)) %>%
                        distinct(human, .keep_all = T) %>%
                        arrange(desc(avg_logFC))
                
                l <- data$avg_logFC
                names(l) <- data$human
                
                res <- fgsea(pathways = pathway, l, minSize = 15, maxSize = 500, nperm = 100000)
                
                res %<>%
                        add_column(cluster = clusters[i], .before = 1)
                
                gsea_res[[i]] <- res
        }
        
        gsea_res <- do.call("rbind", gsea_res)
        
        return(as_tibble(gsea_res))
        
}

plot_GSEA <- function(gsea_res, pattern = "HALLMARK_", p_cutoff = 0.05, levels) {
        
        gsea_res %>%
                mutate(pathway = str_remove(string = pathway, pattern = pattern)) %>%
                mutate(color = as.factor((padj < p_cutoff) * (ifelse(NES > 0, 1, -1)))) %>%
                ggplot(aes(x = factor(pathway), y = factor(cluster, levels = levels))) + 
                geom_point(aes(size = abs(NES), color = color)) +
                scale_color_manual(name = '', 
                                   values = c('dodgerblue1','grey','red'),
                                   labels = c('Down-regulation','Non-significant','Up-regulation')) +
                coord_flip() +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank())
        
}


# Test


gfp.combined <- gfp_combined

gfp.combined$celltype.group <- paste(Idents(object = gfp.combined), gfp.combined$group, sep = "_")
Idents(object = gfp.combined) <- "celltype.group"
Idents(object = gfp.combined) <- "celltype"
        
de <- FindMarkers(gfp.combined, 
                  ident.1 = paste(gfp_levels[1], stages[2], sep = "_"),
                  ident.2 = paste(gfp_levels[1], stages[1], sep = "_"),
                  logfc.threshold = 0,
                  assay = "RNA")

DefaultAssay(gfp.combined) <- "RNA"

test <- find_diff_genes(gfp_combined, gfp_levels, stages, logfc = 0.3)


hist(de$avg_logFC, breaks = 100)
hist(-log10(de$p_val_adj) * sign(de$avg_logFC), breaks = 100)
cor(de$avg_logFC, -log(de$p_val_adj) * sign(de$avg_logFC))



DoHeatmap(object = gfp_combined, assay = 'integrated', features =  get_top_genes(gfp_combined, gfp_markers, 6), 
          group.bar = F, raster = F, draw.lines = F)

stat <- as_tibble(cbind(group = as.character(cd45_combined$group), cluster = as.character(Idents(cd45_combined))))
stat %<>%
        mutate(group = factor(group, levels = stages),
               cluster = factor(cluster, levels = cd45_levels)) %>%
        group_by(group, cluster) %>%
        summarise(n = n()) %>%
        mutate(freq = n / sum(n))

stat %<>%
        mutate(freq = round(freq, 3))

ggplot(stat) + 
        geom_bar(aes(x = group, y = freq, fill = group), stat = "identity") +
        geom_text(aes(x = group, y = freq, label = scales::percent(freq)), vjust = -0.5) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), labels = scales::percent_format(accuracy = 0.1)) +
        facet_wrap(cluster~., scales = "free", ncol = 4) +
        scale_fill_manual(values = c('#7bccc4','#f03b20'), name = "Stages") +
        labs(x = NULL, y = "Proportion")

stat$cluster

