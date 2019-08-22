library(monocle)
library(colormap)
library(colorRamps)
library(Seurat)

# Helper functions

colfunc <- colorRampPalette(c("dodgerblue4","dodgerblue","skyblue1",'white',"#ffff33","#ff7f00","#e41a1c"))
colfunc2 <- colorRampPalette(c("#377eb8",'white',"#e41a1c"))

plotGeneHeatMap <- function(data, diff_table, topn = 50, cluster = 3) {
        
        genes <- head(as.character(diff_table[with(diff_table, order(qval)),]$gene_short_name), n = topn)
        
        plot_pseudotime_heatmap(data[genes,],
                                num_clusters = cluster,
                                cores = 1,
                                show_rownames = T,
                                return_heatmap = F,
                                hmcols = colfunc2(72),
                                cluster_rows = T)

}




# Load data 

gfp.combined <- readRDS('../GFP_V3/gfp.combined.rds')
gfp.data <- GetAssayData(gfp.combined, slot = 'counts')

gfp.pd <- new('AnnotatedDataFrame', data = gfp.combined@meta.data)

gfp.fd <- new('AnnotatedDataFrame', 
              data = data.frame(gene_short_name = row.names(gfp.data), row.names = rownames(gfp.data)))

# Construct monocle cds
gfp.cds <- newCellDataSet(gfp.data,
                              phenoData = gfp.pd,
                              featureData = gfp.fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

gfp.cds <- estimateSizeFactors(gfp.cds)
gfp.cds <- estimateDispersions(gfp.cds)

gfp.cds <- detectGenes(gfp.cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(gfp.cds),
                                    num_cells_expressed >= 10))

print(head(pData(gfp.cds)))

qplot(nCount_RNA, data = pData(gfp.cds), color = stage, geom =
              "density")

qplot(mito, data = pData(gfp.cds), color = celltype, geom =
              "density")

qplot(Size_Factor, data = pData(gfp.cds), color = stage, geom =
              "density")

# Subset

gfp.cds_epi <- gfp.cds[,row.names(subset(pData(gfp.cds),
                                     celltype %in% c('Basal 1','Luminal 1','Basal 2','Luminal 2')))]

gfp.cds_fi <- gfp.cds[,row.names(subset(pData(gfp.cds),
                                         celltype %in% c('Fibroblast')))]

# Myo-epithelial analysis


diff_test_res_epi <- differentialGeneTest(gfp.cds_epi[expressed_genes,],
                                      fullModelFormulaStr = "~stage")
ordering_genes_epi <- row.names(subset(diff_test_res_epi, qval < 10E-70))

gfp.cds_epi <- setOrderingFilter(gfp.cds_epi, ordering_genes_epi)

plot_ordering_genes(gfp.cds_epi)

gfp.cds_epi <- reduceDimension(gfp.cds_epi, max_components = 2,
                            method = 'DDRTree')

gfp.cds_epi <- orderCells(gfp.cds_epi, reverse = T)

plot_cell_trajectory(gfp.cds_epi, color_by = "stage", show_branch_points = F) + 
        scale_color_manual(values = c('#7bccc4','#f03b20'))

plot_cell_trajectory(gfp.cds_epi, color_by = "stage", show_branch_points = F, cell_size = 0.5) + 
        scale_color_manual(values = c('#7bccc4','#f03b20')) +
        facet_wrap(~stage, nrow = 2)

plot_cell_trajectory(gfp.cds_epi, color_by = "celltype", show_branch_points = F) + 
        scale_color_manual(name = '', values = c('#f03b20','#a1d99b','#6baed6','#fed976'))

plot_cell_trajectory(gfp.cds_epi, color_by = "stage", show_branch_points = F) + 
        scale_color_manual(values = c('#a1d99b','#f03b20'), 
                           labels = c('Early','Advanced'),
                           name = 'Stage') +
        theme(title = element_blank()) +
        facet_wrap(~celltype + stage, nrow = 2)

plot_cell_trajectory(gfp.cds_epi, color_by = "Pseudotime", show_branch_points = F) + 
        scale_color_colormap(colormap = colormaps$density, reverse = T) +
        facet_wrap(~celltype, nrow = 2)

plot_cell_trajectory(gfp.cds_epi, color_by = "Pseudotime", show_branch_points = F) + 
        scale_color_colormap(colormap = colormaps$density, reverse = T) +
        facet_wrap(~stage, nrow = 1)

plot_cell_trajectory(gfp.cds_epi, color_by = "Pseudotime", show_branch_points = F) + 
        scale_color_colormap(colormap = colormaps$density, reverse = T)


plotGeneHeatMap(gfp.cds_epi, diff_test_res_epi, topn=50, cluster = 2)



# Fibroblast analysis

diff_test_res_fi <- differentialGeneTest(gfp.cds_fi[expressed_genes,],
                                          fullModelFormulaStr = "~stage")
ordering_genes_fi <- row.names(subset(diff_test_res_fi, qval < 0.05))

gfp.cds_fi <- setOrderingFilter(gfp.cds_fi, ordering_genes_fi)

plot_ordering_genes(gfp.cds_fi)

gfp.cds_fi <- reduceDimension(gfp.cds_fi, max_components = 2,
                               method = 'DDRTree')

gfp.cds_fi <- orderCells(gfp.cds_fi, reverse = F)

plot_cell_trajectory(gfp.cds_fi, color_by = "stage", show_branch_points = F) + 
        scale_color_manual(values = c('#7bccc4','#f03b20')) +
        facet_wrap(~stage, nrow = 1)

plot_cell_trajectory(gfp.cds_fi, color_by = "celltype", show_branch_points = F) + 
        scale_color_manual(values = c('#f03b20','#a1d99b','#6baed6','#fed976')) +
        facet_wrap(~stage, nrow = 1)

plot_cell_trajectory(gfp.cds_fi, color_by = "stage", show_branch_points = F) + 
        scale_color_manual(values = c('#a1d99b','#f03b20')) +
        facet_wrap(~celltype + stage, nrow = 1)

plot_cell_trajectory(gfp.cds_fi, color_by = "Pseudotime", show_branch_points = F) + 
        scale_color_colormap(colormap = colormaps$density, reverse = T) +
        facet_wrap(~celltype, nrow = 2)

plot_cell_trajectory(gfp.cds_fi, color_by = "Pseudotime", show_branch_points = F) + 
        scale_color_colormap(colormap = colormaps$density, reverse = T) +
        facet_wrap(~stage, nrow = 1)



plotGeneHeatMap(gfp.cds_fi, diff_test_res_fi, cluster = 3)


plot_cell_trajectory(gfp.cds_fi, markers = c("Krt14", "Col1a1"), use_color_gradient = TRUE)




# marker

cds_expressed_genes <-  row.names(subset(fData(cds_fi),
                                          num_cells_expressed >= 10))
cds_filtered <- cds_fi[cds_expressed_genes,]
my_genes <- row.names(subset(fData(cds_filtered),
                             gene_short_name %in% c("Col13a1", "Krt14", "Lama2")))
cds_subset <- cds_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "stage")


print(head(fData(cds)))
print(head(pData(cds_epi)))

colnames(pData(cds))



plot_cell_trajectory(cds, color_by = "State")


# Test

mtx <- GetAssayData(gfp.combined, slot = 'counts')

test <- newCellDataSet(mtx,
                          phenoData = gfp.pd,
                          featureData = gfp.fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())

pd <- new("AnnotatedDataFrame", data = sample_sheet)


bks <- seq(-3.1,3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)

scale_color_colormap(discrete = F,colormap = colormaps$viridis, reverse = T)

colormap(colormap = colormaps$viridis, nshades = 72, format = "hex",
         alpha = 1, reverse = FALSE)

colfunc(72)
plot(rep(1,72),col=colfunc(72),pch=19,cex=3)



plot_pseudotime_heatmap(gfp.cds_fi[ordering_genes_fi,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap = F,
                        hmcols = colormap(colormap=colormaps$viridis, nshades=72),
                        cluster_rows = T)


