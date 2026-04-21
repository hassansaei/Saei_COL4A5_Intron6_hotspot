## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
library(BiocParallel)
register(SerialParam())

## ----eval=FALSE---------------------------------------------------------------
# score <- sum(colSums(E[p, ])**2) / length(p)

## ----message=FALSE------------------------------------------------------------
library(GEOquery)
library(limma)

gse200250 <- getGEO("GSE200250", AnnotGPL = TRUE)[[1]]

es <- gse200250
es <- es[, grep("Th2_", es$title)]
es$time <- as.numeric(gsub(" hours", "", es$`time point:ch1`))
es <- es[, order(es$time)]

exprs(es) <- normalizeBetweenArrays(log2(exprs(es)), method="quantile")

es <- es[order(rowMeans(exprs(es)), decreasing=TRUE), ]
es <- es[!duplicated(fData(es)$`Gene ID`), ]
rownames(es) <- fData(es)$`Gene ID`
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

fData(es) <- fData(es)[, c("ID", "Gene ID", "Gene symbol")]

es <- es[head(order(rowMeans(exprs(es)), decreasing=TRUE), 12000), ]
head(exprs(es))

## -----------------------------------------------------------------------------
library(msigdbr)
pathwaysDF <- msigdbr(species="mouse", collection="H")
pathways <- split(as.character(pathwaysDF$ncbi_gene), pathwaysDF$gs_name)

## ----message=FALSE------------------------------------------------------------
library(fgsea)
set.seed(1)
gesecaRes <- geseca(pathways, exprs(es), minSize = 15, maxSize = 500)

## -----------------------------------------------------------------------------
head(gesecaRes, 10)

## ----fig.width=10, fig.height=4, out.width="100%"-----------------------------
plotCoregulationProfile(pathway=pathways[["HALLMARK_E2F_TARGETS"]], 
                        E=exprs(es), titles = es$title, conditions=es$`time point:ch1`)

## ----fig.width=10, fig.height=4, out.width="100%"-----------------------------
plotCoregulationProfile(pathway=pathways[["HALLMARK_HYPOXIA"]], 
                        E=exprs(es), titles = es$title, conditions=es$`time point:ch1`)



## ----fig.width=10, fig.height=6, out.width="100%"-----------------------------
plotGesecaTable(gesecaRes |> head(10), pathways, E=exprs(es), titles = es$title)

## -----------------------------------------------------------------------------
E <- t(base::scale(t(exprs(es)), scale=FALSE))
pcaRev <- prcomp(E, center=FALSE)
Ered <- pcaRev$x[, 1:10]
dim(Ered)

## -----------------------------------------------------------------------------
set.seed(1)
gesecaResRed <- geseca(pathways, Ered, minSize = 15, maxSize = 500, center=FALSE)
head(gesecaResRed, 10)

## ----fig.width=4, fig.height=4------------------------------------------------
library(ggplot2)
ggplot(data=merge(gesecaRes[, list(pathway, logPvalFull=-log10(pval))],
                  gesecaResRed[, list(pathway, logPvalRed=-log10(pval))])) +
    geom_point(aes(x=logPvalFull, y=logPvalRed)) +
    coord_fixed() + theme_classic()

## -----------------------------------------------------------------------------
suppressMessages(library(Seurat))

## ----fig.width=8, fig.height=3.5----------------------------------------------
obj <- readRDS(url("https://alserglab.wustl.edu/files/fgsea/GSE116240.rds"))
obj

newIds <- c("0"="Adventitial MF",
            "3"="Adventitial MF",
            "5"="Adventitial MF",
            "1"="Intimal non-foamy MF",
            "2"="Intimal non-foamy MF",
            "4"="Intimal foamy MF",
            "7"="ISG+ MF",
            "8"="Proliferating cells",
            "9"="T-cells",
            "6"="cDC1",
            "10"="cDC2",
            "11"="Non-immune cells")

obj <- RenameIdents(obj, newIds)

DimPlot(obj) + ggplot2::coord_fixed()

## -----------------------------------------------------------------------------
obj <- SCTransform(obj, verbose = FALSE, variable.features.n = 10000)

## -----------------------------------------------------------------------------
length(VariableFeatures(obj)) # make sure it's a full gene universe of 10000 genes
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE,
                rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)

E <- obj@reductions$pca.rev@feature.loadings

## -----------------------------------------------------------------------------
library(msigdbr)

pathwaysDF <- msigdbr(species="mouse", collection="C2", subcollection = "CP:KEGG_LEGACY")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

## -----------------------------------------------------------------------------
set.seed(1)
gesecaRes <- geseca(pathways, E, minSize = 5, maxSize = 500, center = FALSE, eps=1e-100)

head(gesecaRes, 10)

## ----fig.width=12, fig.height=7, out.width="100%"-----------------------------
topPathways <- gesecaRes[, pathway] |> head(4)
titles <- sub("KEGG_", "", topPathways)

ps <- plotCoregulationProfileReduction(pathways[topPathways], obj,
                                       title=titles,
                                       reduction="tsne")
cowplot::plot_grid(plotlist=ps[1:4], ncol=2)

## ----fig.width=5, fig.height=3.5, out.width="50%"-----------------------------
plotCoregulationProfileReduction(pathways$KEGG_LYSOSOME, 
                               obj,
                               title=sprintf("KEGG_LYSOSOME (pval=%.2g)",
                                             gesecaRes[match("KEGG_LYSOSOME", pathway), pval]),
                               reduction="tsne")

## ----message=FALSE------------------------------------------------------------
library(Seurat)

obj <- readRDS(url("https://alserglab.wustl.edu/files/fgsea/275_T_ST_Seurat.rds"))

## -----------------------------------------------------------------------------
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE, variable.features.n = 10000)

obj <- RunPCA(obj, assay = "SCT", verbose = FALSE,
                rev.pca = TRUE, reduction.name = "pca.rev",
              reduction.key="PCR_", npcs = 50)

E <- obj@reductions$pca.rev@feature.loadings

## -----------------------------------------------------------------------------
library(msigdbr)
pathwaysDF <- msigdbr(species="human", collection="H")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)

## -----------------------------------------------------------------------------
set.seed(1)
gesecaRes <- geseca(pathways, E, minSize = 15, maxSize = 500, center = FALSE)

head(gesecaRes, 10)

## ----fig.width=10, fig.height=7, out.width="100%"-----------------------------

topPathways <- gesecaRes[, pathway] |> head(4)
titles <- sub("HALLMARK_", "", topPathways)

ps <- plotCoregulationProfileSpatial(pathways[topPathways], obj,
                                     title=titles, 
                                     pt.size.factor=2.5)
cowplot::plot_grid(plotlist=ps, ncol=2)

## ----fig.width=5, fig.height=3.5, out.width="50%"-----------------------------
plotCoregulationProfileSpatial(pathways$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
                               obj,
                               pt.size.factor=2.5,
                               title=sprintf("HALLMARK_OXIDATIVE_PHOSPHORYLATION\n(pval=%.2g)",
                                             gesecaRes[
                                                 match("HALLMARK_OXIDATIVE_PHOSPHORYLATION", pathway),
                                                 pval]))

## ----eval=FALSE---------------------------------------------------------------
# fldr <- ""         # Path to downloaded Xenium output
# annf_path <- ""    # Path to downloaded cell annotation CSV
# 
# # Load the Xenium data as a Seurat object
# xobj <- Seurat::LoadXenium(fldr, molecule.coordinates = FALSE)
# 
# # Remove assays not needed for downstream analysis
# for (nm in names(xobj@assays)[-1]) {
#     xobj@assays[[nm]] <- NULL
# }
# 
# # Add spatial coordinates to metadata
# coords <- data.table::as.data.table(xobj@images$fov@boundaries$centroids@coords)
# xobj@meta.data <- cbind(xobj@meta.data, coords)
# 
# # Integrate external cell annotations
# annot <- data.table::fread(annf_path)
# xobj <- xobj[, rownames(xobj@meta.data) %in% annot$cell_id]
# xobj@meta.data$annotation <- annot[match(rownames(xobj@meta.data), annot$cell_id), ]$group
# xobj@meta.data$annotation <- factor(xobj@meta.data$annotation)
# 
# 
# 
# # Filter out low-quality cells
# xobj <- subset(xobj, subset = nCount_Xenium > 50)
# 
# 
# # Optional: Crop the region and subsample the object
# # These steps reduce size and complexity for a more compact vignette
# # Helps keep the vignette lightweight and quick to render
# xobj <- xobj[, xobj$x < 3000 & xobj$y < 4000]
# set.seed(1)
# xobj <- xobj[, sample.int(ncol(xobj), 40000)]
# 
# # Normalize and scale data
# xobj <- NormalizeData(xobj, normalization.method = "LogNormalize", scale.factor = 1000, verbose = FALSE)
# xobj <- FindVariableFeatures(xobj, nfeatures = 2000, verbose = FALSE)
# xobj <- ScaleData(xobj, verbose = FALSE)
# 
# # Perform PCA and UMAP dimensionality reduction
# xobj <- RunPCA(xobj, verbose = FALSE)
# xobj <- RunUMAP(xobj, dims = 1:20)

## ----eval=FALSE---------------------------------------------------------------
# xobj@reductions$pca <- NULL
# xobj@assays$Xenium@layers$scale.data <- NULL
# xobj@assays$Xenium@layers$data <- NULL

## -----------------------------------------------------------------------------
xobj <- readRDS(url("https://alserglab.wustl.edu/files/fgsea/xenium-human-ovarian-cancer.rds"))

xobj <- NormalizeData(xobj, verbose = FALSE)
xobj <- ScaleData(xobj, verbose = FALSE)

## ----fig.width=10, fig.height=5, out.width="100%"-----------------------------
p1 <- ImageDimPlot(xobj, group.by = "annotation",
                   dark.background = FALSE, size = 0.5, 
                   flip_xy = TRUE) + 
    theme(legend.position = "none")

p1 <- suppressMessages(p1 + coord_flip() + scale_x_reverse() + theme(aspect.ratio = 1.0))

p2 <- DimPlot(xobj, group.by = "annotation", raster = FALSE, pt.size = 0.5, stroke.size = 0.0) + 
    coord_fixed()

p1 | p2

## -----------------------------------------------------------------------------
xobj <- RunPCA(xobj, verbose = F, rev.pca = TRUE, reduction.name = "pca.rev", 
               reduction.key="PCR_", npcs=30)

E <- xobj@reductions$pca.rev@feature.loadings

gsets <- msigdbr(species="human", collection="H")
gsets <- split(gsets$gene_symbol, gsets$gs_name)

set.seed(1)
gesecaRes <- geseca(gsets, E, minSize = 15, maxSize = 500, center = FALSE, eps = 1e-100)

head(gesecaRes, 10)

## ----fig.width=10, fig.height=7, out.width="100%"-----------------------------
topPathways <- gesecaRes[, pathway] |> head(6)
titles <- sub("HALLMARK_", "", topPathways)
pvals <- gesecaRes[match(topPathways, pathway), pval]
pvals <- paste("p-value:", formatC(pvals, digits = 1, format = "e"))

titles <- paste(titles, pvals, sep = "\n")

plots <- fgsea::plotCoregulationProfileReduction(
    gsets[topPathways], xobj, 
    reduction = "umap", title = titles,
    raster = TRUE, raster.dpi = c(500, 500),
    pt.size = 2.5
)

plots <- lapply(plots, function(p){
    p + theme(plot.title = element_text(hjust=0), text = element_text(size=10))
})


cowplot::plot_grid(plotlist = plots, ncol = 2)

## ----fig.width=10, fig.height=15, out.width="100%"----------------------------
imagePlots <- plotCoregulationProfileImage(
    gsets[topPathways], object = xobj, 
    dark.background = FALSE, 
    title = titles, 
    size=0.8
)

imagePlots <- lapply(imagePlots, function(p){
    suppressMessages(
        p + coord_flip() + scale_x_reverse() + 
            theme(plot.title = element_text(hjust=0), text = element_text(size=10), aspect.ratio = 1.0)
    )
})

cowplot::plot_grid(plotlist = imagePlots, ncol = 2)

## ----echo=TRUE----------------------------------------------------------------
sessionInfo()

