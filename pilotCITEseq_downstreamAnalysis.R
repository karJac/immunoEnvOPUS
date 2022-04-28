library(DoubletFinder)
library(Matrix)
library(Seurat)
library(readxl)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(parallel)
library(R.devices)
library(pbmcapply)
library(data.table)
library(schex)

wait(10000)

source("~/all/customFunctionsKJ.R")

setwd("/home/rstudio/all/Opus/redo_pilot")



################### HASHTAG OLIGOS DEMULTIPLEXING ###################

h5files <- list.files(pattern = ".h5")

well1 <- Read10X_h5(h5files[[1]])
well2 <- Read10X_h5(h5files[[2]])
well3 <- Read10X_h5(h5files[[3]])

UMI1 <- well1[[2]][24:27, ]
UMI2 <- well2[[2]][24:27, ]
UMI3 <- well3[[2]][24:27, ]


PROT1 <- well1[[2]][c(1:21), ]
PROT2 <- well2[[2]][c(1:21), ]
PROT3 <- well3[[2]][c(1:21), ]


s1 <- CreateSeuratObject(well1[[1]])
s2 <- CreateSeuratObject(well2[[1]])
s3 <- CreateSeuratObject(well3[[1]])

s1[["ADT"]] <- CreateAssayObject(counts = PROT1)
s2[["ADT"]] <- CreateAssayObject(counts = PROT2)
s3[["ADT"]] <- CreateAssayObject(counts = PROT3)

s1[["HTO"]] <- CreateAssayObject(counts = UMI1)
s2[["HTO"]] <- CreateAssayObject(counts = UMI2)
s3[["HTO"]] <- CreateAssayObject(counts = UMI3)

s1 <- NormalizeData(s1, assay = "HTO", normalization.method = "CLR")
s2 <- NormalizeData(s2, assay = "HTO", normalization.method = "CLR")
s3 <- NormalizeData(s3, assay = "HTO", normalization.method = "CLR")

# s1 <- NormalizeData(s1, assay = "ADT", normalization.method = "CLR", margin = 2)
# s2 <- NormalizeData(s2, assay = "ADT", normalization.method = "CLR", margin = 2)
# s3 <- NormalizeData(s3, assay = "ADT", normalization.method = "CLR", margin = 2)

s1 <- HTODemux(s1, assay = "HTO", positive.quantile = 0.98)
s2 <- HTODemux(s2, assay = "HTO", positive.quantile = 0.98)
s3 <- HTODemux(s3, assay = "HTO", positive.quantile = 0.98)

Idents(s1) <- "HTO_maxID"
Idents(s2) <- "HTO_maxID"
Idents(s3) <- "HTO_maxID"

# QC tests

table(s1$HTO_classification.global)

Idents(s1) <- "HTO_maxID"
RidgePlot(s1, assay = "HTO", features = rownames(s1@assays$HTO@data), ncol = 2)

Idents(s1) <- "HTO_classification.global"
VlnPlot(s1, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

png("HTO_Heatmap.png",width=1000, height = 500)
HTOHeatmap(s1, assay = "HTO", ncells = length(colnames(s3)))
dev.off()

##### demultiplexing #####

Idents(s1) <- "HTO_classification"
Idents(s2) <- "HTO_classification"
Idents(s3) <- "HTO_classification"

shamsList <- lapply(1:3, function(i) {
  tmpList <- lapply(1:2, function(j) {
    tmp <- subset(get(paste("s", i, sep = "")), idents = paste("Bag-", j, "-prot", sep = ""))
    return(tmp)
  })
})
shamsList <- unlist(shamsList)
shamsList <- mapply(function(i, j) {
  Idents(j) <- paste("sham", i, sep = "")
  Project(j) <- paste("sham", i, sep = "")
  return(j)
}, 1:length(shamsList), shamsList)

tumorsList <- lapply(1:3, function(i) {
  tmpList <- lapply(3:4, function(j) {
    tmp <- subset(get(paste("s", i, sep = "")), idents = paste("Bag-", j, "-prot", sep = ""))
    return(tmp)
  })
})
tumorsList <- unlist(tumorsList)
tumorsList <- mapply(function(i, j) {
  Idents(j) <- paste("tumor", i, sep = "")
  Project(j) <- paste("tumor", i, sep = "")
  return(j)
}, 1:length(tumorsList), tumorsList)


tumorsADT <- lapply(tumorsList,function(i){
     tmp <- i@assays$ADT@counts
     colnames(tmp) <- colnames(i)
     return(tmp)
})

shamsADT <- lapply(shamsList,function(i){
     tmp <- i@assays$ADT@counts
     colnames(tmp) <- colnames(i)
     return(tmp)
})

# I tak musisz na nowo poszukac dubletow, bo moga byc w obrebach tego samego hasha. Wyszukaj je wiec DoubletFinderem ale nie wyrzucaj od razu, zobacz gdzie sie ukladaja na UMAPIE po clustrowaniu.

samplesList <- c(shamsList[1:6],tumorsList[1:6])

samplesListQC <- pbmclapply(samplesList, function(sample) {
  sample <- PercentageFeatureSet(sample, "^mt-", col.name = "percent_mito")
  sample <- PercentageFeatureSet(sample, "^Rp[sl]", col.name = "percent_ribo")
  sample <- PercentageFeatureSet(sample, "^Hb[^(p)]", col.name = "percent_hb")
  sample <- PercentageFeatureSet(sample, "Pecam1|Pf4", col.name = "percent_plat")

  feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb", "percent_plat")
  sample$orig.ident <- Project(sample)

  nCountFeaturePlot <- FeatureScatter(sample, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)

  # filter the extreme outliers to get properl violin plot
  selected_mito <- WhichCells(sample, expression = percent_mito < 15)
  test <- subset(sample, cells = selected_mito)

  QCplot <- VlnPlot(test, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()



  exp <- sample@assays$RNA@counts

  # mouse_mito <- read_excel("Mouse.MitoCarta3.0.xls") # lepsze przyrownanie niż "^mt-" ale nie skonczylem go pisac
  # mouseMito <- mouse_mito[,3]


  selected_c <- WhichCells(sample, expression = nFeature_RNA > 200)
  # dodaj usuwanie komórek z mniejsza niz 3k? UMI
  selected_f <- rownames(sample)[Matrix::rowSums(sample) > 20]
  sample.filt <- subset(sample, features = selected_f, cells = selected_c)
  # sample.filt[["ADT"]] <- CreateAssayObject(counts = sample@assays$ADT@data[,which(colnames(sample@assays$ADT@data) %in% selected_c)])
  dim(sample.filt)
  remove(sample)

  selected_mito <- WhichCells(sample.filt, expression = percent_mito < 5.5) # normaly ~10-15%
  selected_ribo <- WhichCells(sample.filt, expression = percent_ribo > 3.5) # normaly ~5 (cuz it correlates with high levels of mito genes)
  tryCatch(selected_hb <- WhichCells(sample.filt, expression = percent_hb < 2), error = function(e) NULL)
  tryCatch(selected_plat <- WhichCells(sample.filt, expression = percent_plat < 1), error = function(e) NULL)

  sample.filt <- subset(sample.filt, cells = selected_mito)
  sample.filt <- subset(sample.filt, cells = selected_ribo)
  if (exists("selected_hb") == TRUE) {
    sample.filt <- subset(sample.filt, cells = selected_hb)
  }
  if (exists("selected_plat") == TRUE) {
    sample.filt <- subset(sample.filt, cells = selected_plat)
  }

  dim(sample.filt)

  table(sample.filt$orig.ident)

  QCfilteredPlot <- VlnPlot(sample.filt, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) +
    NoLegend()

  # Compute the relative expression of each gene per cell Use sparse matrix
  # operations, if your dataset is large, doing matrix devisions the regular way
  # will take a very long time.
  par(mar = c(4, 8, 2, 1))
  C <- sample.filt@assays$RNA@counts
  C <- Matrix::t(Matrix::t(C) / Matrix::colSums(C)) * 100
  most_expressed <- order(apply(C, 1, median), decreasing = T)[30:1]
  tmp <- data.frame(t(C[most_expressed, ]))
  tmp2 <- pivot_longer(tmp, cols = colnames(tmp))
  tmp2$name <- factor(tmp2$name, levels = colnames(tmp), ordered = TRUE)
  highestExprsPlot <- ggplot(tmp2, aes(x = value, y = name, fill = name)) +
    geom_boxplot()

  R.devices::suppressGraphics({
    png(paste("QCplot_", Project(sample.filt), ".png", sep = ""), height = 1200, width = 1200)
    print(plot_grid(QCplot, nCountFeaturePlot, QCfilteredPlot, highestExprsPlot,
      nrow = 2,
      labels = c(paste("QC for sample: ", Project(sample.filt)), "nCountFeature", "QC after filtering", "most expressed"),
      vjust = 0.85, greedy = FALSE
    )) # improve annotations of graphs because they are overlapping with figures
    dev.off()
  })
  return(sample.filt)
}, mc.cores = length(samplesList))



samplesListFilt <- pbmclapply(samplesListQC, function(sample.filt) {

  tmp <- sample.filt
  # Filter MALAT1
  sample.filt <- sample.filt[!grepl("Malat1", rownames(sample.filt)), ]
  # Filter Mitocondrial
  # sample.filt <- sample.filt[!grepl("^MT-", rownames(sample.filt)), ]
  # Filter Ribosomal
  # sample.filt <- sample.filt[!grepl("^RP[SL][[:digit:]]", rownames(sample.filt)), ]
  # Filter Ribosomal rRNA
  sample.filt <- sample.filt[!grepl("rRNA", ignore.case = TRUE, rownames(sample.filt)), ]
  # sample.filt[["ADT"]] <- CreateAssayObject(counts = tmp@assays$ADT@data)
  
  return(sample.filt)
}, mc.cores = length(samplesList))


httr::set_config(httr::config(ssl_verifypeer = FALSE))
m.s.genes <- NULL
while(is.null(m.s.genes)){  #infinite loop that will run until converhumangenelist wont crash (it crash because of connction issues with ensmeble)
     try(m.s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes))
}

m.g2m.genes <- NULL
while(is.null(m.g2m.genes)){
try(m.g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes))
}

samplesListNorm <- mclapply(samplesListFilt, function(sample.filt) {

  # Before running CellCycleScoring the data need to be normalized and
  # logtransformed.
  #sample.filt@active.assay <- "RNA"
  sample.filt <- NormalizeData(sample.filt)


  sample.filt <- CellCycleScoring(object = sample.filt, g2m.features = m.g2m.genes,
                                s.features = m.s.genes)
  cellCycleGraph <- VlnPlot(sample.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
                            ncol = 4, pt.size = 0.1)

R.devices::suppressGraphics({
        png(paste("CellCycleplot_",Project(sample.filt),".png",sep=''))
        print(cellCycleGraph)
        dev.off()
})

  sample.filt <- FindVariableFeatures(sample.filt, verbose = F)
  sample.filt <- ScaleData(sample.filt,
    vars.to.regress = c("nFeature_RNA", "percent_mito"),
    verbose = F
  )
  sample.filt <- RunPCA(sample.filt, verbose = F, npcs = 40)
  sample.filt <- RunUMAP(sample.filt, dims = 1:40, verbose = F)

  return(sample.filt)
}, mc.cores = length(samplesListFilt))



samplesAfterDF2 <- pbmclapply(samplesListNorm, function(sample.filt) {
  sweep.res <- paramSweep_v3(sample.filt)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  pK <- as.numeric(as.character(bcmvn$pK))
  BCmetric <- bcmvn$BCmetric
  pK_choose <- pK[which(BCmetric %in% max(BCmetric))]

  par(mar = c(5, 4, 4, 8) + 1, cex.main = 1.2, font.main = 2)
  # plot(x = pK, y = BCmetric, pch = 16,type="b",  #visualtion of BCmetrics (but function finds the maximum automaticly)
  #     col = "blue",lty=1)
  # abline(v=pK_choose,lwd=2,col='red',lty=2)
  # title("The BCmvn distributions")
  # text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

  # define the expected number of doublet cellscells.
  nExp <- round(ncol(sample.filt) * 0.025) # we expect 2.5% doublets (10X doublets prop) - doublets we discarded at demultiplexing
  sample.filt <- doubletFinder_v3(sample.filt, pN = 0.25, pK = pK_choose, nExp = nExp, PCs = 1:10)
  n <- dim(sample.filt@meta.data)[[2]] - 1
  m <- dim(sample.filt@meta.data)[[2]]
  colnames(sample.filt@meta.data)[[n]] <- "pANN"
  colnames(sample.filt@meta.data)[[m]] <- "DoubletFinder"
  
  return(sample.filt)
}, mc.cores = length(samplesListNorm))



# Merge datasets into one single seurat object

for (i in c(1:6)){
     tumorsADT[[i]]<- tumorsADT[[i]][,which(colnames(tumorsADT[[i]]) %in% colnames(samplesAfterDF2[[i+6]]))]
     shamsADT[[i]]<- shamsADT[[i]][,which(colnames(shamsADT[[i]]) %in% colnames(samplesAfterDF2[[i]]))]
}
for (i in c(1:6)){
     samplesAfterDF2[[i]][["ADT"]]<- CreateAssayObject(counts = shamsADT[[i]])
}
for (i in 1:6){
     samplesAfterDF2[[i+6]][["ADT"]]<- CreateAssayObject(counts = tumorsADT[[i]])
}


alldata <- merge(samplesAfterDF2[[1]], c(samplesAfterDF2[2:length(samplesAfterDF2)]))
DefaultAssay(alldata) <- 'RNA'
alldata <- SCTransform(alldata,ncells=dim(alldata)[[2]]) %>% RunPCA(npc=50)
# alldata <- JackStraw(alldata, reduction = "pca", dims = 100, prop.freq = 0.1) # find the amount of PC u want to use in UMAP
# alldata <- ScoreJackStraw(alldata, dims = 1:100)
# JackStrawPlot(alldata, dims = 1:100)
     
DefaultAssay(alldata) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting
VariableFeatures(alldata) <- rownames(alldata[["ADT"]])
alldata <- NormalizeData(alldata,assay = "ADT", normalization.method = "CLR", margin = 2) %>% ScaleData() %>% RunPCA(reduction.name = 'apca')

DefaultAssay(alldata) <- 'SCT'
alldata <- FindMultiModalNeighbors(
     alldata, reduction.list = list("pca", "apca"), 
     dims.list = list(1:47, 1:20) #You have to choose number of PC's 
                                                                       #from jackstrawplot for RNA and i think that all avaible PC's from protein
)

#UMAP on RNA+prot
alldata <- RunUMAP(alldata, nn.name = "weighted.nn", reduction.name = "aa.umap", reduction.key = "wnnUMAP_", seed.use=69,  min.dist = 0.3)
#UMAP on prot
alldata <- RunUMAP(alldata, reduction = 'apca', dims = 1:18, assay = 'ADT', 
                   reduction.name = 'adt.umap', reduction.key = 'adtUMAP_' )
#UMAP on RNA
alldata <- RunUMAP(alldata, reduction = 'pca', dims = 1:47, assay = 'RNA', 
        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_', min.dist = 0.5, n.neighbors = 40, seed.use=69)

p1 <- DimPlot(alldata, reduction = 'rna.umap',    group.by = 'orig.ident', label = FALSE, 
              repel = TRUE, label.size = 2.5, ncol=3) + NoLegend()
p2 <- DimPlot(alldata, reduction = 'adt.umap',    group.by = 'orig.ident', label = FALSE, 
              repel = TRUE, label.size = 2.5) + NoLegend()
p3 <- DimPlot(alldata, reduction = 'wnn.umap', split.by = 'orig.ident', label = FALSE, 
              repel = TRUE, label.size = 2.5, ncol=6) + NoLegend()
p1 + p2 + p3

FeaturePlot(alldata, features=rownames(alldata@assays$ADT), ncol=6, max.cutoff="q99", reduction = 'wnn.umap')

DF.name <- colnames(alldata@meta.data)[grepl("DoubletFinder", colnames(alldata@meta.data))]

png("doublets.png", width = 1200, height = 500)
cowplot::plot_grid(
  ncol = 2, DimPlot(alldata, group.by = "orig.ident",reduction = "wnn.umap") + NoAxes(),
  DimPlot(alldata, group.by = "DoubletFinder",reduction = "wnn.umap") + NoAxes()
)
dev.off()

            
png("doubletsViolin.png", width = 1200, height = 500)
VlnPlot(alldata, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)
dev.off()

# alldata = alldata[, alldata@meta.data[, DF.name] == "Singlet"]




########################### Dimensionality reduction #############################


top20 <- head(VariableFeatures(alldata), 20)

png("volcanoPlot.png")
LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)
dev.off()

plot_grid(
  ncol = 3, DimPlot(alldata,
    reduction = "pca", group.by = "orig.ident",
    dims = 1:2
  ), DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 3:4),
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 5:6)
)
VizDimLoadings(alldata, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)

png("PCAcontributions.png")
ElbowPlot(alldata, reduction = "pca", ndims = 50) # amonut of variance explained by each PC
dev.off()

#### Check if your data needs integration
png("origIdent.png", height = 600, width = 1800)
DimPlot(alldata, reduction = "wnn.umap", group.by = "orig.ident") +
  ggplot2::ggtitle(label = "UMAP_on_PCA") +
  facet_wrap(~orig.ident, ncol=6)
dev.off()


saveRDS(alldata, "alldataBeforeClust.rds")






############################ Clustering ############################

library(pheatmap)
library(enrichR)
library(rafalib)
library(clustree)
library(leiden)


alldata@active.assay <- "SCT"


#alldata <- FindNeighbors(alldata, dims = 1:47, k.param = 60, prune.SNN = 1 / 15)
names(alldata@graphs)

# pheatmap(alldata@graphs$RNA_nn[1:200, 1:200],
#   col = c("white", "black"), border_color = "grey90",
#   legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2
# )

# Clustering with Leiden (algorithm 4)
for (res in c(4)) {
  alldata <- FindClusters(alldata, graph.name = names(alldata@graphs)[[2]], resolution = res, algorithm = 4) # 4 = Leiden ALGORITHM
}

# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only CCA_snn_res.XX - for each different
# resolution you test.

plot_grid(ncol = 3, DimPlot(alldata, reduction = "umap", group.by = "RNA_snn_res.0.5") +
  ggtitle("leibens_0.5"), DimPlot(alldata, reduction = "umap", group.by = "RNA_snn_res.1") +
  ggtitle("leibens_1"), DimPlot(alldata, reduction = "umap", group.by = "RNA_snn_res.2") +
  ggtitle("leibens_2"))

png("SCTclustree2.png", width = 1200)
clustree(alldata@meta.data, prefix = "wsnn_res.")
dev.off()

alldata <- SetIdent(alldata, value = "wsnn_res.1")

VlnPlot(alldata, features="nCount_RNA")

for (i in 1:24){print(paste(i,":"))
     print(median(alldata$nCount_RNA[which(alldata$wsnn_res.0.5==i)]))}

DimPlot(alldata, group.by = "wsnn_res.1",reduction = "wnn.umap", label=T) + NoAxes()

saveRDS(alldata, "SCTafterClustering.rds")


library(tidyr)

sel.clust <- "RNA_snn_res.1.5" # czemu nie dziala przy PLOT HEATMAP
alldata <- SetIdent(alldata, value = "wsnn_res.0.5")
table(alldata@active.ident)


# plot this clustering
plot_grid(ncol = 3, DimPlot(alldata, label = T) + NoAxes(), DimPlot(alldata, group.by = "orig.ident", shuffle = TRUE) + NoAxes())

alldata@active.assay <- "RNA"
genes <- unique(c("Tmem119", "Cd3d", "Cd14",  "P2ry12", "Lgals3", "Cd38", "Nkg7", "Ccl12", "Ccr2", "Cd79a", "Cd19", "Celiac1", "H2-Eb2", "Cd44", "Il1r2", "Cxcr2", "Fpr2", "Fcgr3a", "Cecam8", "Sox2", "Cd4", "Cd8a", "Cd24a", "Itgax","Ptprc","Cd3d","Nkg7", "H2-Eb2", "Cd44", "Il1r2", "Cd4", "Cd8a", "Itgax","Ptprc", "Itgam", "Foxp3", "Pdcd1", "Il2ra"))
eh <- function(alldata,genes=genes){FeaturePlot(alldata, reduction = "wnn.umap", features = genes, order = T, slot = "data", combine = T, ncol = 6)
}
eh(lymph)

FeaturePlot(alldata, reduction = "wnn.umap", features = c("Foxp3", "Cd25"), order = T, slot = "data", combine = T)

FeaturePlot(alldata, reduction = "wnn.umap", features = c("Cd11b", "Spp1", "Arg1"), order = T, slot = "data", combine = T)




table(Idents(alldata), alldata$orig.ident)
table(Idents(alldata), alldata$orig.ident) %>%
  prop.table(margin = 2) %>%
  round(digits = 4) * 100


