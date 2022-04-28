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




########################## Identification ##########################

alldata <- SetIdent(alldata, value = "wsnn_res.1")
metrics <-  c("nCount_RNA", "nFeature_RNA", "G2M.Score","S.Score", "percent_mito", "percent_ribo")

png("alldata_uninterestingVariation.png", width = 2600, height = 1600)
FeaturePlot(alldata, 
            reduction = "wnn.umap", 
            features = metrics,
            order = TRUE,
            min.cutoff = 'q10',
            label = T,
            ncol=3)
dev.off()



FeaturePlot(alldata, reduction = "wnn.umap", features = c("Cd14","Fcgr3","Cd19","Cd38"), order = T, slot = "data", combine = T)


Idents(alldata) <- alldata$annot3
markers <- FindAllMarkers(alldata, only.pos = TRUE, min.pct = 0.5)


#annot csv is made for human, so brutal translation to mouse genes is far from ideal, we lose a lot of genes (for example H2-...)
#that's why we use it only for genes that we didn't successfully map to annot2
annot <- read_csv("genesAnnotations.csv")

require("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = annot$gene_name, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)

annot <- annot[which(annot$gene_name %in% genesV2$HGNC.symbol),]
annot$gene_name <- genesV2[match(annot$gene_name,genesV2$HGNC.symbol),2]

annot2 <- read_tsv("GENE-DESCRIPTION-TSV_MGI_18.tsv", col_names= FALSE)
annot2 <- annot2[,2:3]
colnames(annot2) <- c("gene","desc")


Desc1 <- annot2[(annot2$gene %in% markers0.5$gene), ]
noDesc1 <- markers0.5[!(markers0.5$gene %in% annot2$gene), ]
Desc2 <-
     annot[(annot$gene_name) %in% noDesc$gene, c("gene_name", "description")]
noDesc2 <- noDesc[!(noDesc$gene %in% annot$gene_name), ]
Desc3 <- data.frame(noDesc2$gene, rep("No description available" , length(noDesc2$gene)))
colnames(Desc1) <- NULL
colnames(Desc2) <- NULL
colnames(Desc3) <- NULL
Desc <- rbind(data.frame(Desc1), data.frame(Desc2), data.frame(Desc3))
colnames(Desc) <- c("gene", "desc")

markers0.5$desc <-
     lapply(match(markers0.5$gene, Desc$gene), function(i) {
          return(Desc$desc[[i]])
     })

write_csv(markers0.5, "markers.csv")

rna <- alldata@assays$RNA@data
adt <- alldata@assays$ADT@data

macExp <- rbind(rna,adt)
macExp <- t(macExp)

metadata <- alldata@meta.data




############### SHORT ANNOTATIONS
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
anno<-getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
            filters = "ensembl_gene_id", values = data_norm@Dimnames[1] , mart = ensembl)

macExpFull <- cbind(metadata,umaps,macExp)
saveRDS(macExpFull,"macExp05_10_21.rds")

grey <- lapply(1:24, function(i){
     greyClusters(alldata,i)
})




#Predcition nr1 !!!!
library(scMCA)



test <- t(macExp)
mca_result <- scMCA(scdata = test, numbers_plot = 3)


################ RECLUSTERING #################


alldata <- readRDS("alldataAfterClust.rds")
annot <- read_csv("annot2_pilot.csv")

macExp <- read_csv("macExp_17.11.21.csv")

alldata$annot1 <- annot$cell_annotation_1[match(macExp$wsnn_res.1,annot$wsnn_res.0.5)]
alldata$annot2 <- annot$cell_annotation_2[match(macExp$wsnn_res.1,annot$wsnn_res.0.5)]
alldata$annot3 <- annot$cell_annotation_3[match(macExp$wsnn_res.1,annot$wsnn_res.0.5)]

DimPlot(alldata, group.by = "annot3",reduction = "wnn.umap", label=T) + NoAxes()

Idents(alldata) <- alldata$annot1


#LYMPHO#

lympho <- subset(alldata, idents = "Lymphoid")

#alldatacp <- alldata
alldata <- lympho

alldata@active.assay <- "RNA"

for (res in c(0.1, 0.25, 0.5, 1, 1.5, 2)) {
     alldata <- FindClusters(alldata, graph.name = names(alldata@graphs)[[2]], resolution = res, algorithm = 4) # 4 = Leiden ALGORITHM
}

alldata <- FindMultiModalNeighbors(
     alldata, reduction.list = list("pca", "apca"), 
     dims.list = list(1:47, 1:20), modality.weight.name = "RNA.weight")

lympho <- alldata

macExp3 <- rbind(lympho@assays$RNA@data,lympho@assays$ADT@data)
macExp3 <- t(macExp3)
metadata2 <- lympho@meta.data
umaps <- cbind(lympho@reductions$wnn.umap@cell.embeddings, lympho@reductions$rna.umap@cell.embeddings, lympho@reductions$adt.umap@cell.embeddings)

macExp4 <- cbind(umaps,metadata2,macExp3)



#MYELO#


myelo <- subset(alldatacp, idents = "Myeloid")

alldata <- myelo

alldata@active.assay <- "RNA"

alldata <- FindMultiModalNeighbors(
     alldata, reduction.list = list("pca", "apca"), 
     dims.list = list(1:47, 1:20), modality.weight.name = "RNA.weight")

for (res in c(0.1, 0.25, 0.5, 1, 1.5, 2)) {
     alldata <- FindClusters(alldata, graph.name = names(alldata@graphs)[[2]], resolution = res, algorithm = 4) # 4 = Leiden ALGORITHM
}

myelo <- alldata





#### SUPERCT PREDICTION ####

library(rSuperCT)

tmp<-rSuperCT::PredCellTypes(as.matrix(alldata@assays$SCT@data),species="mouse")





#### Heatmap for bozena

markers <- read_csv("markers_17.11.21.csv")

listMarkers <- split(markers,markers$cluster)

listMarkers <- split(markers,markers$cluster)


eh <- lapply(listMarkers, function(i){
     return(i[1:5,])
})

test <- do.call("rbind",eh)

genNames <- test[,7]

clustNames <- test[,6]

pop <- macExpAll[,c("annot3",genNames)]
splitList <- split(pop,pop$annot3) 

popList <- lapply(splitList,function(pop){ #we need to make a numeric matrix from char df
     pop <- as.matrix(pop)
     pop <- pop[,-1]
     pop <- matrix(as.numeric(pop),
            ncol = ncol(pop))
})

meanList <- lapply(popList,colMeans)


qwe <- do.call(rbind,meanList)


colnames(qwe) <- genNames



pheatmap(qwe, cluster_cols = FALSE, angle_col = 315)

png("heatmap10TopDifferentiallyExpressedGenesPerCluster.png",height =  520,width =  2400)
print(pheatmap(qwe, cluster_cols = FALSE, cluster_rows = FALSE, angle_col = 315, display_numbers = TRUE, scale="none",  border_color = "white", gaps_col = seq(0, length(0:125), 5)))
dev.off()

write.csv(qwe,"pheatmapV1.csv")


markers %>%
     group_by(cluster) %>%
     top_n(n = 10, wt = avg_log2FC) -> top10

png("heatmap10TopDifferentiallyExpressedGenesPerCluster.png",height =  2000,width =  2600)
DoHeatmap(alldata, features = top10$gene, angle = 90) + NoLegend()
dev.off()




png("lymph_markers.png",width=2500,height = 1500)
FeaturePlot(lymph, reduction = "wnn.umap", features = unique(c("Cd4", "Cd8a", "Cd24a", "Itgax","Ptprc","Cd3d","Nkg7", "H2-Eb2", "Cd44", "Il1r2","Itgax","Ptprc", "Itgam", "Foxp3", "Pdcd1","CD8a-prot","CD4-prot","CD3-prot","NK-1-1-prot","PD1_prot","CD25-prot","CD44-prot","CD62L-prot")), order = T, slot = "data", combine = T, ncol = 6 )
dev.off()



# 14,3 - Double Negative T cells
# 13 - cd 4+ Treg
# 8,10,16 - cd 4+
#      11,6,20,9,12,18,4 - cd 8+
#      1,19 - Bcells
# 21,2,7,5,17 - Nkcells
# 15 - Double Negative T cells?


png("aniaLenkiewicz.png",width = 3600, height = 2000)
FeaturePlot(alldataAfterClust, features= c(testseu(alldataAfterClust,"Nlrp"),testseu(alldataAfterClust,"Nlrc"),
                                           testseu(alldataAfterClust,"Casp"),testseu(alldataAfterClust,"Gsdm"),
                                           "Aim2","Park2","Il1b","Il18","Naip","Mefv","Card8","Pycard","Nos1","Nos2","Irf7"), 
                                           ncol=7, reduction = "wnn.umap")
dev.off()
            
            
            
            
            
            
 
              ################## human TME atlas ##################
## https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/Tumor-Immune-Cell-Atlas ##
library(readr)



list.files(path = ".",  pattern="atlas")

ref <- readRDS("down1000_TMEatlas_expMat.rds")               #"TMEatlas_expMat.csv")



##### handmade integration because they fucked it up

library(patchwork)

ref <- readRDS("TICAatlas_ref.rds")

refList <- SplitObject(ref, split.by = "patient")

refList2 <- mclapply(refList, function(i){
     i@assays$integrated <- NULL
     i <- FindVariableFeatures(i, slection.method = "vst", nfeatures = 3000)
     return(i)
}, mc.cores= 80)

features <- SelectIntegrationFeatures(object.list = refList2)

ref.anchors <- FindIntegrationAnchors(object.list = refList2, anchor.features = features)

# bardzo duza roznica w liczbe komorek od 30 - 8000. Odrzucilbym probki posiadajace mniej niz X
#komorek, pytanie brzmi jaki treshold ustawic -- 300komorek?





#prepare reference dataset to Project TILs program

ref <- FindVariableFeatures(nat, assay = "RNA", nfeatures = 5000)

set.seed(1234)
which.assay="RNA"

varfeat <- ref@assays[[which.assay]]@var.features

refdata <- data.frame(t(ref@assays[[which.assay]]@data[varfeat,]))
refdata <- refdata[, sort(colnames(refdata))]

ref.pca <- prcomp(refdata, rank. = 50, scale. = TRUE, center = TRUE, retx=TRUE)

ref.pca$rotation[1:5,1:5]

library(umap)

seed=1234
n.neighbors=30
min.dist=0.3
metric="cosine"
ndim=10

umap.config <- umap.defaults
umap.config$n_neighbors = n.neighbors
umap.config$min_dist = min.dist
umap.config$metric = metric
umap.config$n_components = 2
umap.config$random_state = seed
umap.config$transform_state = seed

ref.umap <- umap::umap(ref.pca$x[,1:ndim], config = umap.config)

colnames(ref.umap$layout) <- c("umap_1","umap_2")
ref.umap

ref@reductions$umap@cell.embeddings <- ref.umap$layout
ref@reductions$pca@cell.embeddings <- ref.pca$x
ref@reductions$pca@feature.loadings <- ref.pca$rotation
colnames(ref@reductions$pca@cell.embeddings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl=TRUE)
colnames(ref@reductions$pca@feature.loadings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl=TRUE)
#Store the complete PCA and UMAP object in @misc
ref@misc$pca_object <- ref.pca
ref@misc$umap_object <- ref.umap
ref@misc$projecTILs="custom_atlas"

ref$functional.cluster <- ref$cell_type

DimPlot(ref, reduction = "umap", pt.size = 0.5, group.by = "functional.cluster") + ggtitle("UMAP by sample")


saveRDS(ref,"TILsRefTMEatlas.rds")






               ################## Project TILs #################
          # https://github.com/carmonalab/ProjecTILs#reference-atlases #

#install.packages("remotes")
#library(remotes)

#remotes::install_github("carmonalab/UCell")
#remotes::install_github("carmonalab/scGate")
#remotes::install_github("carmonalab/ProjecTILs")


library(ProjecTILs)

ref <- readRDS("natCITEseqTILsAtlas.rds")

#See expression of important marker genes across reference subtypes

markers <- c("Cd4","Cd8a","Ccr7","Tcf7","Pdcd1","Havcr2","Tox","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
VlnPlot(ref,features=markers,stack = T,flip = T,assay = "RNA")


lymph <- readRDS("myelo.rds")
library(nichenetr)

macExp <- lymph@assays$RNA@data
humanNames <- convert_mouse_to_human_symbols(rownames(macExp))
rownames(macExp) <- humanNames
macExp <- macExp[-which(is.na(rownames(macExp))),]

macExp <- sumNotUniqueGenes(macExp)


humanLymph <- CreateAssayObject(macExp)
lymph[["human"]] <- humanLymph
lymph@assays$RNA <- NULL
lymph@assays$ADT <- NULL
lymph@active.assay <- "human"
lymph <- DietSeurat(lymph)
lymph <- FindVariableFeatures(lymph, assay = "human", nfeatures = 5000)


VlnPlot(lymph,features=markers,stack = T,flip = T,assay = "human")

query.projected <- make.projection(myelo, ref=ref, filter.cells=FALSE, query.assay="RNA", skip.normalize=TRUE, ncores=1, future.maxSize=10000)

# querydata <- ProjecTILs::query_example_seurat




#query.projected <- make.projection(test, ref=ref)

plot.projection(ref, query.projected)

query.projected <- cellstate.predict(ref=ref, query=query.projected)
table(query.projected$functional.cluster)

plot.statepred.composition(ref, query.projected,metric = "Percent")

plot.states.radar(ref, query=query.projected, min.cells=30)


af <- query.projected


TcellWith7 <- lymph$wsnn_res.0.1[lymph$wsnn_res.0.1 %in% c(1,4,5,6,7)]
TcellWithOut7 <- lymph$wsnn_res.0.1[lymph$wsnn_res.0.1 %in% c(1,4,5,6)]
           


FeaturePlot(af,features="functional.cluster.conf") + # czemu okregi a nie koła?
     labs(title = "functional.cluster.conf", color = "log2(RPKM+1)") +
     scale_colour_gradientn(colours = (viridis(50))) + # rev(inferno(50))
     # breaks=c(min,max))+
     # , labels=c("min", "max"))+
     xlab("wnnUMAP_1") +
     ylab("wnnUMAP_2") +
     theme_feature



clustNames <- names(table(lymph$TILs))

greyClusters <- function(seuObj, idents, redu="wnn.umap") {
     cellsInCluster <- WhichCells(seuObj, idents = idents)
     DimPlot(seuObj,
             reduction=redu,
             label = F, group.by = "ident",
             cells.highlight = list(cellsInCluster), cols.highlight = c("darkred"),
             cols = "grey", pt.size = 0.3
     ) + ggtitle(idents)
}


output <- lapply(clustNames,function(i){
     greyClusters(lymph, i)
})






png("TILsREDandGREY.png", height = 1200, width = 1500)
grid.arrange(grobs=output)
dev.off()


png("lymphAnnot3Nat.png",width=1400,height = 900)
DimPlot(lymph, label = T, reduction="wnn.umap")
dev.off()



#naive like: 14, 16 , 4, 6 -- ale to chyba nie sa naive like a double negative
NL <- c("Tcf7","Ccr7","Pdcd1","Tnfrsf9") #low Pdcd1 Tnfrsf9
#effector-memory cd8+ 5
EF <- c("Tcf7","Gzmk","Pdcd1") #low to intermediate Pdcd1
#Terminally-exhusted 10
TEx <- c("Gzmk","Pdcd1", "Ctla4", "Lag3", "Tigit", "Havcr2","Tox")
#Pre Exhausted 16
TPEx <- c("Tcf7", "Pdcd1", "Ctla4", "Tox" ,"Havcr2","Gzmk") #low Havcr Gzmk
#CD4+ Th1-like 8
Th1L <- c("Ifng","Cxcr3")   #Cxcr3 tak no nie do końca xD
# CD4+ Tfh (CD4+ follicular-helper)
TfH <- c( "Cxcr5", "Tox", "Slamf6") + ("CXCR5","PDCD1","ICOS","BCL6") # CXCR5 jest kluczowy i go nie ma nigdzie 
#Treg 11 
Treg <- c("Foxp3")
#Th17/helper 13
Th17h<-as.character(expression(fpIL17A, IL17F, BATF, IL2RA, DUSP4, IL21R,
                               CTLA4, IRF4, CCL20, IL26, BATF3, IL4R,
                               STAT3, RORA))
#Cd4+ cytotoxic 9
#NKT 4
NKT <- c("Klrb1c","Cd3e")


tmp <- lymph$wsnn_res.1.5
tmp <- as.character(tmp)
tmp[which(tmp == 13)] <- "Th17"
tmp[which(tmp %in% c(1,17))] <- "Bcells"
tmp[which(tmp %in% c(8))] <- "Th1" # a moze exhausted Th1? https://bmcimmunol.biomedcentral.com/articles/10.1186/s12865-019-0309-9 
#high Pdcd1 (PD-1) i Lag3, wtf do tego sa Nkg7 high a cytotoxic nie XD
tmp[which(tmp %in% c(9))] <- "Cytotoxic CD4+" # Powinny być Gzmb+,Prf1+,Fasl+
tmp[which(tmp %in% c(14))] <- "CD4+ PD-1-" # Cxcr3, Il2ra- too
tmp[which(tmp %in% c(11))] <- "Treg"
tmp[which(tmp %in% c(10))] <- "Terminally-exhausted CD8+"
tmp[which(tmp %in% c(7))] <- "Pre-exhausted CD8+"
tmp[which(tmp %in% c(5))] <- "Effector-memory CD8+"
tmp[which(tmp %in% c(6,16))] <- "Naive CD8+"
tmp[which(tmp %in% c(4))] <- "NKT" 
tmp[which(tmp %in% c(2,3,15))] <- "NK"
tmp[which(tmp %in% c(12))] <- "DN Tcells" # Il2ra+
lymph$myLabels <- tmp
Idents(lymph) <- lymph$myLabels

mygenes <- unique(c(NL,EF,TEx,TPEx,Th1L,TfH,Treg))

eh <- function(alldata,genes=genes,ncol=3){FeaturePlot(alldata, reduction = "wnn.umap", features = genes, order = T, slot = "data", combine = T, ncol = ncol)}
eh(hl,CD4mSC, 5)

#Follicular Helper T Cell Markers from RnD systems
Tfh_2 <- c("Btla","Cd40", "Cd57", "B3gat1","Cd84","Cxcr4","Cxcr5","Icos",
          "Il6ra","Il21r","Neprilysin","Cd10","Ox40","Pd-1","Slam1","Cd150")


FeaturePlot(lymph,Tfh_2,ncol=4,reduction="wnn.umap")
FeaturePlot(ref,Tfh_2,ncol=4,reduction="umap")


#Th1 Cell Markers from RnD systems
Th1_2 <- c("Ccr1","Ccr5","Cxcr3","IFN-gamma","Cd119","Il-12","Il-18","Il-27",
           "Stat1","Stat4","T-bet","Tnfaip1","Il2")

FeaturePlot(lymph,Th1_2,ncol=4,reduction="wnn.umap")
FeaturePlot(ref,Th1_2,ncol=4,reduction="umap")

# jutro -- sprawdz markery z atlasu na ktory sie wyjebales

saveRDS(hl,"hl.rds")

#B cells (general) 
Bc<-as.character(expression(MS4A1, CD79A, CD79B, VPREB3, CD19))
#B cells naive 
BcN<-as.character(expression(IGHM, IGHD, CCR7, SELL))
#B cells activated 
BcA<-as.character(expression(CD69, CD83, JUN, IRF4, MYC))
#B cells memory 
BcM<-as.character(expression(CD27, XBP1, CLEC2D, SSPN))
#B cells unswitched memory 
BcUM<-as.character(expression(IGHM, IGHD, CD27, FCRL5))
#B cells proliferative 
BcP<-as.character(expression(TCL1A, CD79B, CD79A, MS4A1, MKI67))
#Plasma B cells 
PBc<-as.character(expression(JCHAIN, IGHGP, IGKC, IGHM, IGHG3,
                        IGKV1D-39, IGHG2, IGLC3, IGHG1,
                        IGHG4, IGHV4-4, XBP1, CD79A, CD38,
                        MZB1, SSR4, DERL3, FKBP11, PRDX4))
#Plasma blast 
PB<-as.character(expression(IGHGP, IGHG4, JCHAIN, IGLC3, IGHG1,
                        MZB1, IGKC, IGHM, XBP1))
#CD4 memory stem cells 
CD4mSC<-as.character(expression(IL7R, CCR7, TCF7, KLRB1, GPR183,
                        SELL, LEF1, CXCR4, RORA, FOXP1,
                        CCL5, CD3E, CXCR6))
#CD4 effector memory RA 
CD4em<-as.character(expression(TXNIP, KLRD1, KLF2, PRF1, TRBC2,
                        TRBC1, TGFBI, SPON2, CCL4, IL7R,
                        TCF7, SPN, ITGAL, IL6ST, ZNF683))
#CD4 activated 
CD4a<-as.character(expression(DUSP2, CD69, JUN, JUNB))
#CD4 transitional memory 
CD4tm<-as.character(expression(CXCL13, TNFRSF4, TIGIT, IL6ST, PASK,
                        KLRB1, CD40LG, TOX, LEF1, ICOS, CD28,
                        TOX2, CCR7, CD247, RORA, PDCD1))
#CD4 resident effector memory 
CD4rem<-as.character(expression(KLRB1, IL7R, LTB, S100A4, RORA,
                        CD40LG, CXCR6, CD69, IL32))

#Check top variable genes dla poszczegolnych clustrow i przyrownaj je do genow dla poszczegolnych populacji

#CD8 activated 
CD8a<-as.character(expression(IFI6, IFIT3, IFI44L, LY6E, IFI44, STAT1,
                        IFI35, IRF7, IFITM1, IFI16, CD38, IFNG,
                        GZMA, GZMK, CXCL13, CD69, RUNX3,
                        IL32, IL2RG))
#CD8 cytotoxic 
CD8c<-as.character(expression(GZMK, DUSP2, CCL5, ANXA1, NKG7,
                        FYN, XCL2, GZMM, EOMES, CCL4,
                        GZMA, CXCR3, KLRG1, XCL1))
#CD8 exhausted 
CD8e<-as.character(expression(CXCL13, LAG3, HAVCR2, TIGIT, PDCD1,
                        PRF1, TOX))
#CD8 IFN activated 
CD8Ia<-as.character(expression(IFITM1, IFITM3, NFKBIA, CD69, IFITM2,
                        FOS, IL32))
#T cells regulatory 
Tcr<-as.character(expression(TNFRSF4, FOXP3, TNFRSF18, BATF, 
                        IL2RA, IL32, TIGIT, LTB, CTLA4, ICOS))
#T cells proliferative 
Tcp<-as.character(expression(HMGB2, MKI67, HMGN2, TOP2A, CXCL13))
#Th17/helper 
Th17h<-as.character(expression(IL17A, IL17F, BATF, IL2RA, DUSP4, IL21R,
                        CTLA4, IRF4, CCL20, IL26, BATF3, IL4R,
                        STAT3, RORA))
#CXCR5 together with other markers such as PD-1, ICOS, and Bcl-6
eh(hl,toupper(TPEx), 6)

# Czy cluster 4 to nie są NKT?

genes=c("Ifng","Ccr4","Ccr6","Cxcr3","Il4","Il9r","Il3ra","Il2","Il10","Adgrg1",
        "Gzmb","Prf1","Fasl")



eh(lymph,genes,5)


lymph@assays$ADT <- NULL
lymph@reductions$apca@assay.used <- "RNA"
lymph@reductions$adt.umap@assay.used <- "RNA"
lymph@reductions$wnn.umap@assay.used <- "RNA"
lymphCD4 <- subset(lymph, idents=c("CD4+ PD-1-","Cytotoxic CD4+","Th1")) #the lymph/myelo objects are bugged somehow, there is problem with ADT assays when we try to subset so thats why we had to delete it
Idents(lymphCD4) <- lymphCD4$myLabels
markersCD4 <- FindAllMarkers(lymphCD4, only.pos = FALSE, min.pct = 0.5, logfc.threshold=1.0)
listCD4[[1]] <-  markersCD4[markersCD4$cluster == "CD4+ PD-1-",]
listCD4[[1]] <- listCD4[[1]][match(sort(listCD4[[1]]$avg_log2FC),listCD4[[1]]$avg_log2FC),]

tmp <- enrichPathway(gen$entrezgene_id[match(listCD4[[1]]$gene, gen$external_gene_name)], organism="mouse")
mybarplot(tmp, showCategory = 20) + ggtitle(paste("control vs ", "CD4+ PD-1-", sep = ""))


lymphCD8 <- subset(lymph, idents=c("Effector-memory CD8+","Pre-exhausted CD8+","Terminally-exhausted CD8+")) 
Idents(lymphCD8) <- lymphCD8$myLabels
markersCD8 <- FindAllMarkers(lymphCD8, only.pos = FALSE, min.pct = 0.5, logfc.threshold=1.0)
listCD8[[1]] <-  markersCD8[markersCD8$cluster == "Terminally-exhausted CD8+",]
listCD8[[1]] <- listCD8[[1]][match(sort(listCD8[[1]]$avg_log2FC),listCD8[[1]]$avg_log2FC),]

tmp <- enrichPathway(gen$entrezgene_id[match(listCD8[[3]]$gene, gen$external_gene_name)], organism="mouse")
mybarplot(tmp <- enrichPathway(gen$entrezgene_id[match(listCD8[[3]]$gene, gen$external_gene_name)], organism="mouse")
, showCategory = 20) + ggtitle(paste("control vs ", "Terminally Exhausted CD8+", sep = ""))



lymphNK <- subset(lymph, idents=c("2","3"))
Idents(lymphNK) <- lymphNK$wsnn_res.1
markersNK <- FindAllMarkers(lymphNK, only.pos = FALSE, min.pct = 0.5, logfc.threshold=0.75)
listNK <- list()
listNK[[2]] <-  markersNK[markersNK$cluster == "3",]
listNK[[2]] <- listNK[[2]][match(sort(listNK[[2]]$avg_log2FC),listNK[[2]]$avg_log2FC),]

tmp <- enrichPathway(gen$entrezgene_id[match(listNK[[1]]$gene, gen$external_gene_name)], organism="mouse")
mybarplot(tmp, showCategory = 20) + ggtitle(paste("control vs ", "NK_clusterNr3", sep = ""))



# Napisz funkcje która będzie liczyła któremu clustrowi najbliżej do danem zestawu genów markerowych
# lepiej zrobic to tak niz patrzec na oko i jeszcze guilty pleasure dopasowywac wyniki do tego co sie chce uzyskac, trzeba testowac liczenie po medianie i sredniej. Do tego DoHeatmap i RidgePloty trzeba testowac tez

FeaturePlot(lymphCD4,"Nkg7",slot="data") + FeaturePlot(lymphCD4,"Nkg7",slot="scale.data")


myfiles <- list.files(pattern = "*.csv")
myfiles <- myfiles[c(2,3,14)]
for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}




RidgePlot(lymphCD4,c("Pdcd1","Lag3","Tnf","Il2","Ifng"),same.y.lims=TRUE, ncol=3)




nat <- readRDS("exportList_integrated_RNA_ADT_20210409.RDS")



pan <- read_csv("panel_99.csv")
ref <- read_csv("cell_types_1_1.csv")







ggplot(pan, aes("UMAP_1", "UMAP_2")) + 
     geom_point(aes(colour = "Cell_annotation3"))

ggplot(pan, aes(x = UMAP_1, y = UMAP_2)) +
     geom_jitter(size = 0.5, alpha = 0.5, aes_string(color = as.name("Cell_annotation3")))





tmp <- nat[[1]]

library(biomaRt)
ensembl <- useDataset("mmusculus_gene_ensembl",mart=useMart("ensembl"))
anno<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),
            filters = "ensembl_gene_id", values = rownames(tmp) , mart = ensembl, uniqueRows = TRUE)
tmp <- tmp[rownames(tmp) %in% anno$ensembl_gene_id,]
rownames(tmp) <- anno$external_gene_name[match(rownames(tmp),anno$ensembl_gene_id)]
tmp <- tmp[!duplicated(rownames(tmp)),]
nats <- CreateSeuratObject(tmp)
nats[["ADT"]] <- CreateAssayObject(nat[[3]])

nats[["umap"]] <- dim.reduc
nats$labels <- pan$Cell_annotation3

nat <- nats






library(rSuperCT)

tmp <- PredCellTypes(as.matrix(nat@assays$RNA@data), species = 'mouse', model = 'generic_38celltypes', results.dir = './models')

saveRDS(nat,"natCITEseq.rds")

library(scMCA)

tmp <- scMCA(as.matrix(alldata@assays$SCT@data))





mylist <- lapply(names(table(myelo$nat))[-19], function(i){
     return(greyClusters(myelo,i))
})

do.call(grid.arrange,c(mylist,ncol=6))

table(myelo$nat,Idents(myelo))  %>%
     prop.table(margin = 2) %>%
     round(digits = 4) * 100 %>%
     as.matrix %>% heatmap





#### DOUBLET FINDER ####
library(DoubletFinder)
sweep.res <- paramSweep_v3(myelo) 
sweep.stats <- summarizeSweep(sweep.res,GT = FALSE) 
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

# define the expected number of doublet cellscells.
nExp <- round(ncol(myelo) * 0.039)  # we expect 3.9% doublets (10X doublets prop)
myelo <- doubletFinder_v3(myelo, pN = 0.25, pK = pK_choose, nExp = nExp, PCs = 1:10)




######### MYELO #########
library(mart)
library(org.Mm.eg.db)

httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- NULL
while(is.null(mart)){
     try(mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"))
}

gen <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description"), values = rownames(myelo), filters = "external_gene_name", mart = mart, useCache = FALSE)

markers <- FindAllMarkers(homMG, only.pos = FALSE, min.pct = 0.5)

tmp <- enrichPathway(gen$entrezgene_id[match(markers$gene[markers$cluster==2], gen$external_gene_name)], organism="mouse")

mybarplot(tmp, showCategory = 20) + ggtitle(paste("all vs ", "2", sep = ""))


kk1 <- enrichGO(
     gene = na.omit(gen$entrezgene_id[match(markers$gene[markers$cluster==2], gen$external_gene_name)]), OrgDb = "org.Mm.eg.db",
     pvalueCutoff = 0.05, pAdjustMethod = "fdr", ont = "MF", readable = TRUE
)

barplot(kk1, showCategory = 30) + ggtitle(paste("control vs ", "2", sep = ""))




myelo@assays$ADT <- NULL
myelo@reductions$apca@assay.used <- "RNA"
myelo@reductions$adt.umap@assay.used <- "RNA"
myelo@reductions$wnn.umap@assay.used <- "RNA"
homMG <- subset(myelo,idents=c(1,2,8,13))



markers <- FindAllMarkers(homMG, only.pos = TRUE, min.pct = 0.5)


for (i in 13){
tmp <- enrichPathway(gen$entrezgene_id[match(markers$gene[markers$cluster==i], gen$external_gene_name)], organism="mouse")

print(mybarplot(tmp, showCategory = 20) + ggtitle(paste("all vs ",as.character(i), sep = "")))


kk1 <- enrichGO(
     gene = na.omit(gen$entrezgene_id[match(markers$gene[markers$cluster==i], gen$external_gene_name)]), OrgDb = "org.Mm.eg.db",
     pvalueCutoff = 0.05, pAdjustMethod = "fdr", ont = "BP", readable = TRUE
)

print(barplot(kk1, showCategory = 60) + ggtitle(paste("control vs ", as.character(i), sep = "")))
}


genes=as.character(expression(Bhlhe41, Ncoa3, Notch2))
for (i in genes){
mean(homMG@assays$RNA@data[i,homMG$wsnn_res.1 == 1]) %>% print
mean(homMG@assays$RNA@data[i,homMG$wsnn_res.1 == 2]) %>% print
mean(homMG@assays$RNA@data[i,homMG$wsnn_res.1 == 8]) %>% print
}


#Czy cluster 13 to sa rozwalone komorki?
#perc mito perc ribo  hist number of reads i brak markerow typu cd14 i tmem

listNK <- list()
n=0
for (i in c(1,2,8,13)){
     n=n+1
listNK[[n]] <-  markers[markers$cluster == i,]
listNK[[n]] <- listNK[[n]][match(sort(listNK[[n]]$avg_log2FC, decreasing = TRUE),listNK[[n]]$avg_log2FC),]
}

mylist <- list()
for (i in 1:4){
     mylist[[i]] <- list()
     for (j in 1:20){
          mylist[[i]][[j]] <- fp(myelo,listNK[[i]]$gene[[j]])
     }
}


png("clust13Umap.png",width = 2500, height = 2000)
CombinePlots(mylist[[4]][1:20])
dev.off()


markers21vs23 <- FindMarkers(myelo,ident.1=23,ident.2=21, min.pct=0.5)

markers21vs23 <- markers21vs23[match(sort(markers21vs23$avg_log2FC, decreasing = TRUE),markers21vs23$avg_log2FC),]

mylist <- list()
     for (i in rownames(markers21vs23)[1:20]){
          mylist[[i]] <- fp(myelo,i)
     }

png("clust23vs21Umap.png",width = 2500, height = 2000)
CombinePlots(mylist[1:20])
dev.off()


markesr23vsAll <- FindMarkers(myelo,ident.1=23,ident.2=c(1:22,24,25), min.pct=0.5)
markesr23vsAll <- markesr23vsAll[match(sort(markesr23vsAll$avg_log2FC, decreasing = TRUE),markesr23vsAll$avg_log2FC),]
markesr23vsAllR <- markesr23vsAll[match(sort(markesr23vsAll$avg_log2FC, decreasing = FALSE),markesr23vsAll$avg_log2FC),]


mylist <- list()
for (i in rownames(markers21vs23)[1:20]){
     mylist[[i]] <- fp(myelo,i)
}


bc <- read_tsv("BreastCancerMatrix.tsv")


pawel <- FindAllMarkers(myelo,only.pos=TRUE,min.pct=0.5)

for (i in c(unique(myelo$wsnn_res.2))){
     n=n+1
     mylist[[n]] <-  pawel[pawel$cluster == i,]
     mylist[[n]] <- mylist[[n]][match(sort(mylist[[n]]$avg_log2FC, decreasing = TRUE),mylist[[n]]$avg_log2FC),]
}

saveRDS(mylist,"allClustersMarkersList.rds")




seurat2cloupe(lymph, reduction="wnn.umap",dims=1:2,metadata="wsnn_res.1",keyword="lymph_res1",opdir=".")


lymphMarkers <- FindAllMarkers(lymph,min.pct = 0.5, logfc.threshold = 0.75)






tmp <- match(nk4tab$gene,lymphMarkers[lymphMarkers$cluster==14,]$gene)

for (i in length(tmp)){
     if (!(tmp[[i]] %>% is.na())){
          if(nk4tab$p_val_adj[[i]] > 0.05){
               nk4tab$p_val_adj[[i]] <- lymphMarkers$p_val_adj[[tmp[[i]]]]
          }
}}



nk4tab$gene[(nk4tab$p_val_adj > 0.05)]



mark14 <- FindMarkers(myelo,ident.1=14,ident.2=c(9,8,11,10,12,4,5,7,6), min.pct=0.5)
mark9  <-FindMarkers(myelo,ident.1=9,ident.2=c(14,8,11,10,12,4,5,7,6), min.pct=0.5)






#### TUMOR CELLS ON UMAP VS SHAM CELLS ON UMAP
DimPlot(alldata, cells.highlight = WhichCells(alldata, idents=c("tumor1","tumor2","tumor3","tumor4","tumor5","tumor6"))) + DimPlot(alldata, cells.highlight = WhichCells(alldata, idents=c("sham1","sham2","sham3","sham4","sham5","sham6")))







################ SCT IDENTIFICATION ##################

#lymph res: 2.5
# CD4+:
# 22 - effector Th1
# 25 - dysfunctional
# 31 - Tregs
# 38 - immunoadolescent Cd28- 
# CD8+:
# 28 - Effector CD8
# 23 - preexhausted cd8
# 14- Effector memory
# 15- Naive like/central memory
# DN:
# 11 - NKT
# 30 -?? Ctla4+ Cd8, Pdcd1 czesciowo+, troche Gzmc
# 35 - Th17



# Macrophages res: 1.5:
#29- nonclassical monocytes
#21- monocytes 
#6 - Mo/Macrophages
#2 - mature macrophages
#9 -?? there is previous proliferating cluster, and they dont have Vegfa expression -- maybe GMB not induced macrophages?
#25-BAM
#

#Microglia res: 1:5
#30 - proliferating
#Microglia res: 0.5:
#9- activated microglia
#13-higher canonical markers
#5-
#1-
#19-?low number of reads microglia

#Dendritic cells res: 1.5:
#14: cDC2
#27: cDC1
#16: Mo/DC
#21: pDCs

#NK cells res: 2.5:
#17-young
#19-partly matured, exposed to GBM and modified by it, lower level of maturation, proangiogenesis, Cxcr4+ anable them migrate into the tumor, low level of (prf1,gzma,gzmb), low level of activating receptor Ncr1, but not high level of inhibitory receptors, lower level of Infg
#12-matured, high level of (prf1,gzma,gzmb), high level of activating receptors, high level of Ifng

#Bcells res1.5:
#7

#Others res2.5
#41: neutrophiles(ly6G+)
#43: mastcells (kit+)
#39: oligodendrocytes (olig1+,olig2+)
#42: jakies monocyty?
#29: eosinophiles???? (Cd24+,Ccr3+):
#cd24:
#https://www.atsjournals.org/doi/pdf/10.1165/rcmb.2015-0146OC
#https://www.frontiersin.org/articles/10.3389/fimmu.2020.594540/full
#https://pubmed.ncbi.nlm.nih.gov/7510927/
#ccr3:
#https://fibrogenesis.biomedcentral.com/articles/10.1186/s13069-015-0028-7
#https://journals.asm.org/doi/10.1128/JVI.72.4.3351-3361.1998

#naive like: 14, 16 , 4, 6 -- ale to chyba nie sa naive like a double negative
NL <- c("Tcf7","Ccr7","Pdcd1","Tnfrsf9") #low Pdcd1 Tnfrsf9
#effector-memory cd8+ 5
EF <- c("Tcf7","Gzmk","Pdcd1") #low to intermediate Pdcd1
#Terminally-exhusted 10
TEx <- c("Gzmk","Pdcd1", "Ctla4", "Lag3", "Tigit", "Havcr2","Tox")
#Pre Exhausted 16
TPEx <- c("Tcf7", "Pdcd1", "Ctla4", "Tox" ,"Havcr2","Gzmk")





### Markers cd8+ ###

Idents(alldata) <- alldata$labels


mark1 <- FindMarkers(alldata,"CM Cd8+",c("Effector Cd8+","Effector-memory","Naive-like Cd8+","NKT","Pre-exhausted Cd8+","preNKT"))
mark2 <- FindMarkers(alldata,"Effector Cd8+",c("CM Cd8+","Effector-memory","Naive-like Cd8+","NKT","Pre-exhausted Cd8+","preNKT"))
mark3 <- FindMarkers(alldata,"Effector-memory",c("CM Cd8+","Effector Cd8+","Naive-like Cd8+","NKT","Pre-exhausted Cd8+","preNKT"))
mark4 <- FindMarkers(alldata,"Naive-like Cd8+",c("CM Cd8+","Effector Cd8+","Effector-memory","NKT","Pre-exhausted Cd8+","preNKT"))
mark5 <- FindMarkers(alldata,"NKT",c("CM Cd8+","Effector Cd8+","Effector-memory","Naive-like Cd8+","Pre-exhausted Cd8+","preNKT"))
mark6 <- FindMarkers(alldata,"Pre-exhausted Cd8+",c("CM Cd8+","Effector Cd8+","Effector-memory","Naive-like Cd8+","NKT","preNKT"))
mark7 <- FindMarkers(alldata,"preNKT",c("CM Cd8+","Effector Cd8+","Effector-memory","Naive-like Cd8+","NKT","Pre-exhausted Cd8+"))

markers <- list(mark1,mark2,mark3,mark4,mark5,mark6,mark7)

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][match(sort(markers[[i]]$avg_log2FC, decreasing = TRUE),markers[[i]]$avg_log2FC),]
}

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][c(1:25,(length(markers[[i]][[1]])-25):(length(markers[[i]][[1]]))),]
}

clusters <- c("CM Cd8+","Effector Cd8+","Effector-memory","Naive-like Cd8+","NKT","Pre-exhausted Cd8+","preNKT")
for (i in 1:(length(markers))){
     markers[[i]]$cluster <- rep(clusters[[i]],dim(markers[[i]])[1])
}

for (i in 1:(length(markers))){
     markers[[i]]$gene <- rownames(markers[[i]])
}

markers <- do.call(rbind,markers)

markers <- markers[which(markers$p_val_adj < 0.05),]

write_csv(markers,"markersEffectorCD8.csv")

library("rio")

myfiles <- list.files(pattern = "markersCD8_new.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}


#############################################################





############### Markers CD4 ################

mark25 <- FindMarkers(alldata, 25,c(38,31,22))

mark25 <- mark25[-which(mark25$p_val_adj > 0.05),]

mark25 <- mark25[match(sort(mark25$avg_log2FC, decreasing = TRUE),mark25$avg_log2FC),]

mark25 <- mark25[c(1:25,(length(mark25[[1]])-25):(length(mark25[[1]]))),]

mark25$gene <- rownames(mark25)

write_csv(mark25,"mark25.csv")

library("rio")

myfiles <- list.files(pattern = "mark25.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}


#############################################################





############### Markers Macrophages ################


Idents(alldata) <- alldata$labels

### only.pos=TRUE do sciezek, do CSV only.pos=FALSE ###
mark1 <- FindMarkers(alldata,"Phag Mf",c("protumor Mf","Mystery"),only.pos = FALSE) 
mark1 <- mark1[which(mark1$p_val_adj < 0.05),]
mark2 <- FindMarkers(alldata,"protumor Mf",c("Phag Mf","Mystery"),only.pos = FALSE)
mark2 <- mark2[which(mark2$p_val_adj < 0.05),]
mark3 <- FindMarkers(alldata,"Mystery",c("Phag Mf","protumor Mf"),only.pos = FALSE)
mark3 <- mark3[which(mark3$p_val_adj < 0.05),]

markers <- list(mark1,mark2,mark3)

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][match(sort(markers[[i]]$avg_log2FC, decreasing = TRUE),markers[[i]]$avg_log2FC),]
}

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][c(1:25,(length(markers[[i]][[1]])-25):(length(markers[[i]][[1]]))),]
}

clusters <- c("Phag Mf","protumor Mf","not induced Mf")
for (i in 1:(length(markers))){
     markers[[i]]$cluster <- rep(clusters[[i]],dim(markers[[i]])[1])
}

for (i in 1:(length(markers))){
     markers[[i]]$gene <- rownames(markers[[i]])
     markers[[i]] <- markers[[i]][which(markers[[i]]$p_val_adj < 0.05),]
   # markers[[i]] <- markers[[i]][which(markers[[i]]$avg_log2FC < 0.25),] #only.NEG
}
     
     
markers <- do.call(rbind,markers)



write_csv(markers,"markersMacrophages.csv")

library("rio")

myfiles <- list.files(pattern = "markersMacrophages.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}





#### Pathways #####

markers <- list(markers[[3]])

#enrichGO
sciezki <- mclapply(markers,function(i){
     kk1 <- enrichKEGG(gene         = gen$entrezgene_id[match(i$gene,gen$external_gene_name)],
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
     # kk1 <- enrichGO(
     #      gene = na.omit(gen$entrezgene_id[match(i$gene, gen$external_gene_name)]), OrgDb = "org.Mm.eg.db",
     #      pvalueCutoff = 0.05, pAdjustMethod = "fdr", ont = "MF", readable = TRUE
     # )
     return(kk1)
},mc.cores=length(markers))

for (i in 1:length(sciezki)){
sciezki[[i]]@result <- sciezki[[i]]@result[sort(sciezki[[i]]@result$Count, decreasing = TRUE, index.return=TRUE)[[2]],]
}

myplots <- lapply(1:length(sciezki),function(i){
     dotplot(sciezki[[i]], showCategory=25) + ggtitle(clusters[i])
})

png("MysteryClusterKEGGdownreg.png",width=600,height=1000) #,width=1200,height=2000)
do.call(ggarrange,myplots)
dev.off()



#KEGG
sciezki <- mclapply(markers,function(i){
     kk2 <- enrichKEGG(gene         = gen$entrezgene_id[match(i$gene,gen$external_gene_name)],
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
     return(kk2)
},mc.cores=length(markers))

for (i in 1:length(sciezki)){
     sciezki[[i]]@result <- sciezki[[i]]@result[sort(sciezki[[i]]@result$Count, decreasing = TRUE, index.return=TRUE)[[2]],]
}

myplots <- lapply(1:length(sciezki),function(i){
     dotplot(sciezki[[i]], showCategory=25) + ggtitle(clusters[i])
})

myplots <- lapply(1:length(sciezki),function(i){
     dotplot(sciezki[[i]], showCategory=25) + ggtitle(clusters[i])
})

png("macrophagesSciezkiKEGGSorted.png",width=1200,height=2000)
do.call(ggarrange,myplots)
dev.off()

#############################################################






############### Markers Active Microglia ################


Idents(alldata) <- alldata$labels


### only.pos=TRUE do sciezek, do CSV only.pos=FALSE ###
mark1 <- FindMarkers(alldata,"ActMg1",c(c("ActMg2","Proliferating Mg","BAMs")),only.pos = TRUE)
mark1 <- mark1[which(mark1$p_val_adj < 0.05),]
mark2 <- FindMarkers(alldata,"ActMg2",c("ActMg1","Proliferating Mg","BAMs"),only.pos = TRUE)
mark2 <- mark2[which(mark2$p_val_adj < 0.05),]
mark3 <- FindMarkers(alldata,"Proliferating Mg",c("ActMg1","ActMg2","BAMs"),only.pos = TRUE)
mark3 <- mark3[which(mark3$p_val_adj < 0.05),]
mark4 <- FindMarkers(alldata,"BAMs",c("Proliferating Mg","ActMg1","ActMg2"),only.pos = TRUE)
mark4 <- mark4[which(mark4$p_val_adj < 0.05),]

markers <- list(mark1,mark2,mark3,mark4)

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][match(sort(markers[[i]]$avg_log2FC, decreasing = TRUE),markers[[i]]$avg_log2FC),]
}

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][c(1:25,(length(markers[[i]][[1]])-25):(length(markers[[i]][[1]]))),]
}

clusters <- c("ActMg1","ActMg2","Proliferating Mg","BAMs")
for (i in 1:(length(markers))){
     markers[[i]]$cluster <- rep(clusters[[i]],dim(markers[[i]])[1])
}

for (i in 1:(length(markers))){
     markers[[i]]$gene <- rownames(markers[[i]])
}

markers <- do.call(rbind,markers)

markers <- markers[which(markers$p_val_adj < 0.05),]

write_csv(markers,"markersActMicroglia.csv")

library("rio")

myfiles <- list.files(pattern = "markersActMicroglia.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}





#### Pathways #####

#markers <- list(mark1,mark2,mark3,mark4)

sciezki <- mclapply(markers,function(i){
     kk1 <- enrichGO(
          gene = na.omit(gen$entrezgene_id[match(i$gene, gen$external_gene_name)]), OrgDb = "org.Mm.eg.db",
          pvalueCutoff = 0.05, pAdjustMethod = "fdr", ont = "MF", readable = TRUE
     )
     return(kk1)
},mc.cores=length(markers))

for (i in 1:length(sciezki)){
     sciezki[[i]]@result <- sciezki[[i]]@result[sort(sciezki[[i]]@result$Count, decreasing = TRUE, index.return=TRUE)[[2]],]
}

myplots <- lapply(1:length(sciezki),function(i){
     barplot(sciezki[[i]], showCategory=25) + ggtitle(clusters[i])
})

png("actMgSciezkiSorted.png",width=1200,height=2000)
do.call(ggarrange,myplots)
dev.off()
#############################################################






############  CellChat analysis #############


cellchat <- createCellChat(alldata) 
cellchat@idents <- as.factor(alldata$labels)
names(cellchat@idents) <- NULL
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)


cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")




mat <- cellchat@net$weight[,]
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
     mat2 <- matrix(0, nrow = 36, ncol = ncol(mat), dimnames = dimnames(mat))
     mat2[i, ] <- mat[i, ]
     netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pathways.show <- c("IFN-II")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")





#########################################
##### Liczymy procent komorek mieloidalnych danej probki w danym klastrze


maxcells <- (myelo$orig.ident %>% table())[5:10]


tmp1 <- (myelo$orig.ident[myelo$labels=="Mo"] %>% table)[-(1:4)]
tmp2 <- myelo$orig.ident[myelo$labels=="Mo/Mf"] %>% table
tmp3 <- myelo$orig.ident[myelo$labels=="Phag Mf"] %>% table
tmp4 <- myelo$orig.ident[myelo$labels=="protumor Mf"] %>% table
tmp5 <- (myelo$orig.ident[myelo$labels=="not induced Mf"] %>% table)[-1]

print("Mo")
tmp1 <- tmp1/maxcells
print("Mo/Mf")
tmp2 <- tmp2/maxcells
print("Phag Mf")
tmp3 <- tmp3/maxcells
print("protumor Mf")
tmp4 <- tmp4/maxcells
print("not induced Mf")
tmp5 <- tmp5/maxcells


tmp <- data.frame(tmp1,tmp2,tmp3,tmp4,tmp5)[,c(2,4,6,8,10)]
colnames(tmp) <- c("Mo","Mo/Mf","Phag Mf","protumor Mf","not induced Mf")
rownames(tmp) <- c("tumor1","tumor2","tumor3","tumor4","tumor5","tumor6")
pheatmap(as.matrix(tmp), cluster_cols = FALSE, cluster_rows = FALSE, angle_col = "315", display_numbers=TRUE,fontsize=20, main="Procent komorek mieloidalnych danej probki w danym klastrze")


########################################

myelo <- AddModuleScore(myelo,name = "Dojrzalosc",features=list(c("Spp1","Ctsa","Ctsb","Mmp2","Adgre1","Msr1","Cd63","F4-80-prot","I-A-I-E-prot")))
myelo <- AddModuleScore(myelo,name = "Mlodosc",features=list(c("Ly-6C-prot","Cebpb","Ly6i","Plac8","Irf7")))
##














########### IDENTYFIKACJA BOZENY ############


markers <- FindAllMarkers(alldata,  only.pos = FALSE, min.pct = 0.5, logfc.threshold = 0.25)
write_csv(markers,"markersFULL.csv")
markers <- markers[markers$p_val_adj<0.05,]
#markers <- markers[markers$avg_log2FC>0,]
markers[markers$p_val==0,"p_val"] <- 0.000000001
markers[markers$p_val_adj==0,"p_val_adj"] <- 0.000000001
cpmarkers <- markers

mylist <- as.factor(names(table(alldata$oszukaneLabels)))

marklist <- splitDataFrame(markers,mylist)

markers <- marklist

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][match(sort(markers[[i]]$avg_log2FC, decreasing = TRUE),markers[[i]]$avg_log2FC),]
}

for (i in 1:(length(markers))){
     markers[[i]] <- markers[[i]][c(1:25,(length(markers[[i]][[1]])-25):(length(markers[[i]][[1]]))),]
}


redumarkers <- lapply(markers, function(i){
     return(i[1:10,])
})

markers2 <-do.call(rbind,redumarkers)

write.csv(markers2,"redumarkers.csv")
library("rio")
myfiles <- list.files(pattern = "markers.csv")
for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}

Idents(alldata) <- factor(as.character(Idents(alldata)), levels = sort(as.character(levels(Idents(alldata)))))

png("heatmapBoz.png",height = 3000, width=2000)
DoHeatmap(alldata,markers2$gene)
dev.off()


fplist <- mclapply(markers,function(i){
     return(fp(alldata,i$gene[1:4],ncol=2,combine=TRUE)) 
},mc.cores=length(markers))

# png("test.png",height=5000,width=5000) total fail ggarrange
# do.call(ggarrange,c(fplist,ncol=6))
# dev.off()


library(mart)
library(org.Mm.eg.db)

httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- NULL
while(is.null(mart)){
     try(mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"))
}

gen <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description"), values = rownames(alldata), filters = "external_gene_name", mart = mart, useCache = FALSE)

sciezki <- mclapply(markers,function(i){
kk1 <- enrichGO(
     gene = na.omit(gen$entrezgene_id[match(i$gene, gen$external_gene_name)]), OrgDb = "org.Mm.eg.db",
     pvalueCutoff = 0.05, pAdjustMethod = "fdr", ont = "MF", readable = TRUE
)
return(kk1)
},mc.cores=length(markers))

clusnames <- Idents(alldata) %>% levels


sciezkiplots <- mclapply(1:37, function(i){
return(barplot(sciezki[[i]], showCategory=25) + ggtitle(clusnames[i]))
})

png("sciezkibozena.png",width = 2400, height = 4500)
do.call(ggarrange,c(sciezkiplots))
dev.off()


