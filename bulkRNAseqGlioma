library(ggplot2)
library(ggrepel)
library(ggpubr)
library(plotrix)
library("DESeq2")
library("org.Mm.eg.db")
library("clusterProfiler")
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library("ReactomePA")
library(ggnewscale)
library(ComplexHeatmap)
library(data.table)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(EnhancedVolcano)
library(biomaRt)
library(venn)
library(lsa)

setwd("/home/rstudio/all/Opus/cellLinesBUlkRNAseq/glioma")

httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- NULL
while(is.null(mart)){
     try(mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"))
}


counts <- read.table("feature_counts_all.txt", head=TRUE, row.names = 1)

colnames(counts) <- sub(pattern = "Aligned.sortedByCoord.out.bam", replacement = "",colnames(counts))

exon_length <- counts$Length
counts <- counts[,-c(1,2,3,4,5)]

corBetweenLines <- lapply(seq(1, length(counts[1,]), 2),function(i){
     cor(counts[,i],counts[,i+1],method="spearman")
})
corBetweenLines <- unlist(corBetweenLines)
plot(corBetweenLines)
min(corBetweenLines)


counts <- data.frame(counts)

summedLinesList <- lapply(seq(1, length(counts[1,]), 2),function(i){ #sum runs from two flowcells
     return(counts[,as.numeric(i)] + counts[,as.numeric(i+1)])
})

countsSummed <- do.call(cbind,summedLinesList)

rownames(countsSummed) <- rownames(counts)
colnames(countsSummed) <- colnames(counts)[seq(1,length(colnames(counts)),2)]
colnames(countsSummed) <- gsub('.{5}$','',colnames(countsSummed))


seqname <- read.csv("namesFromSequencing.csv", header = FALSE)


eh <- colnames(countsSummed)

eh <- lapply(eh,function(i){
     paste(strsplit(i,"_")[[1]][1:length(strsplit(i,"_")[[1]])-1],collapse ="_")
})

eh <- unlist(eh)
seqname <- seqname[match(eh,seqname$V1),]

colnames(countsSummed) <- eh

summedLinesList2 <- lapply(seq(1, length(countsSummed[1,]), 2),function(i){ #sum runs from two separate illumina runs
     return(countsSummed[,as.numeric(i)] + countsSummed[,as.numeric(i+1)])
     
})

countsSummed2 <- do.call(cbind,summedLinesList2)

meh <- lapply(seqname$V2,function(i){
     paste(strsplit(i,"_")[[1]][1:(length(strsplit(i,"_")[[1]])-2)],collapse ="_")
})

colnames(countsSummed2) <- meh[seq(1, length(meh), 2)]



colnamesMod <- read.csv("colnamesModified.csv", header=FALSE)
colnames(countsSummed2) <- colnamesMod[,1]
RNA_all <- countsSummed2
remove(countsSummed2)

RNA_all_matrix <- matrix(as.numeric(unlist(RNA_all)),nrow=nrow(RNA_all))
rownames(RNA_all_matrix) <- rownames(counts)
colnames(RNA_all_matrix) <- colnames(RNA_all)
RNA_all_matrix <- RNA_all_matrix[rowSums(RNA_all_matrix) >= 90,]

gen <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description"), values = rownames(RNA_all), filters = "ensembl_gene_id", mart = mart, useCache = FALSE)


condition <- factor(read.csv("colnames.csv", header=TRUE)[,2])
coldata <- data.frame(row.names=colnames(RNA_all_matrix), condition)
dds <- DESeqDataSetFromMatrix(countData=RNA_all_matrix, colData=coldata, design= ~condition)
dds<- DESeq(dds)


dinoraMat <- RNA_all_matrix[,1:16]
dinoraMat <- dinoraMat[,-5] #odrzucamy outliera
dinoraCol <- as.data.frame(coldata[c(1:4,6:16),])
colnames(dinoraCol) <- "condition"
rownames(dinoraCol) <- dinoraCol$condition
dinoraCol$condition <- factor(dinoraCol$condition)
condition <- factor(c(rep("NF53",4),rep("O_NF53_wt",3),rep("O_NF53_RAE",4),rep("O_NF53_RAH",4)), levels = c("NF53","O_NF53_wt","O_NF53_RAE","O_NF53_RAH"))
dinoraCol$condition <- condition
ddsDinora <- DESeqDataSetFromMatrix(countData=dinoraMat, colData=dinoraCol, design= ~condition)
ddsDinora <- DESeq(ddsDinora)
vstDinora <- DESeq2::vst(ddsDinora)
p <- plotPCA(vstDinora, intgroup = "condition")
rld <- rlog(ddsDinora, blind=T)
rld_mat <- assay(vstDinora)
rld_cor <- cor(rld_mat)
head(rld_cor)
pheatmap(rld_cor, cluster_rows = FALSE, cluster_cols = FALSE, angle = "315")


castroMat <- RNA_all_matrix[,17:28]
castroCol <- as.data.frame(coldata[17:28,])
colnames(castroCol) <- "condition"
rownames(castroCol) <- castroCol$condition
castroCol$condition <- factor(castroCol$condition)
condition <- factor(c(rep("NPA_C54B",3),rep("NPAI_C3",3),rep("Arf_wt",3),rep("Arf_mt",3)), levels = c("NPA_C54B","NPAI_C3","Arf_wt","Arf_mt"))
castroCol$condition <- condition
ddsCastro <- DESeqDataSetFromMatrix(countData=castroMat, colData=castroCol, design= ~condition)
ddsCastro <- DESeq(ddsCastro)
vst <- DESeq2::vst(ddsCastro)
p <- plotPCA(vst, intgroup = "condition")
p
rld <- rlog(ddsCastro, blind=T)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
head(rld_cor)
pheatmap(rld_cor, cluster_rows = FALSE, cluster_cols = FALSE, angle = "315")



     
## regularized log transformation for clustering, transform counts into log2 counts in order to normalize with respect to rows and library size
vst <- vst(dds, blind=FALSE)
# p <- plotPCA(vst, intgroup = "condition", return = TRUE) to work with ggplot2
p <- plotPCA(vst, intgroup = "condition")

rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
head(rld_cor)
pheatmap(rld_cor, cluster_rows = FALSE, cluster_cols = FALSE, angle = "315")



############# DINORA ############

results <- results(ddsDinora, contrast=c("condition", "O_NF53_wt", "NF53"), independentFiltering = TRUE)
m1<-match(results@rownames,gen[,1])
results1<-cbind(hugo = gen[m1,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(results1) <- results@rownames

results <- results(ddsDinora, contrast=c("condition", "O_NF53_RAE", "O_NF53_wt"), independentFiltering = TRUE)
m2<-match(results@rownames,gen[,1])
results2<-cbind(hugo = gen[m2,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(results2) <- results@rownames   
  
results <- results(ddsDinora, contrast=c("condition", "O_NF53_RAH", "O_NF53_wt"), independentFiltering = TRUE)
m3<-match(results@rownames,gen[,1])
results3<-cbind(hugo = gen[m3,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(results3) <- results@rownames

results <- results(ddsDinora, contrast=c("condition", "O_NF53_RAE", "NF53"), independentFiltering = TRUE)
m4<-match(results@rownames,gen[,1])
results4<-cbind(hugo = gen[m4,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(results4) <- results@rownames

results <- results(ddsDinora, contrast=c("condition", "O_NF53_RAH", "NF53"), independentFiltering = TRUE)
m5<-match(results@rownames,gen[,1])
results5<-cbind(hugo = gen[m5,2] , log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(results5) <- results@rownames

results1 <- data.frame(results1)
results1$log2FoldChange <- as.numeric(results1$log2FoldChange)
results1$pvalue <- as.numeric(results1$pvalue)
results1$padj <- as.numeric(results1$padj)
results2 <- data.frame(results2)
results2$log2FoldChange <- as.numeric(results2$log2FoldChange)
results2$pvalue <- as.numeric(results2$pvalue)
results2$padj <- as.numeric(results2$padj)
results3 <- data.frame(results3)
results3$log2FoldChange <- as.numeric(results3$log2FoldChange)
results3$pvalue <- as.numeric(results3$pvalue)
results3$padj <- as.numeric(results3$padj)
results4 <- data.frame(results4)
results4$log2FoldChange <- as.numeric(results4$log2FoldChange)
results4$pvalue <- as.numeric(results4$pvalue)
results4$padj <- as.numeric(results4$padj)
results5 <- data.frame(results5)
results5$log2FoldChange <- as.numeric(results5$log2FoldChange)
results5$pvalue <- as.numeric(results5$pvalue)
results5$padj <- as.numeric(results5$padj)

res_list <- list(results1,results2,results3,results4,results5)

write.table(results5,"O_NF53_wt_vs_O_NF53_RAH", sep="\t",row.names = TRUE,col.names = NA)

saveRDS(res_list,"res_list.rds")

res_list <- readRDS("res_list.rds")




##################  PLOTS  DINORAH  #####################

n=0
for (i in res_list){
     
     results <- results1
     
     m<-match(rownames(results),gen$ensembl_gene_id)
     
     res<-cbind(rownames(results),gen[m,],results$log2FoldChange,results$pvalue,results$padj)
     
     colnames(res)<-c("gens",colnames(gen),"log2FoldChange","pvalue","padj")
     
     
     w<-which(results$padj<0.05 & (results$log2FoldChange> 2) | results$log2FoldChange< -2)
     
     volc_plot <- EnhancedVolcano(res,lab=res$external_gene_name,x="log2FoldChange",y="padj",
                                  pCutoff=0.05, FCcutoff=1, labSize=5, #boxedLabels = TRUE,
                                  drawConnectors = TRUE, widthConnectors = 0.75,
                                  maxoverlapsConnectors = 10,
                                  selectLab = res$external_gene_name[which((res$log2FoldChange > 2 | res$log2FoldChange < -2) & res$padj < 0.01)],
                                  xlim = c(-12,12))  + 
          theme(legend.position="none") + ggtitle("RAH vs O_NF53")
     volc_plot
     
     wp<-which(results$padj<0.05 & results$log2FoldChange> 2)
     wn<-which(results$padj<0.05 & results$log2FoldChange< -2)
     
     kk1p <- enrichGO(gene         = res$entrezgene_id[wp],OrgDb="org.Mm.eg.db",
                     pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP", readable=TRUE)
     kk1n <- enrichGO(gene         = res$entrezgene_id[wn],OrgDb="org.Mm.eg.db",
                     pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP", readable=TRUE)
     
     GENELIST <- res$log2FoldChange
     names(GENELIST)<-res$entrezgene_id
     
     cnetplot(kk1p, foldChange=GENELIST, showCategory =20) + ggtitle("O_NF53 vs NF53") + theme(axis.text.x =element_text(size=18),
                                                                                               axis.text.y =element_text(size=18),
                                                                                               axis.title=element_text(size=22),
                                                                                               plot.title=element_text(size=25,vjust = 1),
                                                                                               legend.text = element_text(size=13),
                                                                                               legend.title = element_text(size=18))
     cnetplot(kk1n, foldChange=GENELIST, showCategory =20) + ggtitle("O_NF53 vs NF53") + theme(axis.text.x =element_text(size=18),
                                                                                               axis.text.y =element_text(size=18),
                                                                                               axis.title=element_text(size=22),
                                                                                               plot.title=element_text(size=25,vjust = 1),
                                                                                               legend.text = element_text(size=13),
                                                                                               legend.title = element_text(size=18))
     #dotplotReducedONF53vsNF53underexpressed.png
     
     
     ### reduction to choosen pathways ###
     kk1p@result <- kk1p@result[c(1,4,6,8,10,11,12,15,20,23,25),]
     kk1n@result <- kk1n@result[c(1,4,5,6,7,8,10,12,14,18,19,20),]
     ###
     
     kk1p_plot <- dotplot(kk1p, showCategory=25) + ggtitle("RAH vs O_NF53 (overexpressed genes)")
     kk1n_plot <- dotplot(kk1n, showCategory=25) + ggtitle("RAH vs O_NF53 (underexpressed genes)")
     kk1p_plot + theme(axis.text.x =element_text(size=18),
          axis.text.y =element_text(size=18),
          axis.title=element_text(size=22),
          plot.title=element_text(size=25,vjust = 1),
          legend.text = element_text(size=13),
          legend.title = element_text(size=18))
     kk1n_plot + theme(axis.text.x =element_text(size=18),
                       axis.text.y =element_text(size=18),
                       axis.title=element_text(size=22),
                       plot.title=element_text(size=25,vjust = 1),
                       legend.text = element_text(size=13),
                       legend.title = element_text(size=18))

     # print(kk1p_plot)
     # print(kk1n_plot)
     
     

     
     
     ############# VENN DIAGRAM ###############
     sumMat <- lapply(seq(1, 16, 4),function(i){ 
          return(cbind(RNA_all_matrix[,as.numeric(i)] + RNA_all_matrix[,as.numeric(i+1)] + RNA_all_matrix[,as.numeric(i+2)] + RNA_all_matrix[,as.numeric(i+3)]))
     })
     sumMat <- cbind(sumMat[[1]],sumMat[[2]],sumMat[[3]],sumMat[[4]])
     colnames(sumMat) <- c("NF53","O_NF53","RAE","RAH")
     
     #list(rownames(results2)[which(results2$log2FoldChange < -2)], rownames(results3)[which(results3$log2FoldChange < -2)]
     

     venn.diagram(
          x = list(rownames(results1)[which(results1$log2FoldChange < -2)], 
                   rownames(results4)[which(results4$log2FoldChange < -2)], 
                   rownames(results5)[which(results5$log2FoldChange < -2)]),
          category.names = c("Opn-Wt" , "Opn-RAE" , "Opn-RAH"),
          sub.cex = 0.5,
          filename = 'venn.png',
          output = TRUE ,
          imagetype="png" ,
          height = 480 , 
          width = 480 , 
          resolution = 300,
          compression = "lzw",
          lwd = 1,
          col=c("#440154ff", '#009900', '#21908dff'),
          fill = c(alpha("#440154ff",0.3), alpha('#009900',0.3), alpha('#21908dff',0.3)),
          cex = 0.6,
          fontfamily = "sans",
          cat.cex = 0.5,
          cat.default.pos = "outer",
          cat.pos = c(-27, 27, 135),
          cat.dist = c(0.055, 0.055, 0.085),
          cat.fontfamily = "sans",
          cat.col = c("#440154ff", '#009900', '#21908dff'),
          rotation = 1
     )
     
     venn.diagram(
          x = list(rownames(results1)[which(results2$log2FoldChange < -2)], 
                   rownames(results4)[which(results3$log2FoldChange < -2)]),
          category.names = c("Opn-RAE" , "Opn-RAH"),
          sub.cex = 0.5,
          filename = 'venn.png',
          output = TRUE ,
          imagetype="png" ,           
          height = 480 , 
          width = 480 , 
          resolution = 300,
          compression = "lzw",
          lwd = 1,
          col=c("#440154ff", '#21908dff'),
          fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
          cex = 0.6,
          fontfamily = "sans",
          cat.cex = 0.5,
          cat.pos = c(-0, 170),
          main.cex = 1,
          cat.dist = c(0.055, 0.055),
          cat.fontfamily = "sans",
          cat.col = c("#440154ff", '#21908dff'))

     
     
     #optional: to plot names of the intersection; output=FALSE
     grid.newpage()
     grid.draw(v)
     #########################
     
     
     stemgenes <- c("Aldh1a1","Abcg2","L1cam","Aldh3a1","Msi1","Prom1","Olig2","Sox2","Pou5f1","Msi2","Nes","Nanog","Fut4")
     diffgenes <- c("Mbp","Cd44","Map2","Gria2","Gfap","S100b","Tubb3","Gria1","Sema3a","Vim","Tubb")
     EMTgenes <- c("Fn1","Malat1","Mcam","Nr2f1","Rras","Rundc3b")
     EMTgenesBartek <- unlist(h2m(c("VIM","CDH1","FN1","ZEB1","CDH2","MMP2","SNAI2","ZEB2","SPARC","SNAI1","CCN2","TWIST1","CDH11","CLDN4","EPCAM","SERPINE1","TGFB1","COL3A1","ESRP1","INHBA","PMP22","WNT5A","ST14","TNC","CLDN7","EMP3","FSTL1","ITGA5","LAMC2","VCAN","CDH3","COL5A2","DCN","DSP","ERBB3","POSTN","RAB25","SPINT2","TGM2","CDS1","COL1A1","ESRP2","GRHL2","HTRA1","MAP7","MMP9","OCLN","PDGFRB","PRSS8","AP1M2","AP1S2","AXL","CD44","COL1A2","COL5A1","FBN1","FGFBP1","FGFR1","FOXC2","IRF6")))
     autophagygenes <- c("Map1lc3b","Ulk1","Ulk2","Becn1","Atg12","Atg7","Atg5","Atg3","Atg10","Atg16l1","Pik3c3","Atg9a","Atg9b","Map1lc3a","Map1lc3c","Wlpl1","Wlpl2","Pik3r4","Atg13","Atg14","Rb1cc1","Atg101","Atg2a","Atg2b","Atg4a","Atg4b","Atg4c","Atg4d")
     #ECMgenes <- tmp[tmp$Annotated.Term=="extracellular matrix","Symbol"]
     #write.csv(ECMgenes,"ECMgenes.csv")
     ECMgenes <- read_csv("ECMgenes.csv")
     ECMgenes <- unique(ECMgenes$x)
     
     
     HugoDinoraMat <- dinoraMat
     rownames(HugoDinoraMat) <- gen$external_gene_name[match(rownames(dinoraMat),gen$ensembl_gene_id)]
     
     
     stemgenes<- stemgenes[stemgenes %in% rownames(HugoDinoraMat)]
     matS <- HugoDinoraMat[stemgenes,]
     for (i in 1:length(rownames(matS))){
          matS[i,] <- ((matS[i,] - mean(matS[i,])) / sd(matS[i,]))
     }
     #matS <- log2(matS+1)
     pheatmap(matS, cluster_rows = TRUE, cluster_cols = FALSE, angle = "315", main="Stem genes")
     
     
     
     autophagygenes<- autophagygenes[autophagygenes %in% rownames(HugoDinoraMat)]
     matS <- HugoDinoraMat[autophagygenes,]
     for (i in 1:length(rownames(matS))){
          matS[i,] <- ((matS[i,] - mean(matS[i,])) / sd(matS[i,]))
     }
     #matS <- log2(matS+1)
     pheatmap(matS, cluster_rows = TRUE, cluster_cols = FALSE, angle = "315", main="autophagy genes")
     
     
     profgenes<- profgenes[profgenes %in% rownames(HugoDinoraMat)]
     matP <- HugoDinoraMat[profgenes,]
     for (i in 1:length(rownames(matP))){
          matP[i,] <- ((matP[i,] - mean(matP[i,])) / sd(matP[i,]))
     }
     #matP <- log2(matP+1)
     pheatmap(matP, cluster_rows = TRUE, cluster_cols = FALSE, angle = "315", main="Differentation genes")
     
     
     EMTgenes<- EMTgenes[EMTgenes %in% rownames(HugoDinoraMat)]
     matEMT <- HugoDinoraMat[EMTgenes,]
     for (i in 1:length(rownames(matEMT))){
          matEMT[i,] <- ((matEMT[i,] - mean(matEMT[i,])) / sd(matEMT[i,]))
     }
     #matEMT <- log2(matEMT+1)
     pheatmap(matEMT, cluster_rows = TRUE, cluster_cols = FALSE, angle = "315", main="EMT genes")
     
     
     EMTgenes2<- EMTgenesBartek[EMTgenesBartek %in% rownames(HugoDinoraMat)]
     matEMT2 <- HugoDinoraMat[EMTgenes2,]
     for (i in 1:length(rownames(matEMT2))){
          matEMT2[i,] <- ((matEMT2[i,] - mean(matEMT2[i,])) / sd(matEMT2[i,]))
          for (j in 1:dim(matEMT2)[[2]]){
               if (is.na(matEMT2[i,j])){
                    matEMT2[i,j] <- 0
               }
               else if (is.infinite(matEMT2[i,j])){
                    matEMT2[i,j] <- 5
               }
          }
     }
     #matEMT2 <- log2(matEMT2+1)
     pheatmap(matEMT2, cluster_rows = TRUE, cluster_cols = FALSE, angle = "315", main="EMT genes")
     

     ECMgenes<- ECMgenes[ECMgenes %in% rownames(HugoDinoraMat)]
     matECM <- HugoDinoraMat[ECMgenes,]
     for (i in 1:length(rownames(matECM))){
          matECM[i,] <- ((matECM[i,] - mean(matECM[i,])) / sd(matECM[i,]))
          for (j in 1:dim(matECM)[[2]]){
               if (is.na(matECM[i,j])){
                    matECM[i,j] <- 0
               }
               else if (is.infinite(matECM[i,j])){
                    matECM[i,j] <- 5
               }
          }
     }
     #matECM <- log2(matECM+1)
     png("EMTgenes.png",width=800,height = 3000)
     pheatmap(matECM, cluster_rows = TRUE, cluster_cols = FALSE, angle = "315", main="ECM genes")
     dev.off()
     
     
     combo <- c(diffgenes,stemgenes)
     combo<- combo[combo %in% rownames(HugoDinoraMat)]
     matS <- HugoDinoraMat[combo,]
     for (i in 1:length(rownames(matS))){
          matS[i,] <- ((matS[i,] - mean(matS[i,])) / sd(matS[i,]))
     }
     #matS <- log2(matS+1)
     pheatmap(matS, cluster_rows = TRUE, cluster_cols = FALSE, angle = "315", main="Diff and stem genes")
     
     
     
     #
     #w<-which(results$padj<0.05 & (results$log2FoldChange>2 | results$log2FoldChange< -2))
     ###Analiza wzbogaceń w ścieżki sygnałowe KEGG
     kk2 <- enrichKEGG(gene         = res$entrezgene_id[w],
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
     kk2_plot <- mybarplot(kk2, showCategory=20) + ggtitle("O_NF53_wt vs NF53")
     print(kk2_plot)
     
     #w<-which(results$padj<0.05 & (results$log2FoldChange>2 | results$log2FoldChange< -2))
     kk3 <- enrichPathway(gene=res$entrezgene_id[w],pvalueCutoff=0.1, organism="mouse",readable=T)
     kk3_plot <- dotplot(kk3, showCategory=20) + ggtitle("NPA_C54B_vs_NPAI_C3")
     print(kk3_plot)
     
     GENELIST <- res$log2FoldChange
     names(GENELIST)<-res$entrezgene_id
     
     

}

#########################################################




library("rio")
myfiles <- list.files(pattern = "*.tsv")

for (i in myfiles){
convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}





############### FOR ANIA LENKIEWICZ #################

testmacrev(RNA_all_matrix,"Nlrp")


genNames <- c(testmacrev(RNA_all_matrix,"Nlrp"),testmacrev(RNA_all_matrix,"Nlrc"),testmacrev(RNA_all_matrix,"Casp"),testmacrev(RNA_all_matrix,"Gsdm"),
  "Aim2","Park2","Il1b","Il18","Naip","Mefv","Card8","Pycard","Nos1","Nos2","Irf7")

#outliers Casp3, Casp8ap2.Casp2,Irf7
outliers <- c("Casp3","Casp8ap2","Casp2","Irf7")
genNamesRed <- genNames[!(genNames %in% outliers)]

redMat <- RNA_all_matrix[which(rownames(RNA_all_matrix) %in% genNames),]
redMat2 <- RNA_all_matrix[which(rownames(RNA_all_matrix) %in% genNamesRed),]
redMat3 <- RNA_all_matrix[which(rownames(RNA_all_matrix) %in% outliers),]

pheatmap(redMat3, display_numbers = FALSE, angle_col = "315", cluster_rows = FALSE, cluster_cols = FALSE)
