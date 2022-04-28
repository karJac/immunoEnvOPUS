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
library(sva)

setwd("/home/rstudio/all/Opus/cellLinesBUlkRNAseq/microglia")

httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- NULL
while(is.null(mart)){
     try(mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"))
}


counts <- read.table("feature_counts_BAM", head=TRUE, row.names = 1)


colnames(counts) <- sub(pattern = "Aligned.sortedByCoord.out.bam", replacement = "",colnames(counts))

exon_length <- counts$Length
counts <- counts[,-c(1,2,3,4,5)]

corBetweenLines <- lapply(seq(1, length(counts[1,]), 2),function(i){
     cor(counts[,i],counts[,i+1],method="spearman")
})
corBetweenLines <- unlist(corBetweenLines)
mean(corBetweenLines)

min(corBetweenLines)

counts <- data.frame(counts)

tmp <- factor(colnames(counts))

summedLinesList <- lapply(seq(1, length(counts[1,]), 4),function(i){ 
     return(counts[,as.numeric(i)] + counts[,as.numeric(i+1)] + counts[,as.numeric(i+2)] + counts[,as.numeric(i+3)]) #sum sequencing of the same sample from 4 different lines on NovaSeq
})

countsSummed <- do.call(cbind,summedLinesList)

rownames(countsSummed) <- rownames(counts)
colnames(countsSummed) <- colnames(counts)[seq(1,length(colnames(counts)),4)]
colnames(countsSummed) <- gsub('.{9}$','',colnames(countsSummed))

countsSummed <- countsSummed[,-1] #outlier removal

############## BATCH EFFECT REMOVAL ################

batch <- c(rep(1,4),rep(2,5),rep(3,5),rep(4,5))

adjusted <- ComBat_seq(countsSummed, batch=batch, group=NULL)

################################

countsSummed <- adjusted



seqname <- read.csv("colnames.csv", header = FALSE)
seqname <- seqname[-1,]

eh <- colnames(countsSummed)

eh <- lapply(eh,function(i){
     paste(strsplit(i,"_")[[1]][c(1,2,3,4,5)],collapse ="_")
})

eh <- unlist(eh)
seqname <- seqname[match(eh,seqname$V1),]

colnames(countsSummed) <- seqname[,2]

RNA_all <- countsSummed

RNA_all_matrix <- matrix(as.numeric(unlist(RNA_all)),nrow=nrow(RNA_all))
rownames(RNA_all_matrix) <- rownames(counts)
colnames(RNA_all_matrix) <- colnames(RNA_all)
RNA_all_matrix <- RNA_all_matrix[rowSums(RNA_all_matrix) >= 160,]


gen <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description"), values = rownames(RNA_all), filters = "ensembl_gene_id", mart = mart, useCache = FALSE)

rownamesHugo <- gen[match(rownames(RNA_all_matrix),gen[,1]),2]
rownames(RNA_all_matrix) <- rownamesHugo

#condition <- factor(gsub('.{2}$','',colnames(RNA_all_matrix)))
coldata <- readRDS("coldata.rds")
coldata <- coldata[-1,]
rownames(coldata) <- colnames(RNA_all_matrix)
dds <- DESeqDataSetFromMatrix(countData=RNA_all_matrix, colData=coldata, design= ~ batch + condition)
dds<- DESeq(dds)

## regularized log transformation for clustering, transform counts into log2 counts in order to normalize with respect to rows and library size
my_vst <- DESeq2::vst(dds, blind=FALSE)



plotPCA(my_vst, intgroup = "batch")
plotPCA(my_vst, intgroup = "condition")



vst_mat <- assay(my_vst)
vst_cor <- cor(vst_mat)
head(vst_cor)
pheatmap(vst_cor, cluster_rows = FALSE, cluster_cols = TRUE, angle = "315", display_numbers = TRUE)

results <- results(dds, contrast=c("condition", "NF53", "control"), independentFiltering = TRUE, pAdjustMethod="BH")
m1<-match(results@rownames,gen[,1])
results1<-cbind(hugo = gen[m1,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(resultsNF53) <- results@rownames

results <- results(dds, contrast=c("condition", "o-NF53-wild-type", "NF53"), independentFiltering = TRUE)
m2<-match(results@rownames,gen[,1])
results2<-cbind(hugo = gen[m2,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(resultsWT) <- results@rownames   

results <- results(dds, contrast=c("condition", "o-NF53-wild-type", "o-NF53-RAE"), independentFiltering = TRUE)
m3<-match(results@rownames,gen[,1])
results3<-cbind(hugo = gen[m3,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(resultsRAE) <- results@rownames

resultso <- results(dds, contrast=c("condition", "o-NF53-wild-type", "o-NF53-RAH"), independentFiltering = TRUE)
m4<-match(results@rownames,gen[,1])
results4<-cbind(hugo = gen[m4,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(resultsRAH) <- results@rownames


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

res_list <- list(resultsNF53,resultsWT,resultsRAE,resultsoRAH) #(resultsNF53,resultsWT,resultsRAE,resultsoRAH)


saveRDS(res_list,"res_list.rds")

#res_list <- readRDS("res_list.rds")
resNames <- c("NF53","o-NF53-wild-type","o-NF53-RAE","o-NF53-RAH") #


plotsList <- mclapply(1:2, function(i) {
  myName <- resNames[i]
  results <- res_list[[i]]

  
  m<-match(rownames(results),gen$external_gene_name)
  
  res<-cbind(rownames(results),gen[m,],results$log2FoldChange,results$pvalue,results$padj)
  
  colnames(res)<-c("gens",colnames(gen),"log2FoldChange","pvalue","padj")
  
  
  w<-which(results$padj<0.05 & (results$log2FoldChange> 0.25 | results$log2FoldChange< -0.25))
  
  volc_plot <- EnhancedVolcano(res,lab=res$external_gene_name,x="log2FoldChange",y="padj",
                               pCutoff=0.05, FCcutoff=0.5, labSize=3, #boxedLabels = TRUE,
                               drawConnectors = TRUE, widthConnectors = 0.75,
                               maxoverlapsConnectors = 10,
                               selectLab = res$external_gene_name[which((res$log2FoldChange > 0.5 | res$log2FoldChange < -0.5) & res$padj < 0.05)],
                               xlim = c(-12,12))  + 
       theme(legend.position="none") + ggtitle("RAH vs o-NF53")
  volc_plot

  
  wp<-which(results$padj<0.05 & results$log2FoldChange> 0.25)
  wn<-which(results$padj<0.05 & results$log2FoldChange< -0.25)
  
  kk1p <- enrichGO(gene         = res$entrezgene_id[wp],OrgDb="org.Mm.eg.db",
                   pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP", readable=TRUE)
  kk1n <- enrichGO(gene         = res$entrezgene_id[wn],OrgDb="org.Mm.eg.db",
                   pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP", readable=TRUE)
  
  
  kk1p_plot <- barplot(kk1p, showCategory=25) + ggtitle("o-NF53 vs NF53(overexpressed genes)")
  kk1n_plot <- barplot(kk1n, showCategory=25) + ggtitle("o-NF53 vs NF53 (underexpressed genes)")
  print(kk1p_plot)
  print(kk1n_plot)
  
  
  GENELIST <- res$log2FoldChange
  names(GENELIST)<-res$entrezgene_id
  
  cnetplot(kk1p, foldChange=GENELIST) + ggtitle("o-NF53 vs NF53")
  cnetplot(kk1n, foldChange=GENELIST) + ggtitle("o-NF53 vs NF53")
  
  
  
  
  
  

  ### Analiza wzbogaceń w ścieżki sygnałowe KEGG
  kk2 <- enrichKEGG(
    gene = na.omit(gen$entrezgene_id[match(rownames(results), gen$external_gene_name)]),
    organism = "mmu",
    pvalueCutoff = 0.05
  )
  kk2_plot <- mybarplot(kk2, showCategory = 20) + ggtitle(paste("control vs ", myName, sep = ""))

  kk3 <- enrichPathway(gene = na.omit(gen$entrezgene_id[match(rownames(results), gen$external_gene_name)]), pvalueCutoff = 0.1, organism = "mouse", readable = T)
  kk3_plot <- mybarplot(kk3, showCategory = 20) + ggtitle(paste("control vs ", myName, sep = ""))

  GENELIST <- res$log2FoldChange
  names(GENELIST) <- res$entrezgene_id

  cplot <- cnetplot(kk1, foldChange = GENELIST) + ggtitle(paste("control vs ", myName, sep = ""))

  mylist <- list(volc_plot + kk1_plot + kk2_plot + kk3_plot + cplot)
  return(mylist)
  #return(results)
}, mc.cores = length(res_list))



### writing tabels ###
lol <- lapply(plotsList, data.frame)

for (i in 1:3){
     write.csv(lol[i],paste(resNames[i],"_vs_NF53",".csv",sep=""))
}

library("rio")
myfiles <- list.files(pattern = "*.csv")

for (i in myfiles){
     convert(i,paste(gsub('.{4}$','',i),".xlsx",sep=""))
}
##################################



png("cnetPlot.png",height=1400,width=1400)
ggarrange(plotsList[[1]][[1]][[5]],plotsList[[2]][[1]][[5]],plotsList[[3]][[1]][[5]], plotsList[[4]][[1]][[5]], ncol=2, nrow=2)
dev.off()

png("enrichPathway.png",height=900,width=2300)
ggarrange(plotsList[[1]][[1]][[4]],plotsList[[2]][[1]][[4]],plotsList[[3]][[1]][[4]], plotsList[[4]][[1]][[4]], ncol=4)
dev.off()


png("enrichKEGG.png",height=900,width=2300)
ggarrange(plotsList[[1]][[1]][[3]],plotsList[[2]][[1]][[3]],plotsList[[3]][[1]][[3]], plotsList[[4]][[1]][[3]], ncol=4)
dev.off()

ggarrange(plotsList[[1]][[1]][[1]],plotsList[[2]][[1]][[1]],
          plotsList[[3]][[1]][[1]],plotsList[[4]][[1]][[1]], ncol=4)
























##### TESTOWE PODEJSCIE -- obiekt dds z tylko 2 condition i bez genow bialkowych #####


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
library(sva)

setwd("/home/rstudio/all/Opus/cellLinesBUlkRNAseq/microglia")

httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- NULL
while(is.null(mart)){
     try(mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl"))
}


counts <- read.table("feature_counts_BAM", head=TRUE, row.names = 1)


colnames(counts) <- sub(pattern = "Aligned.sortedByCoord.out.bam", replacement = "",colnames(counts))

exon_length <- counts$Length
counts <- counts[,-c(1,2,3,4,5)]

corBetweenLines <- lapply(seq(1, length(counts[1,]), 2),function(i){
     cor(counts[,i],counts[,i+1],method="spearman")
})
corBetweenLines <- unlist(corBetweenLines)
mean(corBetweenLines)

min(corBetweenLines)

counts <- data.frame(counts)

tmp <- factor(colnames(counts))

summedLinesList <- lapply(seq(1, length(counts[1,]), 4),function(i){ 
     return(counts[,as.numeric(i)] + counts[,as.numeric(i+1)] + counts[,as.numeric(i+2)] + counts[,as.numeric(i+3)]) #sum sequencing of the same sample from 4 different lines on NovaSeq
})

countsSummed <- do.call(cbind,summedLinesList)

rownames(countsSummed) <- rownames(counts)
colnames(countsSummed) <- colnames(counts)[seq(1,length(colnames(counts)),4)]
colnames(countsSummed) <- gsub('.{9}$','',colnames(countsSummed))

countsSummed <- countsSummed[,-1] #outlier removal

############## BATCH EFFECT REMOVAL ################

batch <- c(rep(1,4),rep(2,5),rep(3,5),rep(4,5))

adjusted <- ComBat_seq(countsSummed, batch=batch, group=NULL)

################################

countsSummed <- adjusted



seqname <- read.csv("colnames.csv", header = FALSE)
seqname <- seqname[-1,]

eh <- colnames(countsSummed)

eh <- lapply(eh,function(i){
     paste(strsplit(i,"_")[[1]][c(1,2,3,4,5)],collapse ="_")
})

eh <- unlist(eh)
seqname <- seqname[match(eh,seqname$V1),]

colnames(countsSummed) <- seqname[,2]

RNA_all <- countsSummed

RNA_all_matrix <- matrix(as.numeric(unlist(RNA_all)),nrow=nrow(RNA_all))
rownames(RNA_all_matrix) <- rownames(counts)
colnames(RNA_all_matrix) <- colnames(RNA_all)
RNA_all_matrix <- RNA_all_matrix[rowSums(RNA_all_matrix) >= 160,] #normally 10reads per sample


gen <- getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description"), values = rownames(RNA_all), filters = "ensembl_gene_id", mart = mart, useCache = FALSE)

rownamesHugo <- gen[match(rownames(RNA_all_matrix),gen[,1]),2]
rownames(RNA_all_matrix) <- rownamesHugo

#condition <- factor(gsub('.{2}$','',colnames(RNA_all_matrix)))
coldata <- readRDS("coldata.rds")
coldata <- coldata[-1,]
rownames(coldata) <- colnames(RNA_all_matrix)

# select only two condition
RNA_all_matrix <- RNA_all_matrix[,c(1,2,6,7,11,12,16,17)]
coldata <- coldata[c(1,2,6,7,11,12,16,17),]
# select only protein coding 

RNA_all_matrix <-  cbind(RNA_all_matrix,gen$gene_biotype[match(rownames(RNA_all_matrix),gen$external_gene_name)])
colnames(RNA_all_matrix)[[9]] <- "biotype"
RNA_all_matrix[,9] <- RNA_all_matrix[,9]=="protein_coding"
RNA_all_matrix <- RNA_all_matrix[,1:8]
RNA_all_matrixChar <- RNA_all_matrix
RNA_all_matrix <- matrix(as.numeric(RNA_all_matrixChar),  
                                    ncol = ncol(RNA_all_matrixChar))
colnames(RNA_all_matrix) <- colnames(RNA_all_matrixChar)
rownames(RNA_all_matrix) <- rownames(RNA_all_matrixChar)


dds <- DESeqDataSetFromMatrix(countData=RNA_all_matrix, colData=coldata, design= ~ batch + condition)
dds<- DESeq(dds)

my_vst <- DESeq2::vst(dds, blind=FALSE)

plotPCA(my_vst, intgroup = "batch")
plotPCA(my_vst, intgroup = "condition")



results <- results(dds, contrast=c("condition", "o-NF53-wild-type", "NF53"), independentFiltering = TRUE, pAdjustMethod="BH")
m1<-match(results@rownames,gen[,1])
results1<-cbind(hugo = gen[m1,2], log2FoldChange = results$log2FoldChange, pvalue = results$pvalue, padj = results$padj)
rownames(results) <- results@rownames

results <- data.frame(results)
results$log2FoldChange <- as.numeric(results$log2FoldChange)
results$pvalue <- as.numeric(results$pvalue)
results$padj <- as.numeric(results$padj)



m<-match(rownames(results),gen$external_gene_name)

res<-cbind(rownames(results),gen[m,],results$log2FoldChange,results$pvalue,results$padj)

colnames(res)<-c("gens",colnames(gen),"log2FoldChange","pvalue","padj")


w<-which(results$padj<0.05 & (results$log2FoldChange> 0.5 | results$log2FoldChange< -0.5))

volc_plot <- EnhancedVolcano(res,lab=res$external_gene_name,x="log2FoldChange",y="padj",
                             pCutoff=0.05, FCcutoff=0.5, labSize=3, #boxedLabels = TRUE,
                             drawConnectors = TRUE, widthConnectors = 0.75,
                             maxoverlapsConnectors = 10,
                             selectLab = res$external_gene_name[which((res$log2FoldChange > 0.5 | res$log2FoldChange < -0.5) & res$padj < 0.05)],
                             xlim = c(-5,5))  + 
     theme(legend.position="none") + ggtitle("o-NF53 vs NF53")
volc_plot


wp<-which(results$padj<0.05 & results$log2FoldChange> 0.25)
wn<-which(results$padj<0.05 & results$log2FoldChange< -0.25)

kk1p <- enrichGO(gene         = res$entrezgene_id[wp],OrgDb="org.Mm.eg.db",
                 pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP", readable=TRUE)
kk1n <- enrichGO(gene         = res$entrezgene_id[wn],OrgDb="org.Mm.eg.db",
                 pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP", readable=TRUE)


kk1p_plot <- barplot(kk1p, showCategory=25) + ggtitle("o-NF53 vs NF53(overexpressed genes)")
kk1n_plot <- barplot(kk1n, showCategory=25) + ggtitle("o-NF53 vs NF53 (underexpressed genes)")
print(kk1p_plot)
print(kk1n_plot)


GENELIST <- res$log2FoldChange
names(GENELIST)<-res$entrezgene_id

cnetplot(kk1p, foldChange=GENELIST) + ggtitle("o-NF53 vs NF53")
cnetplot(kk1n, foldChange=GENELIST) + ggtitle("o-NF53 vs NF53")



















