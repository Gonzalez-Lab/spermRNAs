library(dplyr)
library(org.Mm.eg.db)

#Ribosome-bound RNA at 1 cell embryo
cell1rep1 <- read.delim("liRiboseq_1-cell_rep1.quant.genes.txt", header=TRUE)
cell1rep2 <- read.delim("liRiboseq_1-cell_rep2.quant.genes.txt", header=TRUE)

#create 1 cell embryo riboseq table
embryo_1c <- cbind(rep1=cell1rep1$FPKM, rep2=cell1rep2$FPKM)
embryo_1c <- cbind(embryo_1c, mean=rowMeans(embryo_1c))
embryo_1c <- as.data.frame(cbind(embryo_1c,ensembl=cell1rep1$gene_id))

embryo_1c$symbol <- mapIds(org.Mm.eg.db, keys = substr(embryo_1c$ensembl,1,18), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
embryo_1c[which(embryo_1c$symbol=="N6amt1"),]

#total RNA at 1 cell embryo
cell1rep1T <- read.delim("totalRNA_1-cell_rep1.quant.genes.txt", header=TRUE)
cell1rep2T <- read.delim("totalRNA_1-cell_rep2.quant.genes.txt", header=TRUE)

#create 1 cell embryo total RNA table
embryo_1cT <- cbind(rep1=cell1rep1T$FPKM, rep2=cell1rep2T$FPKM)
embryo_1cT <- cbind(embryo_1cT, mean=rowMeans(embryo_1cT))
embryo_1cT <- as.data.frame(cbind(embryo_1cT,ensembl=cell1rep1T$gene_id))

embryo_1cT$symbol <- mapIds(org.Mm.eg.db, keys = substr(embryo_1cT$ensembl,1,18), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
embryo_1cT[which(embryo_1cT$symbol=="N6amt1"),]

embryo_1cT <- embryo_1cT[which(embryo_1cT$ensembl %in% embryo_1c$ensembl),]

embryo_1cT <- embryo_1cT[match(embryo_1c$symbol, embryo_1cT$symbol),]

RPKM_TE_table <- cbind(totalRNA=embryo_1cT$mean,translRNA=embryo_1c$mean)
RPKM_TE_table <- cbind(RPKM_TE_table, TE=as.numeric(RPKM_TE_table[,2])/as.numeric(RPKM_TE_table[,1]))
RPKM_TE_table <- as.data.frame(RPKM_TE_table)
RPKM_TE_table$symbol <- embryo_1cT$symbol

RPKM_TE_table[which(RPKM_TE_table$symbol=="N6amt1"),]

tranl_genes <- RPKM_TE_table[which(RPKM_TE_table$TE>0),]$symbol


########################################################################################################
####################################################################################################
#MII embryo analysis

#total RNA at MII oocyte
MIIrep1T <- read.delim("totalRNA_mII_rep1.quant.genes.txt", header=TRUE)
MIIrep2T <- read.delim("totalRNA_MII_rep2.quant.genes.txt", header=TRUE)
#
MII <- as.data.frame(cbind(mii.rep1=MIIrep1T$FPKM, mii.rep2=MIIrep2T$FPKM, ensembl=MIIrep1T$gene_id))
MII$symbol <- mapIds(org.Mm.eg.db, keys = substr(MII$ensembl,1,18), keytype = "ENSEMBL", column="SYMBOL")

embryoT <- embryo_1cT[which(embryo_1cT$symbol %in% MII$symbol),]
embryoT <- embryoT[which(embryoT$symbol!="NA"),]

embryoT <- embryoT %>% distinct(symbol, .keep_all= TRUE)

MIIT <- MII[which(MII$symbol %in% embryoT$symbol),]
MIIT <- MIIT[which(MIIT$symbol!="NA"),]

MIIT <- MIIT %>% distinct(symbol, .keep_all= TRUE)

embryoT <- embryoT[match(embryoT$symbol,MIIT$symbol),]


MII.embryo <- cbind(Mii.rep1 = as.numeric(MIIT$mii.rep1),
                    Mii.rep2 = as.numeric(MIIT$mii.rep2),
                    embryo.rep1 = as.numeric(embryoT$rep1),
                    embryo.rep2 = as.numeric(embryoT$rep2))



##
rownames(MII.embryo) <- MIIT$symbol
MII.embryo[which(rownames(MII.embryo)=="Pnrc1"),]

##

library(DESeq2)

count.matrix <-round(MII.embryo,0)
colnames(count.matrix) <- c("MII_1","MII_2","embryo_1","embryo_2")
count.matrix <- as.matrix(count.matrix)

sample.table <- as.data.frame(matrix(c("mii","mii","embryo","embryo"), ncol = 1))
rownames(sample.table) <- c("MII_1","MII_2","embryo_1","embryo_2")
colnames(sample.table) <- c("group")

sample.table$group <- factor(sample.table$group)

levels(sample.table$group) 

#chequeo que las tablas estan bien
all(rownames(sample.table) %in% colnames(count.matrix))
all(rownames(sample.table) == colnames(count.matrix))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds <- DESeqDataSetFromMatrix(countData=count.matrix, 
                               colData=sample.table, 
                               design=~group)

# al objeto dds (que es el dataset) le agregamos los analisis
dds <- DESeq(dds)

res <- results(dds)

table(res$pvalue < 0.05 & res$log2FoldChange < 0) # upregulados en embryo vs MII
table(res$padj < 0.05 & res$log2FoldChange < 0) # upregulados en embryo vs MII


#le agrego los genes a la tabla de resultados
res$symbol <- MIIT$symbol

HighEmbryo <- res[which(res$pvalue<0.05 & res$log2FoldChange < 0),]$symbol

res[which(res$symbol=="Kmt2c"),]

write.csv2(HighEmbryo, "genes high in embryo vs mii.csv")
