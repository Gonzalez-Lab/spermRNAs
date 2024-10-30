library(dplyr)
library(org.Mm.eg.db)

#upload tables
#Ribosome-bound RNA at zygote
zyg.ribo.rep1 <- read.delim("liRiboseq_1-cell_rep1.quant.genes.txt", header=TRUE)
zyg.ribo.rep2 <- read.delim("liRiboseq_1-cell_rep2.quant.genes.txt", header=TRUE)
#total RNA at zygote
zyg.T.rep1 <- read.delim("totalRNA_1-cell_rep1.quant.genes.txt", header=TRUE)
zyg.T.rep2 <- read.delim("totalRNA_1-cell_rep2.quant.genes.txt", header=TRUE)
#total RNA at MII oocyte
MII.T.rep1 <- read.delim("totalRNA_mII_rep1.quant.genes.txt", header=TRUE)
MII.T.rep2 <- read.delim("totalRNA_MII_rep2.quant.genes.txt", header=TRUE)

#################################################################################################
#Translated RNA at Zygote analysys

#create zygote riboseq table
zyg.ribo.RNA <- cbind(rep1=zyg.ribo.rep1$FPKM, rep2=zyg.ribo.rep2$FPKM)
zyg.ribo.RNA <- cbind(zyg.ribo.RNA, mean=rowMeans(zyg.ribo.RNA))
zyg.ribo.RNA <- as.data.frame(cbind(zyg.ribo.RNA,ensembl=zyg.ribo.rep1$gene_id))

zyg.ribo.RNA$symbol <- mapIds(org.Mm.eg.db, keys = substr(zyg.ribo.RNA$ensembl,1,18), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
zyg.ribo.RNA[which(zyg.ribo.RNA$symbol=="N6amt1"),]

#create zygote total RNA table
zyg.T.RNA <- cbind(rep1=zyg.T.rep1$FPKM, rep2=zyg.T.rep2$FPKM)
zyg.T.RNA <- cbind(zyg.T.RNA, mean=rowMeans(zyg.T.RNA))
zyg.T.RNA <- as.data.frame(cbind(zyg.T.RNA,ensembl=zyg.T.rep1$gene_id))

zyg.T.RNA$symbol <- mapIds(org.Mm.eg.db, keys = substr(zyg.T.RNA$ensembl,1,18), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
zyg.T.RNA[which(zyg.T.RNA$symbol=="N6amt1"),]

#select genes detected in both tables
zyg.ribo.RNA <- zyg.ribo.RNA[which(zyg.ribo.RNA$ensembl %in% zyg.T.RNA$ensembl),]

#match the order
zyg.T.RNA <- zyg.T.RNA[match(zyg.ribo.RNA$symbol, zyg.T.RNA$symbol),]

RPKM_TE_table <- cbind(totalRNA = zyg.T.RNA$mean, translRNA = zyg.ribo.RNA$mean)
RPKM_TE_table <- cbind(RPKM_TE_table, TE = as.numeric(RPKM_TE_table[,2])/as.numeric(RPKM_TE_table[,1]))
RPKM_TE_table <- as.data.frame(RPKM_TE_table)
RPKM_TE_table$symbol <- zyg.T.RNA$symbol

RPKM_TE_table[which(RPKM_TE_table$symbol=="N6amt1"),]

#select translated genes in zygote
tranl_genes <- RPKM_TE_table[which(RPKM_TE_table$TE>0),]$symbol

########################################################################################################
####################################################################################################
#MII zygote DESeq2 analysis

#create MII and zygote total RNA counts table

MII.embryo <- cbind(Mii.rep1 = as.integer(MII.T.rep1$expected_count),
                    Mii.rep2 = as.integer(MII.T.rep2$expected_count),
                    embryo.rep1 = as.integer(zyg.T.rep1$expected_count),
                    embryo.rep2 = as.integer(zyg.T.rep2$expected_count))

##
library(DESeq2)

count.matrix <- MII.embryo
colnames(count.matrix) <- c("MII_1","MII_2","embryo_1","embryo_2")
count.matrix <- as.matrix(count.matrix)

sample.table <- as.data.frame(matrix(c("mii","mii","embryo","embryo"), ncol = 1))
rownames(sample.table) <- c("MII_1","MII_2","embryo_1","embryo_2")
colnames(sample.table) <- c("group")

sample.table$group <- factor(sample.table$group)

levels(sample.table$group) 

#check tables
all(rownames(sample.table) %in% colnames(count.matrix))
all(rownames(sample.table) == colnames(count.matrix))

#We make a *DESeqDataSet* (dds) from a count matrix and column data
dds <- DESeqDataSetFromMatrix(countData=count.matrix, 
                               colData=sample.table, 
                               design=~group)

# add analysis to dds object
dds <- DESeq(dds)

#get results
res <- results(dds)

## UP genes in zygote vs MII
table(res$padj < 0.1 & res$log2FoldChange < -0.5) 

#add gene symbol to res table
res$symbol <- mapIds(org.Mm.eg.db, keys = substr(MII.T.rep1$gene_id,1,18), keytype = "ENSEMBL", column="SYMBOL")

#inspect gene of interest
gene <- "Smarca2"
res[which(res$symbol==gene),]

#select genes higher in zygote vs MII from res table
HighEmbryo <- res[which(res$padj<0.1 & res$log2FoldChange < -0.5),]$symbol

#remove NAs
HighEmbryo <- HighEmbryo[which(HighEmbryo!="NA")]
