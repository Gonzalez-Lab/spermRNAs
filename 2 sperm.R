library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggwordcloud)
library(org.Mm.eg.db)

#upload tables
gene.lengths <- read.csv("gencodevM29.genelength.csv", header=TRUE)
gene.lengths$symbol <- mapIds(org.Mm.eg.db, keys = gene.lengths$gene, keytype = "ENSEMBL", column="SYMBOL")

##########################################################################################################
#table 1 is total and head sperm mRNA from GSE81216
sperm_table1 <- read.delim("GSE81216 sperm counts.tab", header=TRUE)

sperm_table1$gene.lengths <- gene.lengths$mean[match(sperm_table1$GeneId,gene.lengths$symbol)]

table1_rpkm <- cbind(total.1=(sperm_table1[,2]/(sum(sperm_table1[,2])/1000000))/(sperm_table1$gene.lengths/1000),
                     total.2=(sperm_table1[,3]/(sum(sperm_table1[,3])/1000000))/(sperm_table1$gene.lengths/1000),
                     head.1=(sperm_table1[,4]/(sum(sperm_table1[,4])/1000000))/(sperm_table1$gene.lengths/1000),
                     head.2=(sperm_table1[,5]/(sum(sperm_table1[,5])/1000000))/(sperm_table1$gene.lengths/1000))

table1_rpkm <- as.data.frame(table1_rpkm) 
table1_rpkm$mean_rpkm <- rowMeans(table1_rpkm)
table1_rpkm$symbol <- sperm_table1$GeneId
table1_rpkm$percentRank <- percent_rank(table1_rpkm$mean_rpkm)

###########################################################################################################
##########################################################################################################
#table 2 is total sperm mRNA from GSE88732
sperm_table2 <- read.delim("GSE88732 sperm counts.tab", header=TRUE)

sperm_table2$gene.lengths <- gene.lengths$mean[match(sperm_table2$GeneId,gene.lengths$symbol)]

table2_rpkm <- cbind(total.1=(sperm_table2[,2]/(sum(sperm_table2[,2])/1000000))/(sperm_table2$gene.lengths/1000),
                     total.2=(sperm_table2[,3]/(sum(sperm_table2[,3])/1000000))/(sperm_table2$gene.lengths/1000),
                     total.3=(sperm_table2[,4]/(sum(sperm_table2[,4])/1000000))/(sperm_table2$gene.lengths/1000),
                     total.4=(sperm_table2[,5]/(sum(sperm_table2[,5])/1000000))/(sperm_table2$gene.lengths/1000))


table2_rpkm <- as.data.frame(table2_rpkm) 
table2_rpkm$mean_rpkm <- rowMeans(table2_rpkm)
table2_rpkm$symbol <- sperm_table2$GeneId
table2_rpkm$percentRank <- percent_rank(table2_rpkm$mean_rpkm)


###########################################################################################################
##########################################################################################################
#table 3 is total sperm mRNA from E-MTAB-5834 control vs MSUS mice
sperm_table3 <- read.delim("E-MTAB-5834 sperm counts.tab", header=TRUE)

sperm_table3$gene.lengths <- gene.lengths$mean[match(sperm_table3$GeneId,gene.lengths$symbol)]

table3_rpkm <- cbind(total.1=(sperm_table3[,2]/(sum(sperm_table3[,2])/1000000))/(sperm_table3$gene.lengths/1000),
                     total.2=(sperm_table3[,3]/(sum(sperm_table3[,3])/1000000))/(sperm_table3$gene.lengths/1000),
                     total.3=(sperm_table3[,4]/(sum(sperm_table3[,4])/1000000))/(sperm_table3$gene.lengths/1000),
                     total.4=(sperm_table3[,5]/(sum(sperm_table3[,5])/1000000))/(sperm_table3$gene.lengths/1000))


table3_rpkm <- as.data.frame(table3_rpkm) 
table3_rpkm$mean_rpkm <- rowMeans(table3_rpkm)
table3_rpkm$symbol <- sperm_table3$GeneId
table3_rpkm$percentRank <- percent_rank(table3_rpkm$mean_rpkm)

###############################################################################################################
#mean percent rank for the three studies
sperm_tables <- as.data.frame(cbind(table1=table1_rpkm$percentRank, 
                                    table2=table2_rpkm$percentRank,
                                    table3=table3_rpkm$percentRank))

#calculate mean percent rank
sperm_tables$mean_percent_rank <- rowMeans(sperm_tables)

#add gene symbol to sperm tables
sperm_tables$symbol <- table1_rpkm$symbol

#select genes with positive detection
sperm_tables <- sperm_tables[which(sperm_tables$mean_percent_rank>0),]

#create column indicating 1: is translated (TE>0), 0: not translated
sperm_tables$TranslEmbryo <- ifelse(sperm_tables$symbol %in% tranl_genes,1,0)

#create column indicating 1: is higher in zygote vs MII, 0: not different
sperm_tables$HigherEmbryo <- ifelse(sperm_tables$symbol %in% HighEmbryo,1,0)

# select genes with TE>0 and higher in zygote
genes_sperm <- sperm_tables[which(sperm_tables$TranslEmbryo==1 & sperm_tables$HigherEmbryo==1),]

#obtain paternal genes 
paternal_genes <- genes_sperm[genes_sperm$mean_percent_rank > 0.70,]$symbol

write.csv2(paternal_genes, "paternal genes.csv")

######################################################################################
#wordcloud plot

embryo.genes <- c("Brd4","Chp1","Kpna2","Ncor1","Pabpn1","Hpf1","Tulp4","Ptges3",
                  "Cdc42se2","Smad1","Cfl2","Ybx1","Ncoa1","Spin1","Mbtd1","Srsf5",
                  "Dync1i2","Dpysl3","Bbs9","Srrm1","Rnps1")

sperm_tables$embryo <- ifelse(sperm_tables$symbol %in%  embryo.genes,1,0)

ggplot(data = sperm_tables[which(sperm_tables$symbol %in% paternal_genes$symbol),], 
       aes(label = symbol, col=as.character(embryo))) + 
  geom_text_wordcloud(rm_outside = TRUE, max_steps = 4,
                      grid_size = 5, eccentricity =0.9)+
  theme_void()+
  scale_color_manual(values = c("black", "#d52954"))
