#Analises de DEGs do augusto CONTROL X TRATAMENTO 24H
library(DESeq2)
getwd()
cts = read.delim("ReadCount_aug_usp.tab", sep="", header = T)
cts = cts[,-(4:6)]
cts = cts[,-(7:9)]
cts = cts[,-(10:12)]
dim(cts)


Name=colnames(cts)
Time=c(rep("case1",3),rep("case2", 3), rep("control",3))
coldata=data.frame(Name,Time)
head(coldata,18)
dim(coldata)

log_cts<- log(cts+1, 10) # This is log base 10 + 1 for "0"
table(rowSums(log_cts>log10(3))>=9)

keep.exprs<- rowSums(log_cts>log10(3))>=9
cts_filt<-cts[keep.exprs,]
dim(cts_filt)
head(cts_filt,10)

dds <- DESeqDataSetFromMatrix(countData = round(cts_filt), colData = coldata, design = ~ Time)
dds$Time <- factor(dds$Time, levels = c("control", "case1","case2"))
head(dds)

dds <- DESeq(dds)
res1 <- results(dds, contrast=c("Time","case1","control"))
res2 <- results(dds, contrast=c("Time","case2","control"))
res3 <- results(dds, contrast=c("Time","case2","case1"))


data.res1 <- as.data.frame(res1)
data.res2 <- as.data.frame(res2)
data.res3 <- as.data.frame(res3)
write.csv(data.res1,"res_aug_CxC1_24H.csv", row.names = T)
write.csv(data.res2,"res_aug_CxC2_24H.csv", row.names = T)
write.csv(data.res3,"res_aug_T1xT2_24H.csv", row.names = T)

##Transformar o ensemblid
#data.res1:

library(biomaRt)
listEnsembl()
datasets <- listDatasets(ensembl)
attr <- listAttributes(ensembl.con)


ensembl.ids <- rownames(data.res1)
ensembl <- useEnsembl(biomart = "genes")

ensembl.con <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")


ids_1 <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'ensembl_gene_id', 
                            'external_transcript_name',
                            'external_gene_name'),
             filters = 'ensembl_transcript_id_version', 
             values = ensembl.ids,
            mart = ensembl.con)

#combinar data.res1 e res: 

res_ordem1 <- ids_1 %>% arrange(ids_1$ensembl_transcript_id_version)
final1 <- cbind(data.res1, res_ordem1)


#data.res2:

library(biomaRt)
ensembl.ids2 <- rownames(data.res2)
ensembl <- useEnsembl(biomart = "genes")
ensembl.con <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")


ids_2 <- getBM(attributes = c('ensembl_transcript_id_version', 
                            'ensembl_gene_id', 
                            'external_transcript_name',
                            'external_gene_name'),
             filters = 'ensembl_transcript_id_version', 
             values = ensembl.ids2,
             mart = ensembl.con)

#combinar data.res2 e ids_2: 

res_ordem2 <- ids_2 %>% arrange(ids_2$ensembl_transcript_id_version)
final2 <- cbind(data.res2, res_ordem2)


#data.res3:

library(biomaRt)
ensembl.ids3 <- rownames(data.res3)
ensembl <- useEnsembl(biomart = "genes")
ensembl.con <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")


ids_3 <- getBM(attributes = c('ensembl_transcript_id_version', 
                              'ensembl_gene_id', 
                              'external_transcript_name',
                              'external_gene_name'),
               filters = 'ensembl_transcript_id_version', 
               values = ensembl.ids3,
               mart = ensembl.con)

#combinar data.res3 e ids_3: 

res_ordem3 <- ids_3 %>% arrange(ids_3$ensembl_transcript_id_version)
final3 <- cbind(data.res3, res_ordem3)







##VOLCANOPLOTS


##VOLCANO PLOT CONTROLE X TRATAMENTO 1 24H

#Final1: Control x Case1 24h:

library(ggplot2)
library(ggrepel)
# Convert directly in the aes()
p <- ggplot(final1, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
# Add more simple "theme"
p <- ggplot(final1, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
final1$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
final1$diffexpressed[final1$log2FoldChange > 0.58 & final1$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
final1$diffexpressed[final1$log2FoldChange < -0.58 & final1$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(final1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
final1$delabel <- NA
final1$delabel[final1$diffexpressed != "NO"] <- final1$external_gene_name[final1$diffexpressed != "NO"]

ggplot(final1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() 
##MARCAR NOMES DOS GENES
  #+ geom_label_repel(aes(label=delabel), max.overlaps = Inf,size=1, segment.size=0.25, nudge_x=0.5, direction="y")


##SALVAR: 
pdf(file="volplot_cxc1_24h.pdf")
 a1 <- ggplot(final1, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal()+
  geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
 
a1 + ggtitle("Controle vs Tratamento 1 24H") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


#Final2: Control x Case2 24h:

library(ggplot2)
library(ggrepel)
# Convert directly in the aes()
p <- ggplot(final2, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
# Add more simple "theme"
p <- ggplot(final2, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
final2$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
final2$diffexpressed[final2$log2FoldChange > 0.58 & final2$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
final2$diffexpressed[final2$log2FoldChange < -0.58 & final2$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(final2, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
final2$delabel <- NA
final2$delabel[final2$diffexpressed != "NO"] <- final2$external_gene_name[final2$diffexpressed != "NO"]


##SALVAR: 
pdf(file="volplot_cxc2_24h.pdf")
a2 <- ggplot(final2, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal()+
  geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

a2 + ggtitle("Controle vs Tratamento 2 24H") + theme(plot.title = element_text(hjust = 0.5))
dev.off()
a2

#Final3: Case1 x Case2 24h:

library(ggplot2)
library(ggrepel)
# Convert directly in the aes()
p <- ggplot(res3, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
# Add more simple "theme"
p <- ggplot(res3, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
res3$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res3$diffexpressed[res3$log2FoldChange > 0.58 & res3$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res3$diffexpressed[res3$log2FoldChange < -0.58 & res3$pvalue < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(res3, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)


##SALVAR: 
pdf(file="volplot_c1xc2_24h.pdf")
a3 <- ggplot(res3, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal()+
  geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

a3 + ggtitle("Tratamento1 vs Tratamento 2 24H") + theme(plot.title = element_text(hjust = 0.5))
dev.off()
a3
