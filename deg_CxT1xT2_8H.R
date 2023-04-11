#Analises de DEGs do augusto CONTROL X TRATAMENTO 8H
library(DESeq2)
getwd()

cts = read.delim("ReadCount_aug_usp.tab", sep="", header = T)
cts = cts[,-(1:3)]
cts = cts[,-(4:6)]
cts = cts[,-(7:9)]
dim(cts)


Name=colnames(cts)
Time=c(rep("case2",3),rep("control", 3), rep("case1",3))
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
head(res)
dim(res)
data.res1 <- as.data.frame(res1)
data.res2 <- as.data.frame(res2)
data.res3 <- as.data.frame(res3)
write.csv(data.res1,"res_aug_CxC1_8H.csv", row.names = T)
write.csv(data.res2,"res_aug_CxC2_8H.csv", row.names = T)
write.csv(data.res3,"res_aug_T1xT2_8H.csv", row.names = T)
