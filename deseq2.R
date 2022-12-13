library("DESeq2")

# input data matrix
data <- as.matrix(read.csv(file = 'X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_merged_gene_counts.csv', row.names = "gene_id"))
colnames(data) = gsub("X", "P", colnames(data))
trt_data <- read.csv('trt.csv', header=FALSE, stringsAsFactors=FALSE, fileEncoding="latin1")
colnames(trt_data) = trt_data[1,1:5]
trt_data=trt_data[-1,]
trt_data<-trt_data[trt_data$`Treatment Group` %in% c("2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)", "2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)"), ]



#input sample metadata
samples <- read.csv('X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_metadata.csv')
samples<-samples[,2:4]
rownames(samples) <- samples$samplename


#input sample metadata
samples_0_3<-samples
samples_0_3$Subject_ID<-substr(samples_0_3$samplename,1,(tail(unlist(gregexpr('_', samples_0_3$samplename)), n=1))-1)
samples_0_3$Subject_ID =gsub("P203_", "", samples_0_3$Subject_ID)
samples_0_3$Subject_ID=gsub("_", "", samples_0_3$Subject_ID)
samples_0_3<-samples_0_3[samples_0_3$Subject_ID %in% trt_data$`Subject ID`, ]
tt<-merge(trt_data,samples_0_3, by.x= "Subject ID", by.y="Subject_ID")
tt_0_3<-tt[tt$day %in% c("0", "3"), ]
tt_0_3$day <- factor(tt_0_3$day)
rownames(tt_0_3)<-tt_0_3$samplename



#subset data
idx <- match(rownames(tt_0_3), colnames(data))
data_0_3<- data[,idx]
all(rownames(tt_0_3) == colnames(data_0_3))


#create a DESEQ object
dds <- DESeqDataSetFromMatrix(countData = round(data_0_3),
                              colData = tt_0_3,
                              design= ~ ptid+day)

#fiter genes
keep <- rowSums(counts(dds) >= 10) >= 10
sum(keep)
dds <- dds[keep,]
nrow(dds)



#estimate dispersion in DEseq obect
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
dds = estimateDispersions( dds )



#plot dispersion estimate

plotDispEsts(dds)
cds <- DESeq(dds)
res <- results(cds)
table(res$padj<0.1)


#find out significant genes
resSigind = res[ which(res$padj < 0.1 & res$log2FoldChange > 0), ]
resSigrep = res[ which(res$padj < 0.1 & res$log2FoldChange < 0), ]
resSig = rbind(resSigind, resSigrep)
rownames(resSigind)


# order the result from deseq based on pvalues
res <- res[order(res$padj), ]

## Merge with normalized count data
resdata <- merge(as.data.frame(res), 
                 as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", 
                 sort = FALSE)
names(resdata)[1] <- "Gene"
head(resdata)


## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## This is Stephen Turner's code:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), points(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot_3_vs_0.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot_3_vs_0")
dev.off()




## Write results
write.csv(resdata, file = "diffexpr-results_63_vs_3.csv", quote = FALSE, row.names = F)

write.csv(resSig, file = "significant_genes-results_3_vs_0.csv", quote = FALSE, row.names = T)




















