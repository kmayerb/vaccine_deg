

library(limma)
library(edgeR) 
library(RColorBrewer)

# input data matrix
data <- round(as.matrix(read.csv(file = 'X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_counts.csv', row.names = "gene_id")))
colnames(data) = gsub("X", "P", colnames(data))



# input treatment data
trt_data <- read.csv('trt.csv', header=FALSE, stringsAsFactors=FALSE, fileEncoding="latin1")
colnames(trt_data) = trt_data[1,1:5]
trt_data=trt_data[-1,]
trt_data<-trt_data[trt_data$`Treatment Group` %in% c("2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)", "2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)"), ]


# input sample meta data
samples <- read.csv('X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_metadata.csv')
samples<-samples[,2:4]
rownames(samples) <- samples$samplename
samples_0_3<-samples
samples_0_3$Subject_ID<-substr(samples_0_3$samplename,1,(tail(unlist(gregexpr('_', samples_0_3$samplename)), n=1))-1)
samples_0_3$Subject_ID =gsub("P203_", "", samples_0_3$Subject_ID)
samples_0_3$Subject_ID=gsub("_", "", samples_0_3$Subject_ID)
samples_0_3<-samples_0_3[samples_0_3$Subject_ID %in% trt_data$`Subject ID`, ]
tt<-merge(trt_data,samples_0_3, by.x= "Subject ID", by.y="Subject_ID")



#subset data
idx2 <- match(rownames(samples_0_3), colnames(data))
data_4<- data[,idx2]
data_4<-as.matrix(data_4[ , colSums(is.na(data_4))==0])



#transforming the data
cpm <- cpm(data_4)
lcpm <- cpm(data_4, log=TRUE)




# fitering  genes
keep <- rowSums(cpm(data_4)>1) >= 2
d0 <- data_4[keep,]
dim(d0)



#creating a group
trt_grp <- substr(tt$`Treatment Group`, 36, nchar(tt$`Treatment Group`)-1) 
day <- tt$day
group <- interaction(trt_grp,day)
group<-factor(group)


nsamples <- ncol(d0)

# creating a DGE object
d0 <- DGEList(d0)
d0 <- calcNormFactors(d0,method = "TMM")


day<-as.character(tt$day)
ptid<-as.character(tt$`Subject ID`)

design <- model.matrix(~0+day+ptid)
colnames(design) <- gsub("group", "", colnames(design))


contr.matrix <- makeContrasts(
  zerovsthree = day3-day0,
  fiftysixvsthree = day56-day0,
  levels = colnames(design))
contr.matrix


v <- voom(d0, design, plot=TRUE)
v


vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
dt <- decideTests(efit)
summary(dt)
three.vs.zero <- topTreat(efit, coef=1, n=Inf)
fifty_six.vs.zero <- topTreat(efit, coef=2, n=Inf)


plotMD(efit, column=1, status=dt[,2], main=colnames(efit)[2], xlim=c(-8,13))
## Write results
write.csv(three.vs.zero, file = "edge_r_diffexpr-results_3_vs_0.csv", quote = FALSE, row.names = F)
write.csv(fifty_six.vs.zero, file = "edge_r_diffexpr-results_56_vs_0.csv", quote = FALSE, row.names = F)


tfit <- treat(vfit, lfc=1)
summary(decideTests(tfit))







