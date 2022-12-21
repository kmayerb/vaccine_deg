
library(tidyverse)
library(RColorBrewer)  
library(limma)
library(edgeR) 
library(RColorBrewer)
library(dplyr)
library('variancePartition')
library('edgeR')
library('BiocParallel')

data_file<-'X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_counts.csv'
treatment_file<-'trt_pubid_2022-DEC-19.csv'
sample_meta_data<-'X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_metadata.csv'
day_subset<-c("0", "3")


data_preprocessing <- function(data_file,sample_meta_data,treatment_file,
                               day_subset)
  {
  data<-read_csv(
          data_file,
          col_names = TRUE
              ) %>% 
          column_to_rownames(var = "gene_id")
  
  
  trt<-read.csv(treatment_file)
  trt<-trt[trt$Treatment.Group %in% 
                    c("2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)", 
                      "2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)"), ]
  trt<-trt%>% 
    rename("Subject_ID" = "Subject.ID")
  
  sample<-read_csv(
    sample_meta_data,
    col_names = TRUE
          )%>% 
    select("samplename","ptid","day") %>% 
    filter(day %in% day_subset) %>%
    mutate(Subject_ID=as.numeric(gsub("_", "",gsub("P203_", "",
                                                   substr(samplename,1,
                             (tail(unlist(gregexpr('_',samplename))
                                   , n=1))-1)))))

  meta_data<-inner_join(sample,trt, by = "Subject_ID") %>% 
             column_to_rownames(var = "samplename")
  idx2 <- match(rownames(meta_data), colnames(data))
  count_matrix<- data[,idx2]
  count_matrix<-as.matrix(count_matrix[,colSums(is.na(count_matrix))==0])
  
  
  
  
  result_data<-list(data,sample,trt,meta_data,count_matrix)
  names(result_data) <- c("original_count_data","original_sample_meta_data",
  "subset_treatment_data","subset_meta_data","final_count_matrix")

  
return(result_data)
}


data_filtering <- function(count_matrix)
{
  cpm <- cpm(count_matrix)
  lcpm <- cpm(count_matrix, log=TRUE)
  # fitering  genes
  keep <- rowSums(cpm(count_matrix)>1) >= 15
  d0 <- count_matrix[keep,]
  nsamples <- ncol(d0)
  col <- brewer.pal(nsamples, "Paired")
  par(mfrow=c(1,2))
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  lcpm <- cpm(d0, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")

  return(d0)
}

data_normalization <- function(d0)
{
  # creating a DGE object
  d0 <- DGEList(d0)
  d0 <- calcNormFactors(d0,method = "TMM")
  return(d0)
  
}



edgeR_preliminaries <- function(input_directory, comparisons, pipeline_input)
  {

  output_directory <- input_directory %>% 
    dplyr::first() %>% 
    paste0(pipeline_input$analysis_type, '/')
  
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  
  
  
  
  
  
  ty<-data_preprocessing(data_file,sample_meta_data,treatment_file,
                         day_subset)
  data1<-data_filtering(ty$final_count_matrix)
  d_ss<-data_normalization(data1)
  day_s=as.character(ty$subset_meta_data$day)
  
  meta_data<-ty$subset_meta_data
  
  design <- model.matrix(~day_s, d_ss$samples)
  # colnames(design) <- gsub("group", "", colnames(design))
  
  v <- voom(d_ss, design, plot=TRUE)
  
  #contr.matrix <- makeContrasts(vax=day3, levels=colnames(design))
  
  corfit <- duplicateCorrelation(v, design, block=meta_data$pubid)
  v <- voom(d_ss, design, block=meta_data$pubid, correlation = corfit$consensus)
  fit <- lmFit(v, design, block=meta_data$pubid, correlation = corfit$consensus)
  fit_test <- eBayes(fit)
  dt <- decideTests(fit_test)
  

  


  
  

}














ty<-data_preprocessing(data_file,sample_meta_data,treatment_file,
                       day_subset)
data1<-data_filtering(ty$final_count_matrix)
d_ss<-data_normalization(data1)
day_s=as.character(ty$subset_meta_data$day)

meta_data<-ty$subset_meta_data

design <- model.matrix(~day_s, d_ss$samples)
# colnames(design) <- gsub("group", "", colnames(design))

v <- voom(d_ss, design, plot=TRUE)

#contr.matrix <- makeContrasts(vax=day3, levels=colnames(design))

corfit <- duplicateCorrelation(v, design, block=meta_data$pubid)
v <- voom(d_ss, design, block=meta_data$pubid, correlation = corfit$consensus)
fit <- lmFit(v, design, block=meta_data$pubid, correlation = corfit$consensus)
fit_test <- eBayes(fit)
dt <- decideTests(fit_test)








