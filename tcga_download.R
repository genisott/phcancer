library(TCGAbiolinks)
library(SummarizedExperiment)


ctype = "ACC"

query <- GDCquery(
  project = paste("TCGA-",ctype,sep="")
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  platform = "Illumina HiSeq", 
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(query, method = "api", files.per.chunk = 10)
data <- assay(GDCprepare(query))

write.table(data,paste("expression_",ctype,".txt",sep=""),sep=",")
