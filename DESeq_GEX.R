library(DESeq2)
library(data.table)

directory <- "C:\\Users\\josep\\Documents\\ETHZ\\Boeva_lab\\CIMP_unsupervised_cancer\\code\\gene_counts"
all_files <- list.files(path=directory)
#### FOR ALL TCGA CANCER TYPES

for(file in all_files){
  print(file)
  cancer <- strsplit(strsplit(file,split="-")[[1]][2],split=".htseq")[[1]][1]
  cts <- as.data.frame(fread(paste(directory,file,sep="\\")))

  file_clust <- paste0("C:\\Users\\josep\\Documents\\ETHZ\\Boeva_lab\\CIMP_unsupervised_cancer\\code\\cluster_memb\\",cancer,"_cluster_memb_new.csv")
  coldata <- as.data.frame(fread(file_clust))
  if(cancer=="LAML"){
    colnames(cts) <- sub("03A", "03", colnames(cts))
  }
  else{
    colnames(cts) <- sub("01A", "01", colnames(cts))
  }
  rownames(cts) <- cts$"Ensembl_ID"
  if("Unnamed: 0" %in% colnames(coldata)){
    rownames(coldata) <- coldata$`Unnamed: 0`
    coldata <- coldata[,c("Unnamed: 0","cluster")]
  }
  else{
    rownames(coldata) <- coldata$`index`
    coldata <- coldata[,c("index","cluster")]
  }
  

  patients <- intersect(colnames(cts),rownames(coldata))
  cts <- cts[,patients]
  coldata <- coldata[patients, ]

  print(all(rownames(coldata) == colnames(cts)))
  for(col in colnames(cts)){
    cts[,col] <- as.integer(cts[,col])
  }

  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ cluster)

  if(cancer %in% c('ACC', 'BRCA', 'CESC', 'COAD', 'GBM', 'ESCA', 'HNSC', 'KIRC', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'PAAD', 'PCPG', 'READ', 'SKCM', 'THYM', 'UCEC')){
    dds$cluster <- factor(dds$cluster, levels = c(1,2,3))
  }
  else{
    dds$cluster <- factor(dds$cluster, levels = c(1,2))
  }


  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]


  dds <- DESeq(dds)
  res <- results(dds)
  res

  res$pvalue[is.na(res$pvalue)] <- 1


  write.csv(as.data.frame(res[which(res$pvalue<0.05),c("log2FoldChange","pvalue")]),file=paste0("C:\\Users\\josep\\Downloads\\",cancer,"_diff_gene_expression.csv"))

}

