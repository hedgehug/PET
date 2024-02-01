suppressWarnings({library(DESeq2)})
countdata <- read.table("example/example_data.txt", header=TRUE, row.names = 1)

countdata <- countdata[ ,1:10]

# Convert to matrix
countdata <- as.matrix(countdata)

# Assign condition
condition <- factor(c(rep("cond1", 4),rep("cond2", 3),rep("cond3", 3)))

# Analysis with DESeq2 ----------------------------------------------------

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)
SAVE_ALL <- FALSE
if (SAVE_ALL){
  # save results of all comparisons
  for (idx_1 in 1:(nlevels(dds$condition)-1)) {
    for (idx_2 in (idx_1+1):nlevels(dds$condition)){
      print(paste('Comparing', levels(dds$condition)[idx_1], 'vs.', levels(dds$condition)[idx_2]))
      res <- results(dds, contrast=c("condition", levels(dds$condition)[idx_1], levels(dds$condition)[idx_2]))
      # table(res$padj<0.05)
      # Order by adjusted p-value
      res <- res[order(res$padj), ]
      resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
      names(resdata)[1] <- "Gene"
      write.csv(res, file=paste("OUT_DIR/", levels(dds$condition)[idx_1], '.vs.', levels(dds$condition)[idx_2],'.csv', sep =''))
      print(paste('Result saved to ', paste("OUT_DIR/", levels(dds$condition)[idx_1], '.vs.', levels(dds$condition)[idx_2],'.csv', sep ='')))
  }
  }
}else{
  # save results for one comparison
  # Get differential expression results
  print(paste('Comparing', "cond2", 'vs.', "cond1"))
  res <- results(dds, contrast=c("condition", "cond2", "cond1"))
  # table(res$padj<0.05)
  ## Order by adjusted p-value
  res <- res[order(res$padj), ]
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  write.csv(res, file="example/cond2.vs.cond1.deseq_result.txt")
  print(paste('Result saved to ', "example/cond2.vs.cond1.deseq_result.txt"))
}
