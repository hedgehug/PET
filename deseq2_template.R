suppressWarnings({library(DESeq2)})
countdata <- read.table("data_matrix_file_name", header=TRUE, row.names = 1)

countdata <- countdata[ ,1:Sample_num]

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(c(rep("group1_id", group1_num), rep("group2_id", group2_num))))

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
# table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
write.csv(res, file="OUTPUT.csv")