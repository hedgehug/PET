# this script is used to run PIANO
# reference from: https://bioconductor.org/packages/release/bioc/html/piano.html

# Load piano and example data:
library(piano)

# provide the prerank file
pval_table <- read.table('rank_file_name')
pvals <- setNames(pval_table$V2,pval_table$V1)
abs_pvals <- unlist(lapply(pvals, abs))
original_pval <- unlist(lapply(abs_pvals, function(i) 10^(-i)))

directions <- lapply(pval_table$V2, sign)
directions <- unlist(directions)
directions <- setNames(directions,pval_table$V1)

pathways <- loadGSC('pathway_file')

# Run gene-set analysis:
gsares <- runGSA(geneLevelStats=original_pval,
                 directions = directions,
                 gsc = pathways)

# save results
res <- list(names(pathways$gsc), gsares$pAdjDistinctDirUp, gsares$pAdjDistinctDirDn)
my_df <- data.frame(res)
write.table(res, file = "result_file_name", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
