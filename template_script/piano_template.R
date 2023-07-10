# Load piano and example data:
library(piano)

pval_table <- read.table('rank_file')
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
write.table(res, file = "result_file", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
