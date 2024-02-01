# Load piano and example data:
library(piano)

pval_table <- read.table('PVAL_FILE')
pvals <- setNames(pval_table$V2,pval_table$V1)
abs_pvals <- unlist(lapply(pvals, abs))
original_pval <- unlist(lapply(abs_pvals, function(i) 10^(-i)))

directions <- lapply(pval_table$V2, sign)
directions <- unlist(directions)
directions <- setNames(directions,pval_table$V1)

tval_table <- read.table('TVAL_FILE')
tvals <- setNames(pval_table$V2,pval_table$V1)


pathways <- loadGSC('PATHWAY_FILE')

# Run gene-set analysis:
gsares_1 <- runGSA(geneLevelStats=original_pval,
                 directions = directions, 
                 gsc = pathways, geneSetStat="mean")

gsares_2 <- runGSA(geneLevelStats=original_pval,
                 directions = directions,
                 gsc = pathways, geneSetStat="median")

gsares_3 <- runGSA(geneLevelStats=original_pval,
                   directions = directions,
                   gsc = pathways, geneSetStat="sum")

gsares_5 <- runGSA(geneLevelStats=original_pval,
                   directions = directions,
                   gsc = pathways, geneSetStat="fisher")

gsares_6 <- runGSA(geneLevelStats=original_pval,
                   directions = directions,
                   gsc = pathways, geneSetStat="stouffer")

gsares_7 <- runGSA(geneLevelStats=original_pval,
                   directions = directions,
                   gsc = pathways, geneSetStat="reporter")


gsares_8 <- runGSA(geneLevelStats=original_pval,
                   directions = directions, nPerm = 1000,
                   gsc = pathways, geneSetStat="wilcoxon")

gsares_9 <- runGSA(geneLevelStats=original_pval,
                   directions = directions,
                   gsc = pathways, geneSetStat="tailStrength")

# t-val statts needed

gsares_4 <- runGSA(geneLevelStats=tvals,
                   gsc = pathways, geneSetStat="maxmean")

gsares_10 <- runGSA(geneLevelStats=tvals,
                   gsc = pathways, geneSetStat="gsea")


gsares_11 <- runGSA(geneLevelStats=tvals,
                    gsc = pathways, geneSetStat="fgsea")

gsares_12 <- runGSA(geneLevelStats=tvals,
                    gsc = pathways, geneSetStat="page")


combres <- consensusScores(resList=list(gsares_1, gsares_2, gsares_3, gsares_4, gsares_5, gsares_8,
                                        gsares_6, gsares_7, gsares_9, gsares_10, gsares_11, gsares_12),
                           class='distinct', direction='up')

write.table(combres$rankMat, file='OUTFILE_RANK', sep='\t', quote = FALSE)
write.table(combres$pMat, file='OUTFILE_PVAL', sep='\t', quote = FALSE)

