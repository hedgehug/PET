# this script is used to run Enrichment browser combined and individual base method results, by default, combine all base methods' reults
# reference from: https://bioconductor.org/packages/release/bioc/html/EnrichmentBrowser.html
library(EnrichmentBrowser)

gmt.file<-file.path('PATHWAY_FILE') 
hsa.gs<-getGenesets(gmt.file)

data.dir <- 'Enrichment_browser//data/'
exprs.file<-file.path(data.dir,"TARGET.exprs.tab") 
cdat.file<-file.path(data.dir,"TARGET.colData.tab") 
rdat.file<-file.path(data.dir,"TARGET.rowData.tab") 
se<-readSE(exprs.file,cdat.file,rdat.file)


se<-deAna(se,de.method="DESeq2", filter.by.expr=FALSE)

method_idx <- 1
res_list <- list()

# run individual method in enrichment browser and write the result to file
while (method_idx<=length(EnrichmentBrowser::sbeaMethods())) {
  if (method_idx == 5){
    method_idx <- method_idx + 1
    next
  }
  sbea.res<-sbea(method=EnrichmentBrowser::sbeaMethods()[method_idx],se=se,gs=hsa.gs, alpha = 1, out.file='individual_method_output_file_name' )
  res_list[[length(res_list)+1]] <- sbea.res
  # [] <- append(res_list, sbea.res)
  method_idx <- method_idx + 1
}

# combine the results and write to file
comb.res <- combResults(res_list)
res <- gsRanking(comb.res, signif.only = FALSE)
write.table(res, 'output_file_name', sep = '\t', quote = FALSE)
