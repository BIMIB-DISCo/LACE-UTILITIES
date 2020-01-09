# set the working directory
setwd("/data_directory/BimiB/share/LACE/MELANOMA/real_data/3_differential_expression_analysis")

# read data
load(file="raw_data/genes_ensample_mapping.RData")
load(file="results_C_3/pvalues_significant.RData")

summary_results = pvalues_significant[["fdr"]]
names(summary_results) = c("before treatment","4d on treatment","28d on treatment","57d on treatment","all timepoints")

res2 = genes_ensample_mapping[which(genes_ensample_mapping$ensembl_gene_id%in%names(summary_results[["4d on treatment"]])),]
values = pvalues_significant[["fdr"]][[2]]
pvalues = NULL
res = res2
for(i in 1:nrow(res)) {
    pvalues = c(pvalues,as.numeric(values[res[i,"ensembl_gene_id"]]))
}
res = cbind(res,pvalues)
res2 = res[sort.int(res$pvalues,index.return=TRUE)$ix,]
rownames(res2) = 1:nrow(res2)

res3 = genes_ensample_mapping[which(genes_ensample_mapping$ensembl_gene_id%in%names(summary_results[["28d on treatment"]])),]
values = pvalues_significant[["fdr"]][[3]]
pvalues = NULL
res = res3
for(i in 1:nrow(res)) {
    pvalues = c(pvalues,as.numeric(values[res[i,"ensembl_gene_id"]]))
}
res = cbind(res,pvalues)
res3 = res[sort.int(res$pvalues,index.return=TRUE)$ix,]
rownames(res3) = 1:nrow(res3)

res5 = genes_ensample_mapping[which(genes_ensample_mapping$ensembl_gene_id%in%names(summary_results[["all timepoints"]])),]
values = pvalues_significant[["fdr"]][[5]]
pvalues = NULL
res = res5
for(i in 1:nrow(res)) {
    pvalues = c(pvalues,as.numeric(values[res[i,"ensembl_gene_id"]]))
}
res = cbind(res,pvalues)
res5 = res[sort.int(res$pvalues,index.return=TRUE)$ix,]
rownames(res5) = 1:nrow(res5)

summary_results[["4d on treatment"]] = res2
summary_results[["28d on treatment"]] = res3
summary_results[["all timepoints"]] = res5

save(summary_results,file="summary/summary_results.RData")
