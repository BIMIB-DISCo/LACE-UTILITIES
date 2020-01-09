# set the working directory
setwd("/data_directory/BimiB/share/LACE/MELANOMA/real_data/3_differential_expression_analysis")

library("ggplot2")

# read data
load(file="raw_data/inference_clones_mapping.RData")
load(file="raw_data/RNA_data.RData")
RNA_data_raw = RNA_data
load(file="raw_data/RNA_data_final.RData")
RNA_data_normalized = RNA_data_final
load(file="summary/summary_results.RData")

# data for the analysis
RNA_data_raw = RNA_data_raw[,summary_results[["28d on treatment"]][1:156,"ensembl_gene_id"]]
RNA_data_normalized = RNA_data_normalized[,summary_results[["28d on treatment"]][1:156,"ensembl_gene_id"]]
RNA_data_normalized = RNA_data_normalized[sort(rownames(RNA_data_normalized)),]
RNA_data_raw = RNA_data_raw[rownames(RNA_data_normalized),]
RNA_data_normalized = RNA_data_normalized[,sort(colnames(RNA_data_normalized))]
RNA_data_raw = RNA_data_raw[,colnames(RNA_data_normalized)]

# select only cells in time point 3
RNA_data_normalized = RNA_data_normalized[rownames(RNA_data_normalized)[which(rownames(RNA_data_normalized)%in%inference_clones_mapping$ExpID[which(inference_clones_mapping$Age=="28d on treatment")])],]
RNA_data_raw = RNA_data_raw[rownames(RNA_data_raw)[which(rownames(RNA_data_raw)%in%inference_clones_mapping$ExpID[which(inference_clones_mapping$Age=="28d on treatment")])],]

# read cluster assignments
clones = inference_clones_mapping$Clone
clones[which(clones%in%c(2,3,4,6))] = "PRAME_mut"
clones[which(clones%in%c(1,5))] = "PRAME_wt"
g0_PRAME_mut = rownames(RNA_data_raw)[which(rownames(RNA_data_raw)%in%inference_clones_mapping$ExpID[which(clones=="PRAME_mut")])]
g1_PRAME_wt = rownames(RNA_data_raw)[which(rownames(RNA_data_raw)%in%inference_clones_mapping$ExpID[which(clones=="PRAME_wt")])]

# fold change option 1: using raw data, log of differences in mean
fold = rep(NA,ncol(RNA_data_raw))
names(fold) = colnames(RNA_data_raw)
data = (RNA_data_raw + 1)
for(i in 1:ncol(data)) {
    g0_PRAME_mut_val = mean(data[g0_PRAME_mut,i])
    g1_PRAME_wt_val = mean(data[g1_PRAME_wt,i])
    res = log((g0_PRAME_mut_val/g1_PRAME_wt_val),base=2)
    fold[i] = res
}
fold1 = fold

# fold change option 2: using raw data, differences in means of logs
fold = rep(NA,ncol(RNA_data_raw))
names(fold) = colnames(RNA_data_raw)
data = log((RNA_data_raw+1),base=2)
for(i in 1:ncol(data)) {
    g0_PRAME_mut_val = mean(data[g0_PRAME_mut,i])
    g1_PRAME_wt_val = mean(data[g1_PRAME_wt,i])
    res = log((g0_PRAME_mut_val/g1_PRAME_wt_val),base=2)
    fold[i] = res
}
fold2 = fold

# fold change option 3: using normalized data, log of differences in mean
fold = rep(NA,ncol(RNA_data_normalized))
names(fold) = colnames(RNA_data_normalized)
data = (RNA_data_normalized + 1)
for(i in 1:ncol(data)) {
    g0_PRAME_mut_val = mean(data[g0_PRAME_mut,i])
    g1_PRAME_wt_val = mean(data[g1_PRAME_wt,i])
    res = log((g0_PRAME_mut_val/g1_PRAME_wt_val),base=2)
    fold[i] = res
}
fold3 = fold

# fold change option 4: using normalized data, differences in means of logs
fold = rep(NA,ncol(RNA_data_normalized))
names(fold) = colnames(RNA_data_normalized)
data = log((RNA_data_normalized+1),base=2)
for(i in 1:ncol(data)) {
    g0_PRAME_mut_val = mean(data[g0_PRAME_mut,i])
    g1_PRAME_wt_val = mean(data[g1_PRAME_wt,i])
    res = log((g0_PRAME_mut_val/g1_PRAME_wt_val),base=2)
    fold[i] = res
}
fold4 = fold

# save final results
significant_features_fold_change = cbind(summary_results[["28d on treatment"]][1:156,],summary_results[["28d on treatment"]][1:156,"pvalues"],summary_results[["28d on treatment"]][1:156,"pvalues"],summary_results[["28d on treatment"]][1:156,"pvalues"],summary_results[["28d on treatment"]][1:156,"pvalues"])
colnames(significant_features_fold_change)[5:8] = c("fold_change_v1","fold_change_v2","fold_change_v3","fold_change_v4")
for(i in 1:nrow(significant_features_fold_change)) {
    significant_features_fold_change[i,"fold_change_v1"] = fold1[[significant_features_fold_change[i,"ensembl_gene_id"]]]
    significant_features_fold_change[i,"fold_change_v2"] = fold2[[significant_features_fold_change[i,"ensembl_gene_id"]]]
    significant_features_fold_change[i,"fold_change_v3"] = fold3[[significant_features_fold_change[i,"ensembl_gene_id"]]]
    significant_features_fold_change[i,"fold_change_v4"] = fold4[[significant_features_fold_change[i,"ensembl_gene_id"]]]
}

significant_features_fold_change_T3 = significant_features_fold_change
save(significant_features_fold_change_T3,file="summary/significant_features_fold_change_T3.RData")

# genes overexpressed in PRAME_mut (fold change >3)
up_PRAME_mut = significant_features_fold_change$hgnc_symbol[which(significant_features_fold_change$fold_change_v3>(3))]
up_PRAME_mut = significant_features_fold_change[which(significant_features_fold_change$hgnc_symbol%in%up_PRAME_mut),]
up_PRAME_mut = up_PRAME_mut[which(up_PRAME_mut$pvalues<0.10),]
rownames(up_PRAME_mut) = 1:nrow(up_PRAME_mut)
save(up_PRAME_mut,file="summary/up_PRAME_mut.RData")

# genes overexpressed in PRAME_wt (fold change < -3)
down_PRAME_mut = significant_features_fold_change$hgnc_symbol[which(significant_features_fold_change$fold_change_v3<(-3))]
down_PRAME_mut = significant_features_fold_change[which(significant_features_fold_change$hgnc_symbol%in%down_PRAME_mut),]
down_PRAME_mut = down_PRAME_mut[which(down_PRAME_mut$pvalues<0.10),]
rownames(down_PRAME_mut) = 1:nrow(down_PRAME_mut)
save(down_PRAME_mut,file="summary/down_PRAME_mut.RData")

# boxplot for genes overexpressed in PRAME_mut

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000156515"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000156515"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000162616"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000162616"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000063241"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000063241"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000144354"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000144354"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000151092"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000151092"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))

# boxplot for genes overexpressed in PRAME_wt

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000023516"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000023516"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000111837"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000111837"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000236104"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000236104"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))
