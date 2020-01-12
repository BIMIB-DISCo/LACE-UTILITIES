# set the working directory
setwd("/data_directory/BimiB/share/LACE/MELANOMA/real_data/3_differential_expression_analysis")

library("ggplot2")

# read data
load(file="raw_data/inference_clones_mapping.RData")
load(file="raw_data/RNA_data.RData")
RNA_data_raw = RNA_data
load(file="raw_data/RNA_data_final.RData")
RNA_data_normalized = RNA_data_final

# data for the analysis
RNA_data_raw = RNA_data_raw[,"ENSG00000185686",drop=FALSE]
RNA_data_normalized = RNA_data_normalized[,"ENSG00000185686",drop=FALSE]
RNA_data_normalized = RNA_data_normalized[sort(rownames(RNA_data_normalized)),,drop=FALSE]
RNA_data_raw = RNA_data_raw[rownames(RNA_data_normalized),,drop=FALSE]
RNA_data_normalized = RNA_data_normalized[,sort(colnames(RNA_data_normalized)),drop=FALSE]
RNA_data_raw = RNA_data_raw[,colnames(RNA_data_normalized),drop=FALSE]

# read cluster assignments
clones = inference_clones_mapping$Clone
clones[which(clones%in%c(2,3,4,6))] = "PRAME_mut"
clones[which(clones%in%c(1,5))] = "PRAME_wt"
g0_PRAME_mut = rownames(RNA_data_raw)[which(rownames(RNA_data_raw)%in%inference_clones_mapping$ExpID[which(clones=="PRAME_mut")])]
g1_PRAME_wt = rownames(RNA_data_raw)[which(rownames(RNA_data_raw)%in%inference_clones_mapping$ExpID[which(clones=="PRAME_wt")])]

# boxplot for PRAME in the two groups

g0_PRAME_mut_exp = as.numeric(RNA_data_normalized[g0_PRAME_mut,"ENSG00000185686"])
g1_PRAME_wt_exp = as.numeric(RNA_data_normalized[g1_PRAME_wt,"ENSG00000185686"])
VALUE = c(g0_PRAME_mut_exp,g1_PRAME_wt_exp)
GROUP = c(rep("PRAME mutated",length(g0_PRAME_mut_exp)),rep("PRAME wt",length(g1_PRAME_wt_exp)))
p <- ggplot(data.frame(VALUE_LOG2=log(VALUE,base=2),GROUP=GROUP), aes(x=GROUP,y=VALUE_LOG2,fill=GROUP)) + geom_boxplot()
p
print(mean(g0_PRAME_mut_exp))
print(mean(g1_PRAME_wt_exp))
print(median(g0_PRAME_mut_exp))
print(median(g1_PRAME_wt_exp))
print(t.test(log(g0_PRAME_mut_exp+1,base=2),log(g1_PRAME_wt_exp+1,base=2)))
