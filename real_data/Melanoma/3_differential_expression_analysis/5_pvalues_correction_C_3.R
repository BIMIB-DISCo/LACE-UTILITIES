# set the working directory
setwd("/data_directory/BimiB/share/LACE/MELANOMA/real_data/3_differential_expression_analysis")

# load the results
load(file="anova_C_3/pvalues_all_timepoints.RData")
pvalues_all_timepoints = pvalues
load(file="anova_C_3/pvalues.RData")

# select only valid features
pvalues[[1]] = pvalues[[1]]
pvalues[[2]] = pvalues[[2]]
pvalues[[3]] = pvalues[[3]]
pvalues[[4]] = pvalues[[4]]
pvalues_all_timepoints = pvalues_all_timepoints

# pvalue threshold
thr = 0.20

# set the seed
set.seed(8765943)

# pvalue correction
pvalues_original = pvalues
pvalues_original[[5]] = pvalues_all_timepoints
pvalues_bonferroni = list()
pvalues_bonferroni[[1]] = p.adjust(as.numeric(pvalues_original[[1]]),method="bonferroni")
pvalues_bonferroni[[2]] = p.adjust(as.numeric(pvalues_original[[2]]),method="bonferroni")
pvalues_bonferroni[[3]] = p.adjust(as.numeric(pvalues_original[[3]]),method="bonferroni")
pvalues_bonferroni[[4]] = p.adjust(as.numeric(pvalues_original[[4]]),method="bonferroni")
pvalues_bonferroni[[5]] = p.adjust(as.numeric(pvalues_original[[5]]),method="bonferroni")
pvalues_fdr = list()
pvalues_fdr[[1]] = p.adjust(as.numeric(pvalues_original[[1]]),method="fdr")
pvalues_fdr[[2]] = p.adjust(as.numeric(pvalues_original[[2]]),method="fdr")
pvalues_fdr[[3]] = p.adjust(as.numeric(pvalues_original[[3]]),method="fdr")
pvalues_fdr[[4]] = p.adjust(as.numeric(pvalues_original[[4]]),method="fdr")
pvalues_fdr[[5]] = p.adjust(as.numeric(pvalues_original[[5]]),method="fdr")

# select significant features
pvalues_values = list()
pvalues_values[["original"]] = pvalues_original
pvalues_values[["bonferroni"]] = pvalues_bonferroni
pvalues_values[["fdr"]] = pvalues_fdr
pvalues_significant_bonferroni = list()
pvalues_significant_bonferroni = pvalues_bonferroni
for(i in 1:5) {
    curr_gene = names(pvalues_original[[i]])[which(pvalues_significant_bonferroni[[i]]<thr)]
    curr_val = pvalues_significant_bonferroni[[i]][which(pvalues_significant_bonferroni[[i]]<thr)]
    names(curr_val) = curr_gene
    pvalues_significant_bonferroni[[i]] = curr_val
}
pvalues_significant_fdr = list()
pvalues_significant_fdr = pvalues_fdr
for(i in 1:5) {
    curr_gene = names(pvalues_original[[i]])[which(pvalues_significant_fdr[[i]]<thr)]
    curr_val = pvalues_significant_fdr[[i]][which(pvalues_significant_fdr[[i]]<thr)]
    names(curr_val) = curr_gene
    pvalues_significant_fdr[[i]] = curr_val
}
pvalues_significant = list()
pvalues_significant[["bonferroni"]] = pvalues_significant_bonferroni
pvalues_significant[["fdr"]] = pvalues_significant_fdr

# save the results
save(pvalues_values,file="results_C_3/pvalues_values.RData")
save(pvalues_significant,file="results_C_3/pvalues_significant.RData")
