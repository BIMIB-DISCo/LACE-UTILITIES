# set the working directory
setwd("/data_directory/BimiB/share/LACE/MELANOMA/real_data/3_differential_expression_analysis")

# set the seed
set.seed(7654043)

# read the data
load(file="raw_data/inference_clones_mapping.RData")
load(file="raw_data/RNA_data_final.RData")
RNA_data = RNA_data_final

# log2 transformation
RNA_data = RNA_data + 1
RNA_data = log(RNA_data,base=2)

# perform anova for each gene for expression data
RNA_data = RNA_data[inference_clones_mapping$ExpID,]
clones = inference_clones_mapping$Clone
clones[which(clones%in%c(2,3,4,6))] = 3
clones[which(clones==5)] = 2
data = cbind(RNA_data,clones)
colnames(data)[1:(ncol(data)-1)] = paste0("Feature_",1:(ncol(data)-1))
colnames(data)[ncol(data)] = "Cluster"
data = as.data.frame(data)
pvalues = rep(NA,(ncol(data)-1))
names(pvalues) = colnames(RNA_data)
cat("Performing ANOVA for expression data...\n")
all_pvalues = list()
cont = 0
cat(cont,"\n")
for(time_point in c("before treatment","4d on treatment","28d on treatment","57d on treatment")) {
    cont = cont + 1
    curr_data = data[inference_clones_mapping$ExpID[which(inference_clones_mapping$Age==time_point)],]
    curr_pvalues = pvalues
    for(i in 1:(ncol(curr_data)-1)) {
        curr_aov = aov(as.formula(paste0("Cluster~Feature_",i)),data=curr_data)
        curr_pvalues[i] = summary(curr_aov)[[1]][["Pr(>F)"]][1]
    }
    all_pvalues[[cont]] = curr_pvalues
    cat(cont/4,"\n")
}
pvalues = all_pvalues

# save the results
save(pvalues,file="anova_C_3/pvalues.RData")
