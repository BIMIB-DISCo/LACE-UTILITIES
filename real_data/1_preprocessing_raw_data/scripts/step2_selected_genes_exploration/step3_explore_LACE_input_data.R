# set working directory
baseDir = "/data_directory/BimiB/share/LACE/MELANOMA/"
setwd(baseDir)

depth_minimum = 3 # minimum depth to set values to NA
missing_values_max = 0.35 # maximum number of considered missing data per gene

# load data
snpMut_filt_freq <- readRDS(file=paste0("snpMut_filt_freq.rds"))
mutations = as.matrix(read.table("final_data_mut.txt"))
depth = as.matrix(read.table("final_data_depth.txt"))

# set NA values
depth = depth[,colnames(mutations)]
mutations[which(mutations>1,arr.ind=TRUE)] = 1
mutations[which(depth<=depth_minimum,arr.ind=TRUE)] = NA

mycellsdata = list()
mycellsdata[["T1_before_treatment"]] = t(mutations[,sort(unique(snpMut_filt_freq$scID[which(snpMut_filt_freq$Time=="before treatment")]))])
mycellsdata[["T2_4_days_treatment"]] = t(mutations[,sort(unique(snpMut_filt_freq$scID[which(snpMut_filt_freq$Time=="4d on treatment")]))])
mycellsdata[["T3_28_days_treatment"]] = t(mutations[,sort(unique(snpMut_filt_freq$scID[which(snpMut_filt_freq$Time=="28d on treatment")]))])
mycellsdata[["T4_57_days_treatment"]] = t(mutations[,sort(unique(snpMut_filt_freq$scID[which(snpMut_filt_freq$Time=="57d on treatment")]))])

# evaluate number of missing data per time point
t1 = NULL
t2 = NULL
t3 = NULL
t4 = NULL
for(i in 1:ncol(mycellsdata[["T1_before_treatment"]])) {
    t1 = c(t1,length(which(is.na(mycellsdata[["T1_before_treatment"]][,i])))/nrow(mycellsdata[["T1_before_treatment"]]))
    t2 = c(t2,length(which(is.na(mycellsdata[["T2_4_days_treatment"]][,i])))/nrow(mycellsdata[["T2_4_days_treatment"]]))
    t3 = c(t3,length(which(is.na(mycellsdata[["T3_28_days_treatment"]][,i])))/nrow(mycellsdata[["T3_28_days_treatment"]]))
    t4 = c(t4,length(which(is.na(mycellsdata[["T4_57_days_treatment"]][,i])))/nrow(mycellsdata[["T4_57_days_treatment"]]))
}
valid_genes = sort(unique(colnames(mycellsdata[["T1_before_treatment"]])[which(t1<=missing_values_max&t2<=missing_values_max&t3<=missing_values_max&t4<=missing_values_max)]))
valid = NULL
for(i in valid_genes) {
    valid = c(valid,unique(snpMut_filt_freq$Gene[grep(i,snpMut_filt_freq$UIDsnp)]))
}
valid = paste0(valid,"_",valid_genes)

# get names of valid genes
valid_genes_names = NULL
for(i in valid) {
    valid_genes_names = c(valid_genes_names,strsplit(i,split="_")[[1]][[1]])
}
print(length(valid_genes_names))
print(length(unique(valid_genes_names)))
print(sort(valid_genes_names[which(duplicated(valid_genes_names))]))

# make final list of candidate selected variants
snpMut_filt_freq_reduced = unique(snpMut_filt_freq[which(snpMut_filt_freq$Gene%in%valid_genes_names),c("Gene","Chr","PosStart","PosEnd","REF","ALT")])
snpMut_filt_freq_reduced = snpMut_filt_freq_reduced[order(snpMut_filt_freq_reduced[,1],snpMut_filt_freq_reduced[,2],snpMut_filt_freq_reduced[,3],snpMut_filt_freq_reduced[,4],snpMut_filt_freq_reduced[,5],snpMut_filt_freq_reduced[,6]),]
write.table(snpMut_filt_freq_reduced,file="snpMut_filt_freq_reduced.txt",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
