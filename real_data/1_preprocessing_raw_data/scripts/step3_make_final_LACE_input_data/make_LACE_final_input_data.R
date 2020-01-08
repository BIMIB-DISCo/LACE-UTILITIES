# set working directory
baseDir = "/data_directory/BimiB/share/LACE/MELANOMA/"
setwd(baseDir)

# load required libraries
library("TRONCO")

# list of manually verified mutations
verified_genes = c("ARPC2","CCT8","COL1A2","CYCS","HNRNPC","PCBP1","PRAME","RPL5")

depth_minimum = 3 # minimum depth to set values to NA
minumum_median_total = 10 # minimum median depth for total reads
minumum_median_mutation = 4 # minimum median depth for reads supporting mutations

# load data
cells_aggregate_info <- readRDS(file=paste0("cells_aggregate_info.rds"))
snpMut_filt_freq <- readRDS(file=paste0("snpMut_filt_freq.rds"))
depth = as.matrix(read.table("final_data_depth.txt"))

# select only verified mutations
snpMut_filt_freq = snpMut_filt_freq[which(snpMut_filt_freq$Gene%in%verified_genes),c("scID","Time","Gene","Chr","PosStart","PosEnd","REF","ALT","MutType","depth","Allele_Ratio")]
snpMut_filt_freq = snpMut_filt_freq[order(snpMut_filt_freq[,3],snpMut_filt_freq[,4],snpMut_filt_freq[,5],snpMut_filt_freq[,6],snpMut_filt_freq[,7],snpMut_filt_freq[,8],snpMut_filt_freq[,9],snpMut_filt_freq[,1],snpMut_filt_freq[,2],snpMut_filt_freq[,10],snpMut_filt_freq[,11]),]

# compute frequency of each mutation
distinct_mutations = unique(snpMut_filt_freq[,c("Gene","Chr","PosStart","PosEnd","REF","ALT","ALT","ALT","ALT","ALT","ALT","ALT")])
colnames(distinct_mutations)[7] = "FreqT1"
colnames(distinct_mutations)[8] = "FreqT2"
colnames(distinct_mutations)[9] = "FreqT3"
colnames(distinct_mutations)[10] = "FreqT4"
colnames(distinct_mutations)[11] = "MedianDepth"
colnames(distinct_mutations)[12] = "MedianDepthMut"
cells_timepoints = unique(snpMut_filt_freq[,c("scID","Time")])
distinct_mutations$FreqT1 = NA
distinct_mutations$FreqT2 = NA
distinct_mutations$FreqT3 = NA
distinct_mutations$FreqT4 = NA
distinct_mutations$MedianDepth = NA
distinct_mutations$MedianDepthMut = NA
t1 = as.numeric(table(cells_timepoints$Time)["before treatment"])
t2 = as.numeric(table(cells_timepoints$Time)["4d on treatment"])
t3 = as.numeric(table(cells_timepoints$Time)["28d on treatment"])
t4 = as.numeric(table(cells_timepoints$Time)["57d on treatment"])
for(i in 1:nrow(distinct_mutations)) {
    curr = snpMut_filt_freq[which(snpMut_filt_freq$Gene%in%distinct_mutations[i,"Gene"]&snpMut_filt_freq$Chr%in%distinct_mutations[i,"Chr"]&snpMut_filt_freq$PosStart%in%distinct_mutations[i,"PosStart"]&snpMut_filt_freq$PosEnd%in%distinct_mutations[i,"PosEnd"]&snpMut_filt_freq$REF%in%distinct_mutations[i,"REF"]&snpMut_filt_freq$ALT%in%distinct_mutations[i,"ALT"]),]
    distinct_mutations[i,"FreqT1"] = as.numeric(table(curr$Time)["before treatment"]) / t1
    distinct_mutations[i,"FreqT2"] = as.numeric(table(curr$Time)["4d on treatment"]) / t2
    distinct_mutations[i,"FreqT3"] = as.numeric(table(curr$Time)["28d on treatment"]) / t3
    distinct_mutations[i,"FreqT4"] = as.numeric(table(curr$Time)["57d on treatment"]) / t4
    distinct_mutations[i,"MedianDepth"] = as.numeric(median(curr$depth))
    distinct_mutations[i,"MedianDepthMut"] = as.numeric(median(curr$Allele_Ratio))
}
distinct_mutations = distinct_mutations[-c(3,5,6,8),]
distinct_mutations = distinct_mutations[-sort(unique(c(which(distinct_mutations$MedianDepth<minumum_median_total),which(distinct_mutations$MedianDepthMut<minumum_median_mutation)))),]
rownames(distinct_mutations) = 1:nrow(distinct_mutations)
valid_distinct_mutations = distinct_mutations
valid_distinct_mutations_values = NULL
for(i in 1:nrow(valid_distinct_mutations)) {
    valid_distinct_mutations_values = c(valid_distinct_mutations_values,paste0(valid_distinct_mutations[i,c("Chr","PosStart")],collapse="_"))
}
save(valid_distinct_mutations,file="valid_distinct_mutations.RData")

# make final mutations data structures
load("valid_clones_mapping.RData")
mutations = array(0,c(length(unique(valid_clones_mapping$Run)),6))
rownames(mutations) = sort(unique(valid_clones_mapping$Run))
colnames(mutations) = paste0(valid_distinct_mutations$Gene,"_",valid_distinct_mutations_values,"_",valid_distinct_mutations$REF,"_",valid_distinct_mutations$ALT)
for(i in 1:nrow(valid_distinct_mutations)) {
    curr_gene = valid_distinct_mutations[i,"Gene"]
    curr_chr = valid_distinct_mutations[i,"Chr"]
    curr_start = valid_distinct_mutations[i,"PosStart"]
    curr_end = valid_distinct_mutations[i,"PosEnd"]
    curr_ref = valid_distinct_mutations[i,"REF"]
    curr_alt = valid_distinct_mutations[i,"ALT"]
    curr_mutant_cells = which(cells_aggregate_info$Gene==curr_gene&cells_aggregate_info$Chr==curr_chr&cells_aggregate_info$PosStart==curr_start&cells_aggregate_info$PosEnd==curr_end&cells_aggregate_info$REF==curr_ref&cells_aggregate_info$ALT==curr_alt)
    curr_mutant_cells = cells_aggregate_info$scID[curr_mutant_cells]
    mutations[curr_mutant_cells[which(curr_mutant_cells%in%rownames(mutations))],i] = 1
}

# set NA values
depth = t(depth)
depth = depth[,valid_distinct_mutations_values]
depth = depth[rownames(mutations),]
colnames(depth) = colnames(mutations)
mutations[which(depth<=depth_minimum,arr.ind=TRUE)] = NA # missing values rate equals to 573/3792, that is approx 15%

# make final data
cells_aggregate_info = cells_aggregate_info[which(cells_aggregate_info$scID%in%rownames(mutations)),]
mycellsdata = list()
mycellsdata[["T1_before_treatment"]] = mutations[sort(unique(cells_aggregate_info$scID[which(cells_aggregate_info$Time=="before treatment")])),]
mycellsdata[["T2_4_days_treatment"]] = mutations[sort(unique(cells_aggregate_info$scID[which(cells_aggregate_info$Time=="4d on treatment")])),]
mycellsdata[["T3_28_days_treatment"]] = mutations[sort(unique(cells_aggregate_info$scID[which(cells_aggregate_info$Time=="28d on treatment")])),]
mycellsdata[["T4_57_days_treatment"]] = mutations[sort(unique(cells_aggregate_info$scID[which(cells_aggregate_info$Time=="57d on treatment")])),]
D = mycellsdata
save(D,file="D.RData")

# make Oncoprint
clusters = array(NA,c((nrow(D[[1]])+nrow(D[[2]])+nrow(D[[3]])+nrow(D[[4]])),1))
clusters[1:nrow(D[[1]]),1] = "T1_before_treatment"
clusters[(nrow(D[[1]])+1):(nrow(D[[1]])+nrow(D[[2]])),1] = "T2_4_days_treatment"
clusters[(nrow(D[[1]])+nrow(D[[2]])+1):(nrow(D[[1]])+nrow(D[[2]])+nrow(D[[3]])),1] = "T3_28_days_treatment"
clusters[(nrow(D[[1]])+nrow(D[[2]])+nrow(D[[3]])+1):(nrow(D[[1]])+nrow(D[[2]])+nrow(D[[3]])+nrow(D[[4]])),1] = "T4_57_days_treatment"
rownames(clusters) = c(rownames(D[[1]]),rownames(D[[2]]),rownames(D[[3]]),rownames(D[[4]]))
data = rbind(D[[1]],D[[2]],D[[3]],D[[4]])
data[which(is.na(data))] = 0
data = import.genotypes(data)
data = annotate.stages(data,clusters)
oncoprint(data,excl.sort=FALSE,group.by.stage=TRUE)
