# set the working directory
setwd("/data_directory/BimiB/share/LACE/MELANOMA/real_data/3_differential_expression_analysis")

# read raw data
inference_clones_mapping = read.table(file="raw_data/Inference_clones.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
inference_clones_mapping = inference_clones_mapping[order(inference_clones_mapping[,"Run"],inference_clones_mapping[,"GSM"],inference_clones_mapping[,"ExpID"],inference_clones_mapping[,"Age"],inference_clones_mapping[,"Clone"]),]
RNA_data = read.table(file="raw_data/GSE116237_scRNAseq_expressionMatrix.txt",header=TRUE,sep=",",row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
RNA_data = as.matrix(t(RNA_data))
RNA_data = RNA_data[sort(rownames(RNA_data)),]
RNA_data = RNA_data[,sort(colnames(RNA_data))]
RNA_data = RNA_data[inference_clones_mapping$ExpID,]

# save results
save(inference_clones_mapping,file="raw_data/inference_clones_mapping.RData")
save(RNA_data,file="raw_data/RNA_data.RData")
