# set the working directory
setwd("/data_directory/BimiB/share/LACE/MELANOMA/real_data/3_differential_expression_analysis")

# load the required libraries
library("biomaRt")

# load the results
load(file="raw_data/RNA_data.RData")

# convert ensample names to genes id
features = sort(unique(colnames(RNA_data)))
mart = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl")
G_list = getBM(filters="ensembl_gene_id",attributes= c("ensembl_gene_id","hgnc_symbol","gene_biotype"),values=features,mart=mart)
genes_ensample_mapping = G_list[order(G_list[,"ensembl_gene_id"],G_list[,"hgnc_symbol"],G_list[,"gene_biotype"]),]
genes_ensample_mapping = genes_ensample_mapping[-which(genes_ensample_mapping$gene_biotype%in%names(table(genes_ensample_mapping$gene_biotype))[grep("pseudogene",names(table(genes_ensample_mapping$gene_biotype)))]),] # remove pseudogenes
genes_ensample_mapping = genes_ensample_mapping[-which(genes_ensample_mapping$gene_biotype%in%c("Mt_rRNA","Mt_tRNA")),] # remove mitochondrial genes
genes_ensample_mapping = genes_ensample_mapping[-grep("MT-",genes_ensample_mapping$hgnc_symbol),] # remove mitochondrial genes
genes_ensample_mapping = genes_ensample_mapping[-which(genes_ensample_mapping$hgnc_symbol==""),] # remove features with no hgnc symbol
rownames(genes_ensample_mapping) = 1:nrow(genes_ensample_mapping)
save(genes_ensample_mapping,file="raw_data/genes_ensample_mapping.RData")

# build final RNA_data
RNA_data = RNA_data[,genes_ensample_mapping$ensembl_gene_id]
RNA_data = RNA_data[sort(rownames(RNA_data)),]
RNA_data = RNA_data[,sort(colnames(RNA_data))]

# normalization by library size
RNA_data = RNA_data / rowSums(RNA_data)
RNA_data = RNA_data * 10^6
RNA_data_final = RNA_data
save(RNA_data_final,file="raw_data/RNA_data_final.RData")
