##################################################################
# GET DEPTH FOR EACH SELECTED MUTATION FROM BAMs
##################################################################

# set working directory
baseDir = "/data_directory/BimiB/share/LACE/MELANOMA/"
setwd(baseDir)

# load required libraries
library("foreach")
library("doParallel")

# load data
snpMut_filt_freq <- readRDS(file=paste0("snpMut_filt_freq.rds"))

listMut <- strsplit(snpMut_filt_freq$UIDsnp, "_")
snpMut_filt_freq$Mut_pos <- unlist(lapply(listMut, function(m) paste0(m[2], "_", m[3])))

listMutFormatted <- unlist(lapply(unique(strsplit(snpMut_filt_freq$Mut_pos, "_")), function(m) paste0(m[1], ":", m[2], "-", m[2])))

BAMlist <- paste(dir(path = paste0(baseDir, "picard"), pattern = "s4.bam$", full.names = TRUE), collapse=" ")

numCores <- detectCores()
registerDoParallel(numCores)

res <- foreach(i = 1:length(listMutFormatted), .combine=c) %dopar% {
    system(paste0("samtools depth -a -Q 20 -r ",listMutFormatted[i]," ", BAMlist), intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
}

BAM = gsub(x = paste(unlist(strsplit(x = gsub(x=BAMlist, pattern=paste0(baseDir, "picard/"), ""), split=".bam")), collapse = "\t"), " ","")
header = paste0("CHR\tPOS\t",BAM)


##################################################################
# MAKE DEPTH MATRIX
##################################################################

depth_tbl <- read.table(text = paste0(header, "\n", paste(res, collapse='\n')), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(depth_tbl) <- paste(depth_tbl$CHR, depth_tbl$POS, sep = "_")
depth_tbl <- depth_tbl[,-c(1,2)]

snpMut_filt_freq$Mut_pos <- paste(snpMut_filt_freq$Chr, snpMut_filt_freq$PosStart, sep = "_")

mut_tbl <- as.data.frame.matrix(table(snpMut_filt_freq$Mut_pos, snpMut_filt_freq$scID))

cellInfo <- readRDS("data_info.rds")

for(i in 1:ncol(depth_tbl)) {
    run <- colnames(depth_tbl)[i]
    colnames(depth_tbl)[i] <- as.character(cellInfo$Run[cellInfo$Run == gsub(".bbq_s4","",run)])
}

depth_tbl <- depth_tbl[ order(row.names(depth_tbl)), order(colnames(depth_tbl))]
mut_tbl <- mut_tbl[ order(row.names(mut_tbl)), order(colnames(mut_tbl))]

# write final results to text files
write.table(x = depth_tbl, file = paste0(baseDir, "final_data_depth.txt"), sep = '\t', row.names = TRUE, col.names = TRUE)
write.table(x = mut_tbl, file = paste0(baseDir, "final_data_mut.txt"), sep = '\t', row.names = TRUE, col.names = TRUE)
