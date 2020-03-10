# Draw B
draw.B <- function( B, mut_label = colnames(B)[-1] ) {
    
    Broot <- Node$new('r')
    Broot$mut <- B[1,]
    
    nClone <- nrow(B)
    Clones <- list(Broot)
    
    for(rP in 1:(nrow(B)-1)) {
        
        for(rC in ((rP+1):nrow(B))) {
            
            if(all(Clones[[rP]]$mut[1:rP]==B[rC,1:rP])&&(sum(Clones[[rP]]$mut)==(sum(B[rC,])-1))) {

                mutName <- paste(mut_label[which(B[rC,-1]==1)],collapse="")
                Clones[[rC]] <- Clones[[rP]]$AddChild(mutName)
                Clones[[rC]]$mut <- B[rC,]
                
            }
            
        }
        
    }
    
    return(Broot)
    
}

# Draw B ts
draw.B.ts <- function(B, C) {
    
    # library(timescape)
    
    # Generate tree_edges
    adj_matrix <- as.adj.matrix.C_C(B)
    idx_woP <- which(colSums(adj_matrix)==0)
    
    tree_edg <- data.frame(source = character(), target = character())
    
    if(length(idx_woP)>1){
        for(i in 1:length(idx_woP)) {
            
            tree_edg <- rbind(tree_edg, data.frame(source = 'r', target = as.character(idx_woP[i])))
        }    
    }
    
    
    for(i in 1:nrow(adj_matrix)) {
        for(c in which(adj_matrix[i,]==1)){
            tree_edg <- rbind(tree_edg, data.frame(source = as.character(i), target = as.character(c)))
        }
    }
    
    cln_prev <- data.frame(timepoint = character(), clone_id = character(), clonal_prev = character())
    
    cln_prev <- expand.grid(timepoint = 1:length(C), clone_id = 1:nrow(adj_matrix))
    cln_prev <- cln_prev[order(cln_prev$timepoint),]
    cln_prev$clonal_prev <- 0
    # Generate clonal_prev
    for(c in 1:nrow(cln_prev)) {
        
        cln_prev$clonal_prev[c] <- sum(C[[cln_prev$timepoint[c]]]==cln_prev$clone_id[c]) / length(C[[cln_prev$timepoint[c]]])
        
    }
    
    
    timescape(clonal_prev = cln_prev, tree_edges = tree_edg, height=260, alpha=15)
    
}
