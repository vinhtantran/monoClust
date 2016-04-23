perm.test <- function(object, data, auto.pick = FALSE, sig.val = 0.05, ...){
    
    if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")
    
    frame <- object$frame
    node <- as.numeric(row.names(frame))
    
    # It would be better to create a jump table for reference of tree walking
    jump.table <- cbind(node, frame[,c("bipartvar", "cut")])
    jump.table$right <- jump.table$left <- NA
    jump.table$p.value <- NA
    
    for (i in 2:nrow(jump.table)) {
        parent.node <- jump.table$node[i] %/% 2
        is.left <- jump.table$node[i] %% 2 == 0
        
        if (is.left) {
            jump.table$left[which(jump.table$node == parent.node)] <- i
        } else
            jump.table$right[which(jump.table$node == parent.node)] <- i
    }
    
    # Now tracing the tree to find the cluster
    # Tracing the tree by Node-Left-Right algorithm
    assign(".Jump_Table", jump.table, envir = .GlobalEnv)
    assign(".Data", data, envir = .GlobalEnv)
    recursive.walk(1, 1:nrow(data), auto.pick, sig.val)
    
    jump.table <- .Jump_Table
    rm(list = c(".Jump_Table", ".Data"), envir = globalenv())
    frame$p.value <- jump.table$p.value
    if (auto.pick) {
        a <- frame[which(frame$p.value > sig.val),]
        last.split <- min(a$split.order)
        object$numofclusters <- last.split
    }
    object$frame <- frame
    return(object)
}

recursive.walk <- function(current, members, auto.pick, sig.val) {
    # Check stopping at leaf node
    if (is.na(.Jump_Table$bipartvar[current])) return(0)
    
    # Node
    split.var <- .Jump_Table$bipartvar[current]
    split.value <- .Jump_Table$cut[current]
    
    data.temp <- .Data
    data.temp[,split.var] <- NULL
    distmat.reduced <- as.matrix(daisy(data.frame(data.temp)))
    
    members.L <- members[which(.Data[members, split.var] < split.value)]
    members.R <- setdiff(members, members.L)
    
    dist.mat.twogroup <- distmat.reduced[c(members.L, members.R),c(members.L, members.R)]
    fmem2 <- factor(c(rep(1, length(members.L)), rep(2, length(members.R))))
    result <- adonis(dist.mat.twogroup ~ fmem2)
    pvalue.adj <- (as.numeric(row.names(.Jump_Table[current,])) %/% 2 + 1) * result$aov.tab[1,6]
    .Jump_Table$p.value[current] <<- ifelse(pvalue.adj > 1, 1, pvalue.adj)
    
    if (auto.pick && (pvalue.adj > sig.val)) return(0)    
    # Left
    recursive.walk(.Jump_Table$left[current], members.L, auto.pick, sig.val)
    
    # Right
    recursive.walk(.Jump_Table$right[current], members.R, auto.pick, sig.val)
}