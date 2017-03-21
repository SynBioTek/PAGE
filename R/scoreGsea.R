#' Gene Expression Signature Similarity based on GSEA
#' 
#' Gene Expression Signature Similarity based on  Gene Set Enrichment Analysis
#' 
#' @param dat gene rank matrix or data.frame
#' @param top top genes selected
#' @param ncore the number of cores used
#' @param type avg or max
#' @export scoreGsea
#'  
#' @examples 
#' data(Psoriasis_Etanercept_LogFC)
#' Psoriasis_Etanercept_rank <- apply(-Psoriasis_Etanercept_LogFC, 2, rank)
#' 
#' sim <- scoreGsea(Psoriasis_Etanercept_rank, type="avg")
#' sim <- scoreGseam(Psoriasis_Etanercept_rank, type="avg", ncore=2)
#' 
scoreGsea <- function(dat, top=250, type=c("avg", "max") ) {
    
    type <- match.arg(type, c("avg", "max"))
    
    ES <- matrix(-2, ncol=ncol(dat), nrow=ncol(dat) )
    for (i in 1:ncol(dat) ) {
        for (j in 1:ncol(dat) ) {
            ES[i,j] <- quickenrichmentscore(which(dat[,j]<=top), which(dat[,j]>=nrow(dat)-top+1), dat[,i])
        }
    }
    if (type=="avg") {
        sim <- (ES + t(ES))/2
    } else if (type=="max") {
        sim <- pmax(ES, t(ES))
    }
    rownames(sim) <- colnames(sim) <- colnames(dat)
    
    return (sim)
    
}

#' @inheritParams scoreGsea
#' @import methods
#' @import parallel
#' @import foreach
#' @import doParallel
#' @export scoreGseam
#' @rdname scoreGsea
#' 
scoreGseam <- function(dat, top=250, type=c("avg", "max"), ncore=2) {
    
    type <- match.arg(type, c("avg", "max"))
    
    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)

    ES <- foreach (j = 1:ncol(dat), .combine = cbind) %:% 
        foreach(i = 1:ncol(dat), .combine = c) %dopar% {
        quickenrichmentscore(which(dat[,j]<=top), which(dat[,j]>=nrow(dat)-top+1), dat[,i])
    }
    
    parallel::stopCluster(cl)
    
    if (type=="avg") {
        sim <- (ES + t(ES))/2
    } else if (type=="max") {
        sim <- pmax(ES, t(ES))
    }
    rownames(sim) <- colnames(sim) <- colnames(dat)
    
    return (sim)
    
    
}



