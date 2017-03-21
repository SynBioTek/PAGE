#' quickenrichmentscore
#' 
#' quickenrichmentscore
#' 
#' @param S rank of up-regualted genes
#' @param S1 rank of down-regulated genes
#' @param List a rank vector
#' @return similarity based on averaged GSEA Enrichment score
#' @export
#' @author Yang Cao
#' @source http://bioconductor.org/packages/GeneExpressionSignature/
#' @examples 
#' data(Psoriasis_Etanercept_LogFC)
#' dat <- apply(-Psoriasis_Etanercept_LogFC, 2, rank)
#' SignatureLength <- 250
#' 
#' # caculate similarity only for some columns
#' sim_mat <- matrix(-2, ncol=ncol(dat), nrow=ncol(dat))
#' for (i in 11:15) {
#'     for (j in 1:10) {
#'         sim_mat[j,i] <- (quickenrichmentscore(which(dat[,j]<=SignatureLength), 
#'         which(dat[,j]>=nrow(dat)-SignatureLength+1), dat[,i]) + 
#'         quickenrichmentscore(which(dat[,i]<=SignatureLength), 
#'         which(dat[,i]>=nrow(dat)-SignatureLength+1), dat[,j]) )/2
#'     }
#' 
#' }
#' 
quickenrichmentscore <- function(S,S1,List){
	Or_list=List;
	#List=sort(List);
	#List=as.matrix(List)
	Rank=order(Or_list)
	Rank=as.matrix(Rank)
	N=length(Rank)
	Nh=length(S)
	tmp=matrix(0,N)
	for (i in 1:Nh){
		tmp[Or_list[S[i]]]=1
	}
	hitCases=cumsum(tmp)
	missCases=cumsum(1-tmp)
	NR=length(S)
	Phit=hitCases/NR
	Pmiss=missCases/(N-Nh)
	abshm=abs(Phit-Pmiss)
	abshm=as.matrix(abshm)
	#m=apply(abshm,2,max);
	t=apply(abshm,2,which.max)
	ESUP = Phit[t]-Pmiss[t];
	RS = Phit-Pmiss;

	Or_list2=Or_list
	Rank=order(Or_list2)
	Rank=as.matrix(Rank) 
	N=length(Rank)
	Nh=length(S1)
	tmp=matrix(0,N)
	for (i in 1:Nh){
		tmp[Or_list[S1[i]]]=1
	}
	hitCases=cumsum(tmp)
	missCases=cumsum(1-tmp) 
	NR=length(S1)
	Phit=hitCases/NR
	Pmiss=missCases/(N-Nh)
	abshm=abs(Phit-Pmiss)
	abshm=as.matrix(abshm) 
	#m=apply(abshm,2,max);
	t=apply(abshm,2,which.max)
	ESDOWN = Phit[t]-Pmiss[t];
	RS = Phit-Pmiss;
	ES = (ESUP - ESDOWN)/2;
	return (ES)
}
