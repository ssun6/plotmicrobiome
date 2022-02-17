#' Parse metagenomics outputs from kraken, metaplan, humann
#' @param
#' @keywords format, wgs
#' @export
#' @examples
#'

format_pathway <- function(taxa_file = NULL,sep="\t",reads_cutoff=0,normalization=T,rarefy=F,rarefy_num=1000) {
  if(sep=="\t"){
    tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T,check.names=FALSE,quote="")
  }else{
    tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T,check.names=FALSE)
  }

  tab_all=na.omit(tab_all)
  tab_all=as.matrix(tab_all[,order(colnames(tab_all))])
  options(warn=-1)
  tab_all1=apply(tab_all,2,as.numeric)
  options(warn=0)
  rownames(tab_all1)=rownames(tab_all)
  if(!is.null(reads_cutoff)){
    tab_all1=tab_all1[,which(colSums(tab_all1)>reads_cutoff)]
  }
  tab_all1=tab_all1[!rowSums(tab_all1)==0,]
  if(normalization){
    if(rarefy){
      tab_all1 <- rarefy(tab_all1, rarefy_num)
    }else{
      tab_all1=t(t(tab_all1)/colSums(tab_all1))*mean(colSums(tab_all1))
    }
  }else{
    tab_all1=tab_all1
  }
  tab_all1=tab_all1[order(rowSums(tab_all1),decreasing = T),]
  return(tab_all1)
}
