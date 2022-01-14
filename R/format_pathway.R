#' Parse metagenomics outputs from kraken, metaplan, humann
#' @param
#' @keywords format, wgs
#' @export
#' @examples
#'

format_pathway <- function(taxa_file = NULL,sep="\t",reads_cutoff=0,rarefy=F,rarefy_num=1000) {
  if(sep=="\t"){
    tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T,check.names=FALSE,quote="")
  }else{
    tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T,check.names=FALSE)
  }

  tab_all=as.matrix(tab_all[,order(colnames(tab_all))])
  if(!is.null(reads_cutoff)){
    tab_all=tab_all[,which(colSums(tab_all)>reads_cutoff)]
  }
  tab_all=tab_all[!rowSums(tab_all)==0,]
  if(rarefy){
    tab_all <- rarefy(tab_all, rarefy_num)
  }else{
    tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))
  }
  tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]
  return(tab_all)
}
