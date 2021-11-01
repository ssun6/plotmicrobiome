#' Parse metagenomics outputs from kraken, metaplan, humann
#' @param
#' @keywords format, wgs
#' @export
#' @examples
#'

format_wgs <- function(taxa_file = NULL,sep="\t",method="metaphlan") {
  tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T,quote="",check.names=FALSE)
  tab_all=tab_all[,order(colnames(tab_all))]
  tab_all=tab_all[grep("Bacteria",rownames(tab_all)),]
  tab_all=tab_all[!rowSums(tab_all)==0,]
  tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))
  tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]
  rownames(tab_all)=taxa_edit(rownames(tab_all))
  tab_all=tab_all[which(!grepl("__--__--__",rownames(tab_all))),]

  if (method=="kraken"){
    taxa_l=c("d__","p__","c__","o__","f__","g__","s__","t__")
    miss_row=vector()
    j=1
    for (i in 1:nrow(tab_all)){
      n1=strsplit(rownames(tab_all)[i],"--")[[1]]
      n2=length(n1)
      if(!grepl(taxa_l[n2],n1[n2])){
        miss_row[j]=i
        j=j+1
      }
    }
    if(length(miss_row)!=0){
      tab_all=tab_all[-miss_row,]
    }
  }else if (method=="metaphlan"){
    tab_all=tab_all
  }
  return(tab_all)
}
