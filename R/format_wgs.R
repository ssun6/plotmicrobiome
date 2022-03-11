#' Parse metagenomics outputs from kraken, metaplan, humann
#' @param
#' @keywords format, wgs
#' @export
#' @examples
#'

format_wgs <- function(taxa_file = NULL,sep="\t",method="metaphlan",normalization=T,reads_cutoff=0,rarefy=F,rarefy_num=1000) {
  tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T,quote="",check.names=FALSE)
  tab_all=as.matrix(tab_all[,order(colnames(tab_all))])
  if(!is.null(reads_cutoff)){
    tab_all=tab_all[,which(colSums(tab_all)>reads_cutoff)]
  }
  tab_all=tab_all[!rowSums(tab_all)==0,]
  if (normalization){
    if(rarefy){
      tab_all <- rarefy(tab_all, rarefy_num)
    }else{
      tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))
    }
  }else{
    tab_all=tab_all
  }


  tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]
  rownames(tab_all)=taxa_edit(rownames(tab_all))
  tab_all=tab_all[which(!grepl("__--__--__",rownames(tab_all))),]

  if (method=="kraken"){
    rownames(tab_all)=gsub("--k__Fungi","",rownames(tab_all))
    taxa_l=c("d__","p__","c__","o__","f__","g__","s__","t__")
    miss_row=vector()
    j=1
    name_new1=vector()
    for (i in 1:nrow(tab_all)){
      n1=strsplit(rownames(tab_all)[i],"--")[[1]]
      n2=length(n1)
      if(!grepl(taxa_l[n2],n1[n2])){
        miss_row[j]=i
        match1=match(paste0(sapply(strsplit(n1,"__"),"[[",1),"__"),taxa_l)
        name_new=vector()
        name_new[match1]=n1
        name_new[which(is.na(name_new))]=taxa_l[which(is.na(name_new))]
        name_new1[j]=paste(name_new,collapse = "--")
        j=j+1
      }
    }
    if(length(miss_row)!=0){
      rownames(tab_all)[miss_row]=name_new1
    }
  }else if (method=="metaphlan"){
    tab_all=tab_all
  }
  return(tab_all)
}
