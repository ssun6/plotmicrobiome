#' Parse metagenomics outputs from kraken, metaplan, humann
#' @param
#' @keywords format, wgs
#' @export
#' @examples
#'

format_wgs <- function(taxa_file = NULL,sep="\t",method="metaphlan",normalization=T,reads_cutoff=0,rarefy=F,rarefy_num=1000) {
  tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T,check.names=FALSE,quote="",comment.char = "")
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


  tab_all=tab_all[order(rownames(tab_all),decreasing = F),]
  rownames(tab_all)=taxa_edit(rownames(tab_all))
  tab_all=tab_all[which(!grepl("__--__--__",rownames(tab_all))),]

  if (method=="kraken"){
    #Bacteria and Archaea are domains, Fungi is a kingdom, Virus has many kingdoms
    #The taxonomic structuring are different, the Eukaryota level is removed and the kingdom levels of virus are removed to be consistent across branches
    tab_all1=tab_all[!grepl("d__Eukaryota",rownames(tab_all)) & !grepl("d__Viruses",rownames(tab_all)),]
    tab_all2=tab_all[grep("d__Eukaryota",rownames(tab_all)),][-1,]
    tab_all3=tab_all[grep("d__Viruses",rownames(tab_all)),]

    taxa_l=c("d__","p__","c__","o__","f__","g__","s__","t__")
    miss_row=vector()
    j=1
    name_new1=vector()
    for (i in 1:nrow(tab_all1)){
      n1=strsplit(rownames(tab_all1)[i],"--")[[1]]
      n2=length(n1)
      if(!grepl(taxa_l[n2],n1[n2])){
        print(i)
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
      rownames(tab_all1)[miss_row]=name_new1
    }

    rownames(tab_all2)=gsub("d__Eukaryota--","",rownames(tab_all2))
    taxa_l=c("k__","p__","c__","o__","f__","g__","s__","t__")
    miss_row=vector()
    j=1
    name_new1=vector()
    for (i in 1:nrow(tab_all2)){
      n1=strsplit(rownames(tab_all2)[i],"--")[[1]]
      n2=length(n1)
      if(!grepl(taxa_l[n2],n1[n2])){
        print(i)
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
      rownames(tab_all2)[miss_row]=name_new1
    }

    taxa_l=c("d__","p__","c__","o__","f__","g__","s__","t__")
    miss_row=vector()
    j=1
    name_new1=vector()
    for (i in 1:nrow(tab_all3)){
      n1=strsplit(rownames(tab_all3)[i],"--")[[1]]
      n2=length(n1)
      if(!grepl(taxa_l[n2],n1[n2]) | grepl("k__",rownames(tab_all3)[i])){
        print(i)
        miss_row[j]=i
        match1=match(paste0(sapply(strsplit(n1,"__"),"[[",1),"__"),taxa_l)
        name_new=vector()
        if(any(is.na(match1))){
          name_new[na.omit(match1)]=n1[-which(is.na(match1))]
          match1=na.omit(match1)
        }else{
          name_new[match1]=n1
        }
        name_new[which(is.na(name_new))]=taxa_l[which(is.na(name_new))]
        name_new1[j]=paste(name_new,collapse = "--")
        j=j+1
      }
    }
    if(length(miss_row)!=0){
      rownames(tab_all3)[miss_row]=name_new1
    }
    if(any(name_new1=="d__Viruses")){
      tab_all3=tab_all3[-miss_row[which(name_new1=="d__Viruses")],]
    }
    tab_all=rbind(tab_all1,tab_all2,tab_all3)

  }else if (method=="metaphlan"){
    tab_all=tab_all
  }
  return(tab_all)
}
