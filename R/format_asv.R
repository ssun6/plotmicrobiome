#' Parse various 16S outputs
#' biom file, multiple biom files at different taxonomic levels, one or separate tsv or csv files
#' @param
#' @keywords format, 16s
#' @export
#' @examples
#'
library(rhdf5)
library(rbiom)
taxa_edit=function(list1){
  list1=gsub("\\.D","--D",list1)
  list1=gsub("\\.__","--__",list1)
  list1=gsub("\\|",";",list1)
  list1=gsub("; ","--",list1)
  list1=gsub(";","--",list1)
  list1=gsub(" ","_",list1)
  list1=gsub("\\[","",list1)
  list1=gsub("\\]","",list1)
  list1=gsub(":","",list1)
  list1=gsub("D_0","d",list1)
  list1=gsub("D_1","p",list1)
  list1=gsub("D_2","c",list1)
  list1=gsub("D_3","o",list1)
  list1=gsub("D_4","f",list1)
  list1=gsub("D_5","g",list1)
  list1=gsub("D_6","s",list1)
  list1=gsub("\\(","_",list1)
  list1=gsub("\\)","_",list1)
  list1=gsub("\\'","",list1)
  return(list1)
}

format_asv <- function(taxa_file = NULL,sep="\t",onefile=F,biom=F,ASV=F) {
  if (onefile){
    if (biom & ASV){
      biom= rbiom::read.biom(taxa_file,tree=FALSE)
      tab=as.matrix(biom$counts)
      tax_l=biom$taxonomy
      tax_l[tax_l==""]="__"
      tax1=apply(tax_l[,1:2],1,function(i){paste(i,collapse=";")})
      tab_all=t(data.frame(sapply(by(tab,tax1,colSums),identity)))
      tab_all=tab_all[!rowSums(tab_all)==0,]
      tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))

      for (n in 3:7){
        tax1=apply(tax_l[,1:n],1,function(i){paste(i,collapse=";")})
        tab_n=t(data.frame(sapply(by(tab,tax1,colSums),identity)))
        tab_n=tab_n[!rowSums(tab_n)==0,]
        tab_n=t(t(tab_n)/colSums(tab_n))*mean(colSums(tab_n))
        tab_all=rbind(tab_all,tab_n)
      }
      tax_asv=paste(apply(tax_l,1,function(i){paste(i,collapse=";")}),rownames(biom$taxonomy),sep=";")
      tab_asv=tab[!rowSums(tab)==0,]
      tab_asv=t(t(tab_asv)/colSums(tab_asv))*mean(colSums(tab_asv))
      rownames(tab_asv)=tax_asv
      tab_all=rbind(tab_all,tab_asv)

      tab_all=tab_all[,order(colnames(tab_all))]
      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]

    }else if (!biom & ASV){
      message ("If the taxa abundance table was converted from biom file, please remove # from header")
      tab1=read.table(file=taxa_file,sep=sep,row.names=1,header = T)
      tax1=as.character(tab1[,ncol(tab1)])
      tab=tab1[,-ncol(tab1)]
      tax_l=matrix(nrow=nrow(tab),ncol=7)
      for (i in 1:nrow(tab)){
        n=length(strsplit(tax1[i],"; ")[[1]])
        tax_l[i,1:n]=strsplit(tax1[i],"; ")[[1]][1:n]
      }
      tax_l[is.na(tax_l)]="__"
      tax1=apply(tax_l[,1:2],1,function(i){paste(i,collapse=";")})
      tab_all=t(sapply(by(tab,tax1,colSums),identity))
      tab_all=tab_all[!rowSums(tab_all)==0,]
      tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))

      for (n in 3:7){
        tax1=apply(tax_l[,1:n],1,function(i){paste(i,collapse=";")})
        tab_n=t(data.frame(sapply(by(tab,tax1,colSums),identity)))
        tab_n=tab_n[!rowSums(tab_n)==0,]
        tab_n=t(t(tab_n)/colSums(tab_n))*mean(colSums(tab_n))
        tab_all=rbind(tab_all,tab_n)
      }
      tax_asv_name=paste(apply(tax_l,1,function(i){paste(i,collapse=";")}),rownames(tab),sep=";")
      tab_asv=tab[!rowSums(tab)==0,]
      tab_asv=t(t(tab_asv)/colSums(tab_asv))*mean(colSums(tab_asv))
      rownames(tab_asv)=tax_asv_name
      tab_all=rbind(tab_all,tab_asv)
      tab_all=tab_all[,order(colnames(tab_all))]
      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]

    }else if (biom & !ASV){
      biom= rbiom::read.biom(taxa_file,tree=FALSE)
      tab_all=as.matrix(biom$counts)
      tab_all=tab_all[,order(colnames(tab_all))]
      tab_all=tab_all[!rowSums(tab_all)==0,]
      tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))
      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]

    }else{
      message ("If the taxa abundance table was converted from biom file, please remove # from header")
      tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T)
      tab_all=tab_all[,order(colnames(tab_all))]
      tab_all=tab_all[!rowSums(tab_all)==0,]
      tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))
      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]
    }
  }else{
    if (biom){
      file_list=list.files(taxa_file,pattern = ".biom")
      for (f1 in file_list){
        biom= rbiom::read.biom(paste0(taxa_file,"/",f1),tree=FALSE)
        tab=as.matrix(biom$counts)
        tab1=tab[,order(colnames(tab))]
        tab1=tab1[!rowSums(tab1)==0,]

        tab1=t(t(tab1)/colSums(tab1))*mean(colSums(tab1))
        tab1=tab1[order(rowSums(tab1),decreasing = T),]
        if (f1==file_list[1]){
          tab_all=tab1
        }else{
          tab_all=rbind(tab_all,tab1)
        }
        tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]
      }
    }else{
      message ("If the taxa abundance table was converted from biom file, please remove # from header")
      file_list=list.files(taxa_file)
      for (f1 in file_list){
        tab=read.table(file=paste0(taxa_file,"/",f1),sep=sep,row.names=1,header = T)
        tab1=tab[,order(colnames(tab))]
        tab1=tab1[!rowSums(tab1)==0,]

        tab1=t(t(tab1)/colSums(tab1))*mean(colSums(tab1))
        tab1=tab1[order(rowSums(tab1),decreasing = T),]
        if (f1==file_list[1]){
          tab_all=tab1
        }else{
          tab_all=rbind(tab_all,tab1)
        }
      }
      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]
    }
  }
  rownames(tab_all)=taxa_edit(rownames(tab_all))
  return(tab_all)
}

