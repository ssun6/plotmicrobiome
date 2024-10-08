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

format_asv <- function(taxa_file = NULL,sep="\t",biom=T,ASV=T,normalization=T,reads_cutoff=0,rarefy=F,rarefy_num=1000) {
    if (biom & ASV){
      biom= rbiom::read.biom(taxa_file,tree=FALSE)
      tab1=as.matrix(biom$counts)
      if(!is.null(reads_cutoff)){
        tab1=tab1[,which(colSums(tab1)>reads_cutoff)]
      }

      if(normalization){
        if(rarefy){
          tab <- rbiom::rarefy(tab1, rarefy_num)
        }else{
          tab=t(t(tab1)/colSums(tab1))*mean(colSums(tab1))
        }
      }else{
        tab=tab1
      }


      #match the taxonomy of rarefied table
      tax_l=biom$taxonomy
      tax_l[tax_l==""]="__"
      tax_l=tax_l[match(rownames(tab),rownames(tab1)),]
      flag1=any(grepl("p__",tax_l)) & any(grepl("c__",tax_l)) & any(grepl("o__",tax_l))
      if(!flag1){
       stop("No taxonomic structure in this table. Please select 'table without taxonomic structure' as data type or check formats at https://github.com/ssun6/plotmicrobiome")
      }


      tax1=apply(tax_l[,1:2],1,function(i){paste(i,collapse=";")})
      tab_all=t(sapply(by(tab,tax1,colSums),identity))
      tab_all=tab_all[!rowSums(tab_all)==0,]

      for (n in 3:7){
        tax1=apply(tax_l[,1:n],1,function(i){paste(i,collapse=";")})
        tab_n=t(sapply(by(tab,tax1,colSums),identity))
        tab_n=tab_n[!rowSums(tab_n)==0,]
        tab_all=rbind(tab_all,tab_n)
      }
      tax_asv=paste(apply(tax_l,1,function(i){paste(i,collapse=";")}),rownames(tab),sep=";") #formatted taxa
      tab_asv=tab #count table
      rownames(tab_asv)=tax_asv
      tab_all=rbind(tab_all,tab_asv)

      tab_all=tab_all[,order(colnames(tab_all))]
      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]

    }else if (!biom & ASV){
      tab1=read.table(file=taxa_file,sep=sep,row.names=1,header = T,check.names=FALSE,quote="",comment.char = "")
      tax1=as.character(tab1[,ncol(tab1)])
      tab1=as.matrix(tab1[,-ncol(tab1)])

      if(!is.null(reads_cutoff)){
        tab1=tab1[,which(colSums(tab1)>reads_cutoff)]
      }

      if(normalization){
        if(rarefy){
          tab <- rbiom::rarefy(tab1, rarefy_num)
        }else{
          tab=t(t(tab1)/colSums(tab1))*mean(colSums(tab1))
        }
      }else{
        tab=tab1
      }


      tax_l=matrix(nrow=nrow(tab1),ncol=7)
      for (i in 1:nrow(tab)){
        n=length(strsplit(tax1[i],"; ")[[1]])
        if(n>7){
          n=7
        }
        tax_l[i,1:n]=strsplit(tax1[i],"; ")[[1]][1:n]
      }
      tax_l[is.na(tax_l)]="__"
      tax_l=tax_l[match(rownames(tab),rownames(tab1)),]

      #flag for taxa structure
      flag1=any(grepl("p__",tax1)) & any(grepl("c__",tax1)) & any(grepl("o__",tax1))
      if(!flag1){
        stop("No taxonomic structure in this table. Please select 'table without taxonomic structure' as data type or check formats at https://github.com/ssun6/plotmicrobiome")
      }


      tax1=apply(tax_l[,1:2],1,function(i){paste(i,collapse=";")})

      tab_all=t(sapply(by(tab,tax1,colSums),identity))
      tab_all=tab_all[!rowSums(tab_all)==0,]

      for (n in 3:7){
        tax1=apply(tax_l[,1:n],1,function(i){paste(i,collapse=";")})
        tab_n=t(sapply(by(tab,tax1,colSums),identity))
        tab_n=tab_n[!rowSums(tab_n)==0,]
        tab_all=rbind(tab_all,tab_n)
      }
      tax_asv_name=paste(apply(tax_l,1,function(i){paste(i,collapse=";")}),rownames(tab),sep=";")
      tab_asv=tab
      rownames(tab_asv)=tax_asv_name
      tab_all=rbind(tab_all,tab_asv)
      tab_all=tab_all[,order(colnames(tab_all))]
      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]

    }else if (biom & !ASV){
      biom= rbiom::read.biom(taxa_file,tree=FALSE)
      tab_all=as.matrix(biom$counts)
      if(!is.null(reads_cutoff)){
        tab_all=tab_all[,which(colSums(tab_all)>reads_cutoff)]
      }
      tab_all=tab_all[,order(colnames(tab_all))]
      tab_all=tab_all[!rowSums(tab_all)==0,]
      if(normalization){
        if(rarefy){
          tab_all <- rbiom::rarefy(tab_all, rarefy_num)
        }else{
          tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))
        }
      }

      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]
      #flag for taxa structure
      flag1=any(grepl("p__",rownames(tab_all))) & any(grepl("c__",rownames(tab_all))) & any(grepl("o__",rownames(tab_all)))
      if(!flag1){
        stop("No taxonomic structure in this table. Please select 'table without taxonomic structure' as data type or check formats at https://github.com/ssun6/plotmicrobiome")
      }

    }else{
      tab_all=read.table(file=taxa_file,sep=sep,row.names=1,header = T,check.names=FALSE,quote="",comment.char = "")
      tab_all=as.matrix(tab_all[,order(colnames(tab_all))])
      if(!is.null(reads_cutoff)){
        tab_all=tab_all[,which(colSums(tab_all)>reads_cutoff)]
      }
      tab_all=tab_all[!rowSums(tab_all)==0,]
      if(normalization){
        if(rarefy){
          tab_all <- rbiom::rarefy(tab_all, rarefy_num)
        }else{
          tab_all=t(t(tab_all)/colSums(tab_all))*mean(colSums(tab_all))
        }
      }

      tab_all=tab_all[order(rowSums(tab_all),decreasing = T),]
      #flag for taxa structure
      flag1=any(grepl("p__",rownames(tab_all))) & any(grepl("c__",rownames(tab_all))) & any(grepl("o__",rownames(tab_all)))
      if(!flag1){
        stop("No taxonomic structure in this table. Please select 'table without taxonomic structure' as data type or check formats at https://github.com/ssun6/plotmicrobiome")
      }
  }
  rownames(tab_all)=taxa_edit(rownames(tab_all))
  tab_all=tab_all[which(rowSums(tab_all)!=0),]
  return(tab_all)
}

