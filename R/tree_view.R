#' Build taxonomic tree of the data and highlight the significant taxa
#' @keywords tree, highlight
#' @export
#' @examples
#'
#'
library(ggplot2)
library(ggpubr)
library(ggtree)
tree_view <- function(taxa_table = NULL, metadata=NULL,fdrs=NULL,test_metadata=NULL,test_metadata_continuous=F,single_parent_branch_removal=F,single_child_branch_removal=F,fdr_cutoff=0.1,node_size_breaks=c(0,0.01,0.05,0.5,5),taxa_removal="",palette_highlight=c("red","blue","orange","green"),node_size_limits=c(0,10),prevalence_cutoff=0.1,abundance_cutoff=0,domain="Bacteria") {

  if (prevalence_cutoff>0){
    taxa_table=taxa_table[which(apply(taxa_table,1,function(i){length(which(i!=0))})>=ncol(taxa_table)*prevalence_cutoff),]
  }

  if (abundance_cutoff>0){
    taxa_table=taxa_table[which(rowMeans(taxa_table)>abundance_cutoff),]
  }

  if (domain=="Bacteria"){
    taxa_table=taxa_table[grepl("__Bacteria",rownames(taxa_table)),]
  }else if (domain=="Archaea"){
    taxa_table=taxa_table[grepl("__Archaea",rownames(taxa_table)),]
  }else if (domain=="Eukaryota"){
    taxa_table=taxa_table[grepl("__Eukaryota",rownames(taxa_table)),]
  }

  metadata=metadata[which(!is.na(metadata[,test_metadata])),]
  map_s=metadata[intersect(colnames(taxa_table),rownames(metadata)),]
  taxa_table=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  if(!test_metadata_continuous){
    metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
    map_s[[test_metadata]]=droplevels(factor(map_s[[test_metadata]]))
    map_s[[test_metadata]]=gsub("/","_",map_s[[test_metadata]])
    map_s[[test_metadata]]=gsub(" ","_",map_s[[test_metadata]])
    fdrs[,4]=gsub("/","_",fdrs[,4])
    fdrs[,4]=gsub(" ","_",fdrs[,4])
  }

  fdrs_n=fdrs[!grepl("__--__--__",rownames(fdrs)),]
  names1=rownames(taxa_table)
  ln=sapply(strsplit(rownames(taxa_table),"--"),length)
  taxa_table=taxa_table[which(ln<8),]

  names7=names1[ln==7]

  l1=list()
  for (n in 7:3){
    l2=list()
    if(n==7){
      for (m in names1[which(ln==n-1)]){
        names2=names1[ln==n]
        l1[m]=paste0("(",paste(names2[grep(m,names2)],collapse=","),")",m)
      }
    }else{
      for (m in names1[which(ln==n-1)]){
        names2=names1[ln==n]
        l2[m]=paste0("(",paste(l1[names2[grep(m,names2)]],collapse=","),")",m)
      }
      l1=l2
    }
  }
  if(domain=="Bacteria"){
    l2=paste0("(",paste(l1,collapse=","),")Bacteria;")
  }else if (domain=="Archaea"){
    l2=paste0("(",paste(l1,collapse=","),")Archaea;")
  }else if (domain=="Eukaryota"){
    l2=paste0("(",paste(l1,collapse=","),")Eukaryota;")
  }else{
    stop("Please select Bacteria, Archaea or Eukaryota for domain.")
  }


  tree <- ape::read.tree(text = l2)

  fdrs_n=fdrs_n[which(!is.na(fdrs_n[,3])),]
  names_inter=intersect(rownames(fdrs_n),rownames(taxa_table))
  taxa_table1=taxa_table[match(names_inter,rownames(taxa_table)),]
  fdrs_n=fdrs_n[match(names_inter,rownames(fdrs_n)),]

  fdrs=as.numeric(fdrs_n[,3])
  tab_p=t(t(taxa_table1)/colSums(taxa_table))*100
  if(test_metadata_continuous){
    tab1_3=tab_p
    tab1_4=fdrs_n[,4]
    tab1_4[fdrs>fdr_cutoff]=NA
  }else{
    tab1_3=sapply(by(t(tab_p),map_s[[test_metadata]],colMeans),identity)
    tab1_4=fdrs_n[,4]
    tab1_4[fdrs>fdr_cutoff]=NA
  }

  tab1_5=apply(taxa_table1,1,function(i){length(which(i!=0))})/ncol(taxa_table1)
  ln1=sapply(strsplit(rownames(taxa_table1),"--"),length)
  a=data.frame(cbind(rownames(taxa_table1),rowMeans(tab1_3),tab1_5,tab1_4,fdrs,ln1))


  colnames(a)=c("taxa","abundance","prevalence","group","fdr","level")
  a$abundance=as.numeric(as.character(a$abundance))
  a$prevalence=as.numeric(as.character(a$prevalence))
  a$level=as.numeric(as.character(a$level))
  p=ggtree::ggtree(tree, layout='circular',edge.length=NULL)

  dd=a[!is.na(a$group),]
  dd=dd[!grepl("--__",dd$taxa) | dd$level>=6 ,]

  if (any(dd$level<6)){
    dd1=dd[dd$level<6,]

    node_num=vector()
    for (i in 1:nrow(dd1)){
      nodem=try(ggtree::MRCA(tree,tree$tip.label[grep(as.character(dd1[i,1]),tree$tip.label)]))
      if(class(nodem)=="try-error"){
        node_num[i]=NA
      }else{
        node_num[i]=nodem
      }
    }

    dd2=data.frame(cbind(as.character(dd1[,1]),node_num))
    colnames(dd2)=c("label","node_num")
    dd2$group=droplevels(factor(a$group[match(dd2$label,a$taxa)]))
    dd2$level=a$level[match(dd2$label,a$taxa)]
    dd2$ext=6-dd2$level

    if (single_parent_branch_removal==T){
      dd2=dd2[order(dd2$level,decreasing = T),]
      dd2=dd2[match(unique(dd2[,2]),dd2[,2]),]
    }

    if (single_child_branch_removal==T){
      dd2=dd2[order(dd2$level,decreasing = F),]
      dd2=dd2[match(unique(dd2[,2]),dd2[,2]),]
    }


    if(test_metadata_continuous){
      num_group=match(levels(factor(dd2$group)),c("positive","negative"))
    }else{
      num_group=match(levels(factor(dd2$group)),names(table(map_s[[test_metadata]])))
    }

    dd2$col=palette_highlight[num_group][factor(dd2$group)]
    dd2=na.omit(dd2)
    dd2$node_num=as.numeric(as.character(dd2$node_num))

    for (i in unique(dd2$level)){
      dd3=dd2[dd2$level==i,]
      ext1=(6-i)*0.6
      p = p + geom_hilight(data=dd3,mapping=aes(fill=group, node=node_num),alpha=0.2,extend=ext1) +
        scale_fill_manual(values=palette_highlight[num_group])
    }

    lett1=c(letters,paste0(rep(letters,each=10),rep(c(0:9),26)))
    dd4=dd2[dd2$level<6,]
    dd4=dd4[order(dd4$label),]

    dd4$letter1=rep("",nrow(dd4))
    if(any(dd4$level==2)){
      dd4$letter1[which(dd4$level==2)]=sapply(strsplit(as.character(dd4$label[dd4$level==2]),"p__"),"[[",2)
    }

    #dd4$letter1[which(dd4$level==3)]=sapply(strsplit(as.character(dd4$label[dd4$level==3]),"c__"),"[[",2)
    dd4$letter1[which(dd4$level>2)]=lett1[1:length(which(dd4$level>2))]

    dd4$label=as.character(dd4$label)
    if(nrow(dd4)==1){
      dd4$tax1=apply(dd4,1,function(i){strsplit(strsplit(dd4$label,"--")[[1]][as.numeric(dd4$level)],"__")[[1]][2]})
    }else{
      dd4$tax1=apply(dd4,1,function(i){strsplit(strsplit(i[1],"--")[[1]][as.numeric(i[4])],"__")[[1]][2]})
    }
    dd4=dd4[which(dd4$tax1!="__"),]

    if (!taxa_removal==""){
      if(any(dd4$tax1%in%taxa_removal)){
        dd4=dd4[-which(dd4$tax1%in%taxa_removal),]
      }
    }

    for (j in unique(dd4$level)){
      dd5=dd4[dd4$level==j,]
      ext1=(6-j)*0.6+0.1
      for (i in 1:nrow(dd5)){
        p=p + ggtree::geom_cladelabel(node = dd5$node_num[i], label = dd5$letter1[i],barsize =0,angle = 'auto',horizontal = F,offset.text = ext1,hjust = "center", fontsize = 4)
      }
    }

    dd6=dd4[dd4$level>2,]
    dd6$text1=paste(dd6$letter1,dd6$tax1)
    dd6$x=rep(1,nrow(dd6))
    dd6$y=c(nrow(dd6):1)
    if(test_metadata_continuous){
      num_group=match(levels(factor(dd6$group)),c("positive","negative"))
    }else{
      num_group=match(levels(factor(dd6$group)),names(table(map_s[[test_metadata]])))
    }


    if (length(table(dd6$col))==1){
      p1=ggplot(dd6,aes(x, y, label = text1))+geom_text(x=0,aes(colour = group),size = 2,hjust = 0)+ scale_color_manual(values=palette_highlight[num_group])+
        scale_y_continuous(limits = c(-20, nrow(dd6)+20))+
        theme(line = element_blank(),
              text = element_blank(),
              title = element_blank(),
              panel.background = element_blank(),
              legend.position = "none")
    }else{
      p1=ggplot(dd6,aes(x, y, label = text1))+geom_text(x=0,aes(colour = group),size = 3,hjust = 0)+ scale_color_manual(values=palette_highlight[num_group])+
        scale_y_continuous(limits = c(-20, nrow(dd6)+20))+
        theme(line = element_blank(),
              text = element_blank(),
              title = element_blank(),
              panel.background = element_blank(),
              legend.position = "none")
    }
  }

  lab1=c(paste0("<",node_size_breaks[2],"%"),paste0(node_size_breaks[2],"-",node_size_breaks[3],"%"),paste0(node_size_breaks[3],"-",node_size_breaks[4],"%")
         ,paste0(node_size_breaks[4],"-",node_size_breaks[5],"%"),paste0(">",node_size_breaks[5],"%"))

  if(test_metadata_continuous){
    num_group=match(levels(factor(dd$group)),c("positive","negative"))
  }else{
    num_group=match(levels(factor(dd$group)),names(table(map_s[[test_metadata]])))
  }
  p=p %<+% dd + geom_point(aes(size=abundance,color=group), alpha=1)+ scale_colour_manual(values=palette_highlight[num_group],na.translate = F)+scale_size_continuous(
    breaks = node_size_breaks,
    labels = lab1,
    limits = node_size_limits,
    guide = "legend",
    range= c(2,8)
  )

  p=p+theme(legend.position = "left")

  if (any(dd$level<6)){
    p2=ggarrange(p, p1, widths = c(3,1))
  }else{
    p2=p
  }
  return(p2)
}
