#' Boxplots of individual taxa with stats
#' @keywords individual taxa, boxplot, stats
#' @export
#' @examples
#'
library(RColorBrewer)
taxa_barplot=function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,one_level=F,test_metadata_order="default",num_taxa=10,taxa_level="phylum",xlab_direction=1,legend_size=1.5,palette_group="default"){

  metadata=metadata[which(!is.na(metadata[,test_metadata])),]
  tab1=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  metadata=metadata[intersect(colnames(tab1),rownames(metadata)),]
  metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
  metadata[,test_metadata]=droplevels(metadata[,test_metadata])
  if(test_metadata_order[1]!="default"){
    metadata[,test_metadata]=factor(metadata[,test_metadata],levels=test_metadata_order)
  }

  if(one_level){
    tab1n=tab1
  }else{
    tax_l=sapply(strsplit(rownames(tab1),"--"),function(i){length(i)})
    level1=c("kingdom","phylum","class","order","family","genus","species","strain")
    level_n=c(1:8)
    tab1n=tab1[which(tax_l==level_n[match(taxa_level,level1)]),]
  }

  tab1n=t(t(tab1n)/colSums(tab1n))*100
  tab1n=tab1n[order(rowSums(tab1n),decreasing = T),]
  tab1n_c=data.frame(sapply(by(t(tab1n),metadata[,test_metadata],colMeans),identity))
  colnames(tab1n_c)=levels(metadata[,test_metadata])

  if(length(palette_group)==1){
    if(palette_group=="default"){
      qual_col_pals = RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual',]
      col11=c("red","blue","green","orange","turquoise1","deeppink","black","royalblue","darkgreen","hotpink","darkcyan","goldenrod1","brown","grey","purple")
      col_vector = c(col11,unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
      palette_group=col_vector[1:(num_taxa+1)]
    }else{
      warning("More colors are needed")
    }
  }else{
    palette_group=palette_group
  }

  if(!one_level){
    tab1n_c=tab1n_c[!grepl("--.*__$",rownames(tab1n_c),perl=T),]
    tab1n_c=tab1n_c[!grepl("--.__uncultured",rownames(tab1n_c)),]
    tab1n_c=tab1n_c[!grepl("--.__uncultured.bacterium",rownames(tab1n_c)),]
    tab1n_c=tab1n_c[!grepl("--.__gut.metagenome",rownames(tab1n_c)),]
    rownam1=sapply(strsplit(rownames(tab1n_c),"--"),"[[",level_n[match(taxa_level,level1)])
    tab1n_c1=aggregate(tab1n_c, by=list(rownam1), FUN=sum)
    rownames(tab1n_c1)=tab1n_c1[,1]
    tab1n_c1=tab1n_c1[,-1]
    tab1n_c=tab1n_c1[order(rowSums(tab1n_c1),decreasing = T),]
    #tab1n_c1=t(data.frame(sapply(by(tab1n_c,rownam1,rowSums),identity)))
    #rownames(tab1n_c)=sapply(strsplit(rownames(tab1n_c),"--"),"[[",level_n[match(taxa_level,level1)])
  }

  if(num_taxa>nrow(tab1n_c)){
    if(all(colSums(tab1n_c)==100)){
      num_taxa=nrow(tab1n_c)
      palette_group=palette_group[1:num_taxa]
      par(mfrow=c(1,2))
      barplot(as.matrix(tab1n_c),col=palette_group,ylab="Percentage (%)",las=xlab_direction)
      plot.new()
      legend("left",rev(rownames(tab1n_c)),col=rev(palette_group[1:num_taxa]),pch=15,bty="n",cex=legend_size)

    }else{
      num_taxa=nrow(tab1n_c)
      other=100-colSums(tab1n_c)
      tab1n_c2=rbind(tab1n_c,other)
      rownames(tab1n_c2)[num_taxa+1]="Other"
      par(mfrow=c(1,2))
      barplot(as.matrix(tab1n_c2),col=palette_group,ylab="Percentage (%)",las=xlab_direction)
      plot.new()
      legend("left",rev(rownames(tab1n_c2)),col=rev(palette_group[1:(num_taxa+1)]),pch=15,bty="n",cex=legend_size)

    }
  }else{
    tab1n_c1=tab1n_c[1:num_taxa,]
    other=100-colSums(tab1n_c1)
    tab1n_c2=rbind(tab1n_c1,other)
    rownames(tab1n_c2)[num_taxa+1]="Other"
    par(mfrow=c(1,2),mar=c(25,5,3,3))
    barplot(as.matrix(tab1n_c2),col=palette_group,ylab="Percentage (%)",las=xlab_direction,cex.lab=1.5,cex.axis = 1.5,cex.names=1.5)
    plot.new()
    legend("left",rev(rownames(tab1n_c2)),col=rev(palette_group[1:(num_taxa+1)]),pch=15,bty="n",cex=legend_size)
  }
}
