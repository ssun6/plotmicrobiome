#' Boxplots of individual taxa with stats
#' @keywords individual taxa, boxplot, stats
#' @export
#' @examples
#'
library(RColorBrewer)
taxa_barplot=function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,num_taxa=10,taxa_level="phylum",xlab_direction=1,legend_size=1,palette_group="default"){

  metadata=metadata[which(!is.na(metadata[,test_metadata])),]
  tab1=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  metadata=metadata[intersect(colnames(tab1),rownames(metadata)),]
  metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
  metadata[,test_metadata]=droplevels(metadata[,test_metadata])

  tax_l=sapply(strsplit(rownames(tab1),"--"),function(i){length(i)})
  level1=c("kingdom","phylum","class","order","family","genus","species","strain")
  level_n=c(1:8)
  tab1n=tab1[which(tax_l==level_n[match(taxa_level,level1)]),]

  tab1n=t(t(tab1n)/colSums(tab1n))*100
  tab1n=tab1n[order(rowSums(tab1n),decreasing = T),]
  tab1n_c=data.frame(sapply(by(t(tab1n),metadata[,test_metadata],colMeans),identity))

  if(length(palette_group)==1){
    if(palette_group=="default"){
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      palette_group=col_vector[1:(num_taxa+1)]
    }else{
      warning("More colors are needed")
    }
  }else{
    palette_group=palette_group
  }

  rownames(tab1n_c)=sapply(strsplit(rownames(tab1n_c),"--"),"[[",level_n[match(taxa_level,level1)])
  if(num_taxa>nrow(tab1n_c)){
    num_taxa=nrow(tab1n_c)
    palette_group=palette_group[1:num_taxa]
    par(mfrow=c(1,2))
    barplot(as.matrix(tab1n_c),col=palette_group,ylab="Percentage (%)",las=xlab_direction)
    plot.new()
    legend("left",rev(rownames(tab1n_c)),col=rev(palette_group[1:num_taxa]),pch=15,bty="n",cex=legend_size)

  }else{
    tab1n_c1=tab1n_c[1:num_taxa,]
    other=colSums(tab1n_c)-colSums(tab1n_c1)
    tab1n_c2=rbind(tab1n_c1,other)
    rownames(tab1n_c2)[num_taxa+1]="Other"
    par(mfrow=c(1,2))
    barplot(as.matrix(tab1n_c2),col=palette_group,ylab="Percentage (%)",las=xlab_direction)
    plot.new()
    legend("left",rev(rownames(tab1n_c2)),col=rev(palette_group[1:(num_taxa+1)]),pch=15,bty="n",cex=legend_size)

  }

}
