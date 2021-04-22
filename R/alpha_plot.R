#' Build boxplots of alpha diversity
#' @keywords alpha-diversity, boxplots
#' @export
#' @examples
#'
#'
alpha_plot=function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,taxa_level="species",method = "wilcoxon",xlab_direction=1,palette_group=c("red","blue","orange","green")){

  tab1=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  metadata=metadata[intersect(colnames(tab1),rownames(metadata)),]
  metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
  metadata[,test_metadata]=droplevels(metadata[,test_metadata])

  tax_l=sapply(strsplit(rownames(tab1),"--"),function(i){length(i)})
  level1=c("kingdom","phylum","class","order","family","genus","species","strain")
  level_n=c(1:8)
  tab1n=tab1[which(tax_l==level_n[match(taxa_level,level1)]),]

  shannon=vegan::diversity(tab1n,index = "shannon", MARGIN = 2, base = exp(1))
  simp=vegan::diversity(tab1n,index = "simpson", MARGIN = 2, base = exp(1))
  invsimp=vegan::diversity(tab1n,index = "invsimpson", MARGIN = 2, base = exp(1))
  num_species=vegan::specnumber(tab1n,  MARGIN = 2)

  alpha_mat=cbind(shannon,simp,invsimp,num_species)
  colnames(alpha_mat)=c("shannon","simp","invsimp","num_species")

  par(mfrow=c(2,2),mar=c(5,5,5,5))
  for (i in 1:4){

    if (method == "wilcoxon"){
      wil_p1=try(wilcox.test(as.numeric(alpha_mat[,i])~metadata[,test_metadata])$p.value)
    }else if (method == "t.test"){
      wil_p1=try(t.test(as.numeric(alpha_mat[,i])~metadata[,test_metadata])$p.value)
    }else if (method == "anova"){
      wil_p1=try(summary(aov(as.numeric(alpha_mat[,i])~metadata[,test_metadata]))[[1]][1,5])
    }else{
      wil_p1=kruskal.test(as.numeric(alpha_mat[,i])~metadata[,test_metadata])$p.value
    }

    if(wil_p1<0.001){
      wil_p=formatC(wil_p1, format = "e", digits = 2)
    }else{
      wil_p=formatC(wil_p1, digits = 2)
    }

    boxplot(alpha_mat[,i]~metadata[,test_metadata],main=paste(colnames(alpha_mat)[i],"P =",wil_p),border=palette_group,xlab=test_metadata,ylab=colnames(alpha_mat)[i],las=xlab_direction)
    stripchart(alpha_mat[,i]~metadata[,test_metadata],vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = palette_group)
  }
}

