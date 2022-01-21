#' Build boxplots of alpha diversity
#' @keywords alpha-diversity, boxplots
#' @export
#' @examples
#'
#'
alpha_plot=function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,test_metadata_order="default",one_level=FALSE,method = "wilcoxon",xlab_direction=1,palette_group=c("red","blue","orange","green")){

  tab1=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  metadata=metadata[intersect(colnames(tab1),rownames(metadata)),]
  metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
  metadata[,test_metadata]=droplevels(metadata[,test_metadata])
  if(test_metadata_order[1]!="default"){
    metadata[,test_metadata]=factor(metadata[,test_metadata],levels=test_metadata_order)
  }

  if(one_level){
    shannon=vegan::diversity(tab1,index = "shannon", MARGIN = 2, base = exp(1))
    simp=vegan::diversity(tab1,index = "simpson", MARGIN = 2, base = exp(1))
    invsimp=vegan::diversity(tab1,index = "invsimpson", MARGIN = 2, base = exp(1))
    num_species=vegan::specnumber(tab1,  MARGIN = 2)

    alpha_mat=cbind(shannon,simp,invsimp,num_species)
    colnames(alpha_mat)=c("Shannon","Simpson","InverseSimpson","Number")

    par(mfrow=c(7,4),mar=c(12,5,5,5))
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

      if (is.na(wil_p1)){
        wil_p=wil_p1
      }else if(wil_p1<0.001){
        wil_p=formatC(wil_p1, format = "e", digits = 2)
      }else{
        wil_p=formatC(wil_p1, digits = 2)
      }

      boxplot(alpha_mat[,i]~metadata[,test_metadata],main=paste(colnames(alpha_mat)[i],"P =",wil_p),border=palette_group,col="white",cex.lab=1.5,cex.axis=1.2,xlab=test_metadata,ylab=colnames(alpha_mat)[i],las=xlab_direction)
      stripchart(alpha_mat[,i]~metadata[,test_metadata],vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = palette_group)
    }
  }else if (!one_level){
    par(mfrow=c(7,4),mar=c(8,5,5,5))
    tax_l=sapply(strsplit(rownames(tab1),"--"),function(i){length(i)})
    level1=c("kingdom","phylum","class","order","family","genus","species","ASV/strain")
    level_n=c(1:8)
    for (j in 2:max(tax_l)){
      tab1n=tab1[which(tax_l==j),]

      shannon=vegan::diversity(tab1n,index = "shannon", MARGIN = 2, base = exp(1))
      simp=vegan::diversity(tab1n,index = "simpson", MARGIN = 2, base = exp(1))
      invsimp=vegan::diversity(tab1n,index = "invsimpson", MARGIN = 2, base = exp(1))
      num_species=vegan::specnumber(tab1n,  MARGIN = 2)

      alpha_mat=cbind(shannon,simp,invsimp,num_species)
      colnames(alpha_mat)=c("Shannon","Simpson","InverseSimpson","Num_taxa")

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

        if (is.na(wil_p1)){
          wil_p=wil_p1
        }else if(wil_p1<0.001){
          wil_p=formatC(wil_p1, format = "e", digits = 2)
        }else{
          wil_p=formatC(wil_p1, digits = 2)
        }

        boxplot(alpha_mat[,i]~metadata[,test_metadata],main=paste(level1[j],colnames(alpha_mat)[i],"P =",wil_p),border=palette_group,col="white",cex.lab=1.5,cex.axis=1.2,xlab=test_metadata,ylab=colnames(alpha_mat)[i],las=xlab_direction)
        stripchart(alpha_mat[,i]~metadata[,test_metadata],vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = palette_group)
      }
    }
  }
}

