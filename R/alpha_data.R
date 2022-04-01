alpha_data=function(taxa_table = NULL,one_level=FALSE){
  
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
    
    alpha_mat1=cbind(shannon,simp,invsimp,num_species)
    colnames(alpha_mat1)=c("Shannon","Simpson","InverseSimpson","Number")
  }else if (!one_level){
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
      colnames(alpha_mat)=paste0(level1[j],c("_Shannon","_Simpson","_InverseSimpson","_Num_taxa"))
      
      if (j==2){
        alpha_mat1=alpha_mat
      }else{
        alpha_mat1=cbind(alpha_mat1,alpha_mat)
      }
    }
  }
  return(alpha_mat1)
}