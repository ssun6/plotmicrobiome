taxa_bar_table=function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,one_level=F,test_metadata_order="default",num_taxa=10,taxa_level="phylum",xlab_direction=1,legend_size=1.5,palette_group="default"){

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

  if(!one_level){
    tab1n_c=tab1n_c[!grepl("--__",rownames(tab1n_c)),]
    tab1n_c=tab1n_c[!grepl("--.__uncultured",rownames(tab1n_c)),]
    tab1n_c=tab1n_c[!grepl("--.__uncultured.bacterium",rownames(tab1n_c)),]
    tab1n_c=tab1n_c[!grepl("--.__gut.metagenome",rownames(tab1n_c)),]
    rownames(tab1n_c)=sapply(strsplit(rownames(tab1n_c),"--"),"[[",level_n[match(taxa_level,level1)])
  }

  if(num_taxa>nrow(tab1n_c)){
    if(all(colSums(tab1n_c)==100)){
      num_taxa=nrow(tab1n_c)
    }else{
      num_taxa=nrow(tab1n_c)
      other=100-colSums(tab1n_c)
      tab1n_c2=rbind(tab1n_c,other)
      rownames(tab1n_c2)[num_taxa+1]="Other"
    }
  }else{
    tab1n_c1=tab1n_c[1:num_taxa,]
    other=100-colSums(tab1n_c1)
    tab1n_c2=rbind(tab1n_c1,other)
    rownames(tab1n_c2)[num_taxa+1]="Other"
  }
  return(tab1n_c2)
}
