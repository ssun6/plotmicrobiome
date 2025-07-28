#' MDS plots and PERMANOVA tests
#' @keywords PCoA, NMDS, PERMANOVA
#' @export
#' @examples
#'
#'
mds_stats=function(taxa_table = NULL, one_level=F, metadata=NULL,test_metadata=NULL,log_norm=F,taxa_level="genus",method_mds="pcoa",distance_type="bray"){
  
  metadata=metadata[which(!is.na(metadata[,test_metadata])),]
  tab1=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  metadata=metadata[intersect(colnames(tab1),rownames(metadata)),]
  metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
  metadata[,test_metadata]=droplevels(metadata[,test_metadata])
  
  tab1=tab1[which(rowSums(tab1)!=0),]
  
  if(log_norm){
    tab1=log10(tab1+1)
  }
  
  results1=list()
  if (one_level){
    ado_p1=as.numeric(unlist(vegan::adonis2(t(tab1)~metadata[,test_metadata],method=distance_type)[1,c(3,5)]))
    results1[[1]]=ado_p1
    if(method_mds=="pcoa"){
      gen_mds=vegan::capscale(t(tab1n)~1,distance=distance_type)
      sample_coords <- vegan::scores(gen_mds,choices=c(1:6), display = "sites")
      results1[[2]]=sample_coords
    }else if (method_mds=="nmds"){
      par(mfrow=c(1,1),mar=c(5,5,5,5))
      gen_mds=vegan::metaMDS(t(tab1n),distance=distance_type)
      sample_coords <- vegan::scores(gen_mds, display = "sites")
      results1[[2]]=sample_coords
    }else{
      stop("Please use pcoa or nmds for method_mds")
    }
  }else{
    tax_l=sapply(strsplit(rownames(tab1),"--"),function(i){length(i)})
    level1=c("kingdom","phylum","class","order","family","genus","species","ASV_or_strain")
    level_n=c(1:8)
    if (taxa_level=="strain" | taxa_level=="ASV"){
      tab1n=tab1[which(tax_l==8),]
    }else{
      tab1n=tab1[which(tax_l==level_n[match(taxa_level,level1)]),]
    }
    
    tab1n=tab1n[which(rowSums(tab1n)!=0),]
    ado_p1=as.numeric(unlist(vegan::adonis2(t(tab1n)~metadata[,test_metadata],method=distance_type)[1,c(3,5)]))
    results1[[1]]=ado_p1
    if(method_mds=="pcoa"){
      gen_mds=vegan::capscale(t(tab1n)~1,distance=distance_type)
      sample_coords <- vegan::scores(gen_mds,choices=c(1:6), display = "sites")
      results1[[2]]=sample_coords
    }else if (method_mds=="nmds"){
      par(mfrow=c(1,1),mar=c(5,5,5,5))
      gen_mds=vegan::metaMDS(t(tab1n),distance=distance_type)
      sample_coords <- vegan::scores(gen_mds, display = "sites")
      results1[[2]]=sample_coords
    }else{
      stop("Please use pcoa or nmds for method_mds")
    }
  }
  return(results1)
}
