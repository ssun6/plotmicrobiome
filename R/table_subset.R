#' Subset data
#' @param
#' @keywords subset
#' @export
#' @examples
#'
table_subset=function(taxa_table = NULL, metadata=NULL,stratify_by_metadata="",stratify_by_value="",prevalence_cutoff=0.1, abundance_cutoff=0,pathway=F) {

  inters_names=intersect(colnames(taxa_table),rownames(metadata))
  tab1=taxa_table[,match(inters_names,colnames(taxa_table))]
  metadata=metadata[match(inters_names,rownames(metadata)),]

  if(pathway){
    tab1=tab1
  }else{
    tab1=tab1[grepl("__Bacteria",rownames(tab1)) & !grepl("__Bacteria--__",rownames(tab1)),]
  }


  tab1_1=tab1
  if (prevalence_cutoff>0){
    tab1_1=tab1_1[which(apply(tab1_1,1,function(i){length(which(i!=0))})>=ncol(tab1)*prevalence_cutoff),]
  }

  if (abundance_cutoff>0){
    tab1_1=tab1_1[which(rowMeans(tab1_1)>abundance_cutoff),]
  }

  if(!stratify_by_metadata==""){
    tab_s=tab1_1[,which(metadata[,match(stratify_by_metadata,colnames(metadata))]==stratify_by_value)]
  }else{
    tab_s=tab1_1
  }
  return(tab_s)
}
