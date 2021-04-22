#' Parse metadata file
#' @param
#' @keywords format, meta
#' @export
#' @examples
#'
#'
meta_format=function(metadata=NULL,metadata_sep=",",meta_sample_name_col=0){
  if(metadata_sep==","){
    map=read.csv(file=metadata,header=T)
  }else{
    map=read.table(file=metadata,sep=metadata_sep,header=T)
  }
  #if rownames are not samples names and the column is specified,
  #change sample names in metadata to that column
  if(meta_sample_name_col>0){
    rownames(map)=map[,meta_sample_name_col]
  }
  map[map==""]=NA
  return(map)
}
