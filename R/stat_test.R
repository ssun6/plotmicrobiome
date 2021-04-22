#' Statistical tests
#' wilcoxon, t-test, anova, kruskal-wallis
#' @param
#' @keywords statistical tests
#' @export
#' @examples
#'
stat_test <- function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,method="wilcoxon") {

  if (!method%in%c("wilcoxon","t.test","kruskal-wallis","anova")){
    print("Please make method is one of wilcoxon,t.test,kruskal-wallis and anova")
  }

  tab_s=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  map_s=metadata[intersect(colnames(taxa_table),rownames(metadata)),]
  metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
  metadata[,test_metadata]=factor(metadata[,test_metadata])

  plm=vector()
  for (n in 1:nrow(tab_s)){

    if (method == "wilcoxon"){
      a1=try(wilcox.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
    }else if (method == "t.test"){
      a1=try(t.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
    }else if (method == "anova"){
      a1=try(summary(aov(as.numeric(tab_s[n,])~map_s[[test_metadata]]))[[1]][1,5])
    }else{
      a1=kruskal.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value
    }

    plm[n]=a1
  }
  pvals=cbind(plm,p.adjust(plm,method="fdr"))
  colnames(pvals)=c("P","FDR")
  rownames(pvals)=rownames(tab_s)
  return(pvals)
}


