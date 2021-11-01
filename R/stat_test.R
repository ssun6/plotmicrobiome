#' Statistical tests
#' wilcoxon, t-test, anova, kruskal-wallis
#' @param
#' @keywords statistical tests
#' @export
#' @examples
#'
stat_test <- function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,method="wilcoxon") {

  tab_s=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  map_s=metadata[intersect(colnames(taxa_table),rownames(metadata)),]
  metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
  metadata[,test_metadata]=factor(metadata[,test_metadata])

  plm=matrix(nrow=nrow(tab_s),ncol=3)
  for (n in 1:nrow(tab_s)){

    if (method == "wilcoxon"){
      plm[n,1]=try(wilcox.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$statistic)
      plm[n,2]=try(wilcox.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
    }else if (method == "t.test"){
      plm[n,1]=try(t.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$statistic)
      plm[n,2]=try(t.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
    }else if (method == "anova"){
      plm[n,1]=try(summary(aov(as.numeric(tab_s[n,])~map_s[[test_metadata]]))[[1]][1,4])
      plm[n,2]=try(summary(aov(as.numeric(tab_s[n,])~map_s[[test_metadata]]))[[1]][1,5])
    }else if (method == "kruskal-wallis"){
      plm[n,1]=try(kruskal.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$statistic)
      plm[n,2]=try(kruskal.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
    }else if (method == "pearson"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]])$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]])$p.value)
    }else if (method == "spearman"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]],method="spearman")$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]],method="spearman")$p.value)
    }else if (method == "kendall"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]],method="kendall")$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]],method="kendall")$p.value)
    }else{
      print("Please select method from wilcoxon,t.test,kruskal-wallis, anova,
            pearson, spearman and kendall")
    }
  }
  plm[,3]=p.adjust(plm[,2],method="fdr")
  colnames(plm)=c("stats","P","FDR")
  rownames(plm)=rownames(tab_s)
  return(plm)
}


