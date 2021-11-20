#' Statistical tests
#' wilcoxon, t-test, anova, kruskal-wallis
#' @param
#' @keywords statistical tests
#' @export
#' @examples
#'
library(nlme)
stat_test = function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,test_metadata_continuous=F,lm_anova=F,model_lm=NULL,random_effect_var=NULL,method="wilcoxon") {

  tab_s=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  map_s=metadata[intersect(colnames(taxa_table),rownames(metadata)),]
  if(!test_metadata_continuous){
    map_s[,test_metadata]=factor(as.character(map_s[,test_metadata]))
  }
  #change names to avoid special characters in linear models
  tab_s1=tab_s
  rownames(tab_s1)=paste0("a",c(1:nrow(tab_s1)))
  tabMeta=cbind(t(tab_s1),map_s)

  plm=matrix(nrow=nrow(tab_s),ncol=4)
  for (n in 1:nrow(tab_s)){

    if (method == "wilcoxon"){
      plm[n,1]=try(wilcox.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$statistic)
      plm[n,2]=try(wilcox.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
      est1=t.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$estimate
      if(is.na(plm[n,1])){
        plm[n,4]=NA
      }else if((est1[1]-est1[2])>0){
        plm[n,4]=gsub("mean in group ","",names(est1)[1])
      }else{
        plm[n,4]=gsub("mean in group ","",names(est1)[2])
      }
    }else if (method == "t.test"){
      plm[n,1]=try(t.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$statistic)
      plm[n,2]=try(t.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
      est1=t.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$estimate
      if(is.na(plm[n,1])){
        plm[n,4]=NA
      }else if((est1[1]-est1[2])>0){
        plm[n,4]=gsub("mean in group ","",names(est1)[1])
      }else{
        plm[n,4]=gsub("mean in group ","",names(est1)[2])
      }
    }else if (method == "anova"){
      plm[n,1]=try(summary(aov(as.numeric(tab_s[n,])~map_s[[test_metadata]]))[[1]][1,4])
      plm[n,2]=try(summary(aov(as.numeric(tab_s[n,])~map_s[[test_metadata]]))[[1]][1,5])
      plm[n,4]=NA
    }else if (method == "kruskal-wallis"){
      plm[n,1]=try(kruskal.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$statistic)
      plm[n,2]=try(kruskal.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
      plm[n,4]=NA
    }else if (method == "pearson"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]])$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]])$p.value)
      if(is.na(plm[n,1])){
        plm[n,4]=NA
      }else if(plm[n,1]>0){
        plm[n,4]="positive"
      }else{
        plm[n,4]="negative"
      }
    }else if (method == "spearman"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]],method="spearman")$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]],method="spearman")$p.value)

      if(is.na(plm[n,1])){
        plm[n,4]=NA
      }else if(plm[n,1]>0){
        plm[n,4]="positive"
      }else{
        plm[n,4]="negative"
      }
    }else if (method == "kendall"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]],method="kendall")$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),map_s[[test_metadata]],method="kendall")$p.value)
      if(is.na(plm[n,1])){
        plm[n,4]=NA
      }else if(plm[n,1]>0){
        plm[n,4]="positive"
      }else{
        plm[n,4]="negative"
      }
    }else if (method == "lm"){
      model = as.formula(paste(rownames(tab_s1)[n],"~",paste(test_metadata,model_lm,sep="+")))
      simpleMod = try(lm(model,data=tabMeta))

      if(class( simpleMod )=="try-error"){
        plm[n,1]  =NA
        plm[n,2]  =NA
      }
      if(lm_anova){
        plm[n,1] = anova(simpleMod)[1,4]
        plm[n,2] = anova(simpleMod)[1,5]
        plm[n,4] = NA
      }else{
        plm[n,1] = summary(simpleMod)$coefficients[2,3]
        plm[n,2] = summary(simpleMod)$coefficients[2,4]
        if(test_metadata_continuous){
          if(plm[n,1]<0){
            plm[n,4]="negative"
          }else{
            plm[n,4]="positive"
          }
        }else{
          if(plm[n,1]<0){
            plm[n,4]=names(table(tabMeta[,test_metadata]))[2]
          }else{
            plm[n,4]=names(table(tabMeta[,test_metadata]))[1]
          }
        }
      }
    }else if (method == "lme"){
      model = as.formula(paste(rownames(tab_s1)[n],"~",paste(test_metadata,model_lm,sep="+")))
      random_e=as.formula(paste0("~1|",random_effect_var))
      mixedMod = try(lme(model,method="REML",random=random_e,data=tabMeta,na.action=na.exclude))
      if(class( mixedMod )=="try-error"){
        plm[n,1]  =NA
        plm[n,2]  =NA
      }
      if(lm_anova){
        plm[n,1] = anova(mixedMod)[2,3]
        plm[n,2] = anova(mixedMod)[2,4]
        plm[n,4] = NA
      }else{
        plm[n,1] = summary(mixedMod)$tTable[2,4]
        plm[n,2] = summary(mixedMod)$tTable[2,5]
        if(test_metadata_continuous){
          if(plm[n,1]<0){
            plm[n,4]="negative"
          }else{
            plm[n,4]="positive"
          }
        }else{
          if(plm[n,1]<0){
            plm[n,4]=names(table(tabMeta[,test_metadata]))[2]
          }else{
            plm[n,4]=names(table(tabMeta[,test_metadata]))[1]
          }
        }
      }
    }
    else{
      print("Please select method from wilcoxon,t.test,kruskal-wallis, anova,
            pearson, spearman, kendall,lm and lme")
    }
  }
  plm[,3]=p.adjust(plm[,2],method="fdr")
  colnames(plm)=c("stats","P","FDR","Enriched")
  rownames(plm)=rownames(tab_s)
  return(plm)
}


