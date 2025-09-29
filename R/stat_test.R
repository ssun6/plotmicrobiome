#' Statistical tests
#' wilcoxon, t-test, anova, kruskal-wallis
#' @param
#' @keywords statistical tests
#' @export
#' @examples
#'
library(nlme)
stat_test = function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,test_metadata_continuous=F,log_norm=T,glm_anova=F,outcome_meta=F,model_glm="",glm_dist="default",glm_ref="",method="wilcoxon",random_effect_var="") {

  tab_s=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  map_s=metadata[intersect(colnames(taxa_table),rownames(metadata)),]
  if(!test_metadata_continuous){
    map_s[,test_metadata]=factor(as.character(map_s[,test_metadata]))
  }
  #change names to avoid special characters in linear models
  if(log_norm){
    tab_s=log10(tab_s+1)
  }
  tab_s1=tab_s
  rownames(tab_s1)=paste0("a",c(1:nrow(tab_s1)))
  tabMeta=cbind(t(tab_s1),map_s)

  plm=matrix(nrow=nrow(tab_s),ncol=5)
  for (n in 1:nrow(tab_s)){

    if (method == "wilcoxon"){
      plm[n,1]=try(wilcox.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$statistic)
      plm[n,2]=try(wilcox.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
      est1=t.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$estimate
      plm[n,5]=est1[2]/est1[1]
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
      plm[n,5]=est1[2]/est1[1]
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
      plm[n,5]=NA
    }else if (method == "kruskal-wallis"){
      plm[n,1]=try(kruskal.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$statistic)
      plm[n,2]=try(kruskal.test(as.numeric(tab_s[n,])~map_s[[test_metadata]])$p.value)
      plm[n,4]=NA
      plm[n,5]=NA
    }else if (method == "pearson"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),as.numeric(as.character(map_s[[test_metadata]])))$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),as.numeric(as.character(map_s[[test_metadata]])))$p.value)
      plm[n,5]=NA
      if(is.na(plm[n,1])){
        plm[n,4]=NA
      }else if(plm[n,1]>0){
        plm[n,4]="positive"
      }else{
        plm[n,4]="negative"
      }
    }else if (method == "spearman"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),as.numeric(as.character(map_s[[test_metadata]])),method="spearman")$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),as.numeric(as.character(map_s[[test_metadata]])),method="spearman")$p.value)
      plm[n,5]=NA
      if(is.na(plm[n,1])){
        plm[n,4]=NA
      }else if(plm[n,1]>0){
        plm[n,4]="positive"
      }else{
        plm[n,4]="negative"
      }
    }else if (method == "kendall"){
      plm[n,1]=try(cor.test(as.numeric(tab_s[n,]),as.numeric(as.character(map_s[[test_metadata]])),method="kendall")$estimate)
      plm[n,2]=try(cor.test(as.numeric(tab_s[n,]),as.numeric(as.character(map_s[[test_metadata]])),method="kendall")$p.value)
      plm[n,5]=NA
      if(is.na(plm[n,1])){
        plm[n,4]=NA
      }else if(plm[n,1]>0){
        plm[n,4]="positive"
      }else{
        plm[n,4]="negative"
      }
    }else if (method == "lr"){
      if(model_glm!=""){
        model = as.formula(paste(test_metadata,"~",paste0(rownames(tab_s1)[n],model_glm)))
      }else{
        model = as.formula(paste(test_metadata,"~",rownames(tab_s1)[n]))
      }
      if (glm_ref!=""){
        tabMeta[[test_metadata]] <- relevel(tabMeta[[test_metadata]], ref = glm_ref)
      }else{
        tabMeta[[test_metadata]] <- relevel(tabMeta[[test_metadata]], ref = names(table(tabMeta[[test_metadata]]))[1])
      }
      if(glm_dist=="default"){
        glm_dist="binomial"
      }else{
        glm_dist=glm_dist
      }

      simpleMod = try(glm(model,data=tabMeta,family=glm_dist))
      plm[n,1] = summary(simpleMod)$coefficients[2,3]
      plm[n,2] = summary(simpleMod)$coefficients[2,4]
      if(plm[n,1]<0){
        plm[n,4]=names(table(tabMeta[,test_metadata]))[1]
      }else{
        plm[n,4]=names(table(tabMeta[,test_metadata]))[2]
      }
    }else if (method == "glm"){
      if(outcome_meta){
        if(model_glm!=""){
          model = as.formula(paste(test_metadata,"~",paste0(rownames(tab_s1)[n],model_glm)))
        }else{
          model = as.formula(paste(test_metadata,"~",rownames(tab_s1)[n]))
        }
        if(test_metadata_continuous){
          if(glm_dist=="default"){
            glm_dist="gaussian"
          }else{
            glm_dist=glm_dist
          }
        }else{
          if (glm_ref!=""){
            tabMeta[[test_metadata]] <- relevel(tabMeta[[test_metadata]], ref = glm_ref)
          }else{
            tabMeta[[test_metadata]] <- relevel(tabMeta[[test_metadata]], ref = names(table(tabMeta[[test_metadata]]))[1])
          }
          if(glm_dist=="default"){
            glm_dist="binomial"
          }else{
            glm_dist=glm_dist
          }
        }
        simpleMod = try(glm(model,data=tabMeta,family=glm_dist))
        if(class( simpleMod )=="try-error"){
          plm[n,1]  =NA
          plm[n,2]  =NA
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
              plm[n,4]=names(table(tabMeta[,test_metadata]))[1]
            }else{
              plm[n,4]=names(table(tabMeta[,test_metadata]))[2]
            }
          }
        }
      }else{
        if(model_glm!=""){
          model = as.formula(paste0(rownames(tab_s1)[n],"~",paste(test_metadata,model_glm)))
        }else{
          model = as.formula(paste(rownames(tab_s1)[n],"~",test_metadata))
        }
        if (!test_metadata_continuous){
          if (glm_ref!=""){
            tabMeta[[test_metadata]] <- relevel(tabMeta[[test_metadata]], ref = glm_ref)
          }else{
            tabMeta[[test_metadata]] <- relevel(tabMeta[[test_metadata]], ref = names(table(tabMeta[[test_metadata]]))[1])
          }
        }
        if(glm_dist=="default"){
          glm_dist="gaussian"
        }else{
          glm_dist=glm_dist
        }
        simpleMod = try(glm(model,data=tabMeta,family=glm_dist))
        if(class( simpleMod )[1]=="try-error"){
          plm[n,1]  =NA
          plm[n,2]  =NA
        }else{
          if(glm_anova ){
            plm[n,1] = summary(aov(simpleMod))[[1]][1,4]
            plm[n,2] = summary(aov(simpleMod))[[1]][1,5]
            plm[n,4] = NA
          }else if (!glm_anova){
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
                plm[n,4]=names(table(tabMeta[,test_metadata]))[1]
              }else{
                plm[n,4]=names(table(tabMeta[,test_metadata]))[2]
              }
            }
          }
        }
      }
    }else if (method=="lme"){
      model = as.formula(paste(rownames(tab_s1)[n],"~",paste0(test_metadata,model_glm)))
      colnames(tabMeta)[match(random_effect_var,colnames(tabMeta))]="subject_ID"
      if (!test_metadata_continuous){
        if (!is.factor(tabMeta[[test_metadata]])){
          tabMeta[[test_metadata]]=factor(tabMeta[[test_metadata]])
        }
        if (glm_ref!=""){
          tabMeta[[test_metadata]] <- relevel(tabMeta[[test_metadata]], ref = glm_ref)
        }else{
          tabMeta[[test_metadata]] <- relevel(tabMeta[[test_metadata]], ref = names(table(tabMeta[[test_metadata]]))[1])
        }
      }
      mixedMod = try(lme(model,method="REML",random=~1|subject_ID,data=tabMeta,na.action=na.exclude),silent = TRUE)
      if(class( mixedMod )=="try-error"){
        plm[n,1]  =NA
        plm[n,2]  =NA
      }else{
        if(glm_anova){
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
              plm[n,4]=names(table(tabMeta[,test_metadata]))[1]
            }else{
              plm[n,4]=names(table(tabMeta[,test_metadata]))[2]
            }
          }
        }
      }
    }
    else{
      print("Please select method from wilcoxon,t.test,kruskal-wallis, anova,
            pearson, spearman, kendall,glm and lme")
    }
  }
  plm[,3]=p.adjust(plm[,2],method="fdr")
  colnames(plm)=c("stats","P","FDR","Enriched","Fold_change")
  rownames(plm)=rownames(tab_s)
  plm=data.frame(plm[order(as.numeric(plm[,2])),])
  return(plm)
}


