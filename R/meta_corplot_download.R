meta_corplot_download=function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,col_metadata=NULL,one_level=F,log_norm=T,fdr_cutoff=0.1,cor_method="spearman",taxa_shown="",palette_group=c("red","blue","orange","green")){
  metadata=metadata[which(!is.na(metadata[[test_metadata]])),]
  if (col_metadata!=""){
    metadata=metadata[which(!is.na(metadata[[col_metadata]])),]
  }

  tab1=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  map1=metadata[intersect(colnames(tab1),rownames(metadata)),]

  cor_mat=matrix(nrow=nrow(tab1),ncol=2)
  for (j in 1:nrow(tab1)){
    cor1=cor.test(tab1[j,],map1[,test_metadata],method=cor_method)
    cor_mat[j,1]=cor1$estimate
    cor_mat[j,2]=cor1$p.value
  }
  cor_mat=data.frame(cor_mat)
  cor_mat$fdr=p.adjust(cor_mat[,2],method="fdr")
  rownames(cor_mat)=rownames(tab1)
  cor_mat=cor_mat[order(cor_mat[,3]),]
  cor_mat=na.omit(cor_mat)

  if(taxa_shown==""){
    tab1=tab1
    cor_mat=cor_mat
  }else{
    tab1=tab1[grep(taxa_shown,rownames(tab1)),]
    cor_mat=cor_mat[grep(taxa_shown,rownames(cor_mat)),]
  }

  k=1
  gplots1=list()
  if(nrow(cor_mat)==1){
    if (cor_mat[1,3]<fdr_cutoff){
      if(log_norm){
        tab_s=log10(tab1+1)
        xlab1="log10 (normalized abundance)"
      }else{
        tab_s=tab1
        xlab1="normalized abundance"
      }
      map1$i=tab_s
      colnames(map1)[match(test_metadata,colnames(map1))]="test_metadata"
      if(col_metadata!=""){
        colnames(map1)[match(col_metadata,colnames(map1))]="col_metadata"
      }

      if(one_level){
        tax_name1=rownames(cor_mat)
      }else{
        tax_name=paste0("p__",strsplit(rownames(cor_mat),"--p__")[[1]][2])
        if(nchar(tax_name)>60&nchar(tax_name)<100){
          tax_s=strsplit(tax_name,"--")[[1]]
          tax_l=length(tax_s)
          tax_li=round(tax_l/2)
          tax_name1=paste0(paste(tax_s[1:tax_li],collapse = ";"),"\n",paste(tax_s[(tax_li+1):tax_l],collapse = ";"))
        }else if(nchar(tax_name)>=100){
          tax_s=strsplit(tax_name,"--")[[1]]
          tax_l=length(tax_s)
          if(tax_l>5){
            tax_name1=paste0(paste(tax_s[1:3],collapse = ";"),"\n",paste(tax_s[4:5],collapse = ";"),"\n",paste(tax_s[6:tax_l],collapse = ";"))
          }else{
            tax_name1=paste0(paste(tax_s[1:2],collapse = ";"),"\n",paste(tax_s[3:4],collapse = ";"),"\n",paste(tax_s[5],collapse = ";"))
          }
        }else{
          tax_name1=tax_name
        }
      }

      if(cor_mat[1,2]<0.001){
        wil_p=formatC(cor_mat[1,2], format = "e", digits = 2)
      }else{
        wil_p=formatC(cor_mat[1,2], digits = 3)
      }

      if(cor_mat[1,3]<0.001){
        wil_fdr=formatC(cor_mat[1,3], format = "e", digits = 2)
      }else{
        wil_fdr=formatC(cor_mat[1,3], digits = 3)
      }

      main1=paste(tax_name1,"\n"," rho =",round(cor_mat[1,1],3),"\n P =", wil_p,"\n FDR =", wil_fdr,"\n")

      if (col_metadata!=""){
        g=ggscatter(map1, x = "i", y = "test_metadata",xlab = xlab1, ylab = test_metadata,
                    legend.title=col_metadata,font.x = c(10, "black"),font.y = c(10,  "black"), color = "col_metadata",palette = palette_group, size = 2,
                    add = "reg.line",add.params = list(color = "darkgrey", fill = "lightgray"),conf.int = TRUE,cor.coef = FALSE )
      }else{
        g=ggscatter(map1, x = "i", y = "test_metadata",xlab = xlab1, ylab = test_metadata,
                    font.x = c(10, "black"),font.y = c(10,  "black"),col = palette_group, size = 2,
                    add = "reg.line",add.params = list(color = "darkgrey", fill = "lightgray"),conf.int = TRUE,cor.coef = FALSE )
      }

      gplots1[[1]]=g+annotate(geom="text", x=min(map1$i)+sd(map1$i)*1.3, y=max(map1$test_metadata)-sd(map1$test_metadata)*0.5, label=main1,color="black",size=3)
    }
  }else{
    for (j in 1:nrow(cor_mat)){
      if (cor_mat[j,3]<fdr_cutoff){
        if(log_norm){
          tab_s=log10(tab1+1)
          xlab1="log10 (normalized abundance)"
        }else{
          tab_s=tab1
          xlab1="normalized abundance"
        }
        map1$i=tab_s[rownames(cor_mat)[j],]
        colnames(map1)[match(test_metadata,colnames(map1))]="test_metadata"
        if(col_metadata!=""){
          colnames(map1)[match(col_metadata,colnames(map1))]="col_metadata"
        }

        if(one_level){
          tax_name1=rownames(cor_mat)[j]
        }else{
          tax_name=paste0("p__",strsplit(rownames(cor_mat)[j],"--p__")[[1]][2])
          if(nchar(tax_name)>60&nchar(tax_name)<100){
            tax_s=strsplit(tax_name,"--")[[1]]
            tax_l=length(tax_s)
            tax_li=round(tax_l/2)
            tax_name1=paste0(paste(tax_s[1:tax_li],collapse = ";"),"\n",paste(tax_s[(tax_li+1):tax_l],collapse = ";"))
          }else if(nchar(tax_name)>=100){
            tax_s=strsplit(tax_name,"--")[[1]]
            tax_l=length(tax_s)
            if(tax_l>5){
              tax_name1=paste0(paste(tax_s[1:3],collapse = ";"),"\n",paste(tax_s[4:5],collapse = ";"),"\n",paste(tax_s[6:tax_l],collapse = ";"))
            }else{
              tax_name1=paste0(paste(tax_s[1:2],collapse = ";"),"\n",paste(tax_s[3:4],collapse = ";"),"\n",paste(tax_s[5:tax_l],collapse = ";"))
            }
          }else{
            tax_name1=tax_name
          }
        }

        if(cor_mat[1,2]<0.001){
          wil_p=formatC(cor_mat[1,2], format = "e", digits = 2)
        }else{
          wil_p=formatC(cor_mat[1,2], digits = 3)
        }

        if(cor_mat[1,3]<0.001){
          wil_fdr=formatC(cor_mat[1,3], format = "e", digits = 2)
        }else{
          wil_fdr=formatC(cor_mat[1,3], digits = 3)
        }

        main1=paste(tax_name1,"\n"," rho =",round(cor_mat[1,1],3),"\n P =", wil_p,"\n FDR =", wil_fdr,"\n")

        if (col_metadata!=""){
          g=ggscatter(map1, x = "i", y = "test_metadata",xlab = xlab1, ylab = test_metadata,
                      legend.title=col_metadata,font.x = c(10, "black"),font.y = c(10,  "black"), color = "col_metadata",palette = palette_group, size = 2,
                      add = "reg.line",add.params = list(color = "darkgrey", fill = "lightgray"),conf.int = TRUE,cor.coef = FALSE )
        }else{
          g=ggscatter(map1, x = "i", y = "test_metadata",xlab = xlab1, ylab = test_metadata,
                      font.x = c(10, "black"),font.y = c(10,  "black"),col = palette_group, size = 2,
                      add = "reg.line",add.params = list(color = "darkgrey", fill = "lightgray"),conf.int = TRUE,cor.coef = FALSE )
        }

        gplots1[[k]]=g+annotate(geom="text", x=min(map1$i)+sd(map1$i)*1.3, y=max(map1$test_metadata)-sd(map1$test_metadata)*0.5, label=main1,color="black",size=3)

        k=k+1
      }
    }
  }

  if (length(gplots1)==0){
    message("No Significant hits! Try to increase the FDR cutoff.")
  }
  ggarrange(plotlist=gplots1,ncol = 2, nrow = 2)
}
