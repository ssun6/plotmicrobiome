taxa_boxplot_download=function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,test_metadata_order="default",fdrs=NULL,log_norm=T,cutoff=0.1,xlab_direction=1,palette_group=c("red","blue","orange","green"),taxa_shown="",one_level=F,xlab="default",ylab="default"){
  tab=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  metadata=metadata[match(intersect(colnames(tab),rownames(metadata)),rownames(metadata)),]
  if(is.factor(metadata[,test_metadata])){
    metadata[,test_metadata]=droplevels(metadata[,test_metadata])
  }else{
    metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
    metadata[,test_metadata]=droplevels(metadata[,test_metadata])
  }

  if(test_metadata_order[1]!="default"){
    metadata[,test_metadata]=factor(metadata[,test_metadata],levels=test_metadata_order)
  }

  inter_tax=intersect(rownames(tab),rownames(fdrs))
  tab=tab[match(inter_tax,rownames(tab)),]
  fdrs=fdrs[match(inter_tax,rownames(fdrs)),]

  if(taxa_shown==""){
    tab1=tab
    fdrs1=fdrs
  }else{
    tab1=tab[grep(taxa_shown,rownames(tab)),]
    fdrs1=fdrs[grep(taxa_shown,rownames(fdrs)),]
  }

  sig_l=length(which(as.numeric(fdrs1[,3])<cutoff))
  tab1n=tab1[order(as.numeric(fdrs1[,2])),][1:sig_l,]
  pval_s=as.numeric(fdrs1[,2][order(as.numeric(fdrs1[,2]))][1:sig_l])
  fdrs_s=as.numeric(fdrs1[,3][order(as.numeric(fdrs1[,2]))][1:sig_l])

  names1=rownames(tab1[order(as.numeric(fdrs1[,2])),])[1:sig_l]

  if(log_norm){
    tab1n=log10(tab1n+1)
  }

  par(mfrow=c(3,3),mar=c(10,5,5,5))
  for (i in 1:nrow(tab1n)){
    if (one_level){
      tax_name=names1[i]

      if(nchar(tax_name)>50&nchar(tax_name)<100){
        tax_s=strsplit(tax_name," ")[[1]]
        tax_l=length(tax_s)
        if(tax_l!=1){
          tax_li=round(tax_l/2)
          tax_name1=paste0(paste(tax_s[1:tax_li],collapse = " "),"\n",paste(tax_s[(tax_li+1):tax_l],collapse = " "))
        }
      }else if(nchar(tax_name)>=100){
        tax_s=strsplit(tax_name," ")[[1]]
        tax_l=length(tax_s)
        if(tax_l!=1){
          tax_l1=round(tax_l/3)
          tax_l2=round(tax_l/3*2)
          tax_name1=paste0(paste(tax_s[1:tax_l1],collapse = " "),"\n",paste(tax_s[(tax_l1+1):tax_l2],collapse = " "),"\n",paste(tax_s[(tax_l2+1):tax_l],collapse = " "))
        }
      }else{
        tax_name1=tax_name
      }
    }else{
      #tax_name=paste0("p__",strsplit(names1[i],"--p__")[[1]][2])
      tax_name=names1[i]
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

    if(pval_s[i]<0.001){
      wil_p=formatC(as.numeric(pval_s[i]), format = "e", digits = 2)
    }else{
      wil_p=formatC(as.numeric(pval_s[i]), digits = 2)
    }

    if(fdrs_s[i]<0.001){
      wil_fdr=formatC(as.numeric(fdrs_s[i]), format = "e", digits = 2)
    }else{
      wil_fdr=formatC(as.numeric(fdrs_s[i]), digits = 2)
    }

    if(xlab!="default"){
      xlab1=xlab
    }else{
      xlab1=test_metadata
    }

    if(ylab!="default"){
      ylab1=ylab
    }else{
      ylab1="normalized abundance"
    }

    if(sig_l==1){
      boxplot(as.numeric(tab1n)~metadata[,test_metadata],main=paste(tax_name1,"\nP =",wil_p,"\nFDR =",wil_fdr),border=palette_group,col="white",xlab=xlab1,ylab=ylab1,cex.main=0.8,las=xlab_direction)
      stripchart(as.numeric(tab1n)~metadata[,test_metadata],vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = palette_group)
    }else{
      boxplot(as.numeric(tab1n[i,])~metadata[,test_metadata],main=paste(tax_name1,"\nP =",wil_p,"\nFDR =",wil_fdr),border=palette_group,col="white",xlab=xlab1,ylab=ylab1,cex.main=0.8,las=xlab_direction)
      stripchart(as.numeric(tab1n[i,])~metadata[,test_metadata],vertical = TRUE,  method = "jitter", add = TRUE, pch = 16, col = palette_group)
    }
  }
}
