p_compare=function(table1, table2,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=0.05,cor_method="spearman"){

  table_m=merge(table1,table2,by=0)
  rownames(table_m)=table_m[,1]
  table_m=table_m[,-1]
  table_m$logP1=-log10(as.numeric(table_m[,p_col1]))
  table_m$logP2=-log10(as.numeric(table_m[,ncol(table1)+p_col1]))
  levels2=names(table(table_m[,ncol(table1)+indicator2]))

  if(!indicator1==""){
    levels1=names(table(table_m[,indicator1]))
    levels1_1=rep(NA,nrow(table_m))
    levels1_1[which(table_m[,indicator1]==levels1[1])]=-1
    levels1_1[which(table_m[,indicator1]==levels1[2])]=1
    table_m$logP1=table_m$logP1*sign(levels1_1)
    if(is.na(sd(sign(levels1_1)))){
      lab1=NA
    }
    else if(sd(sign(na.omit(levels1_1)))==0){
      lab1="data1 -log10(P)"
    }else{
      lab1=paste(levels1[1],"     data1 -log10(P)*direction     ",levels1[2])
    }
  }else{
    lab1="data1 -log10(P)"
  }

  if(!indicator2==""){
    levels2=names(table(table_m[,indicator2+ncol(table1)]))
    levels2_1=rep(NA,nrow(table_m))
    levels2_1[which(table_m[,indicator2+ncol(table1)]==levels2[1])]=-1
    levels2_1[which(table_m[,indicator2+ncol(table1)]==levels2[2])]=1
    table_m$logP2=table_m$logP2*sign(levels2_1)

    if(is.na(sd(sign(levels2_1)))){
      lab1=NA
    }
    else if(sd(sign(levels2_1))==0){
      lab2="data2 -log10(P)"
    }else{
      lab2=paste(levels2[1],"     data2 -log10(P)*direction     ",levels2[2])
    }
  }else{
    lab2="data2 -log10(P)"
  }

  cor1=cor.test(table_m$logP1, table_m$logP2,method=cor_method)
  if(cor1$p.value<0.001){
    cor_p=formatC(cor1$p.value, format = "e", digits = 2)
  }else{
    cor_p=formatC(cor1$p.value, digits = 2)
  }
  main1=paste("Cor =",round(cor1$estimate,3),"P =",cor_p)

  table_m=data.frame(table_m)
  p=ggplot(table_m, mapping=aes_string(x="logP1", y="logP2")) +geom_point(color = point_color)+ theme_classic(base_size = 15) + labs(title=main1,x =lab1 , y =lab2)
  tax_lab=rownames(table_m)
  tax_lab1=sapply(strsplit(tax_lab,"--"),function(i){i[length(i)]})
  tax_lab2=sapply(strsplit(tax_lab,"--"),function(i){i[length(i)-1]})
  tax_lab1[which(tax_lab1=="__")]=paste0(tax_lab2[which(tax_lab1=="__")],"__unclassified")

  lab_cutoff1=-log10(lab_cutoff)
  tax_lab1[which(abs(table_m$logP1)<lab_cutoff1 & abs(table_m$logP2)<lab_cutoff1)]=NA
  p2=p+geom_text_repel(aes(label =tax_lab1),size = 3.5)+geom_vline(xintercept=0, linetype="dotted")+geom_hline(yintercept=0, linetype="dotted")+geom_abline(intercept = 0, slope = 1, linetype="dotted")
  print(p2)
}
