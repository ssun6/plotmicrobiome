#' MDS plots and PERMANOVA tests
#' @keywords PCoA, NMDS, PERMANOVA
#' @export
#' @examples
#'
#'
mds_plot=function(taxa_table = NULL, metadata=NULL,test_metadata=NULL,taxa_level="genus",method_mds="pcoa",xaxis=1,yaxis=2,distance_type="bray",palette_group=c("red","blue","orange","green")){

  metadata=metadata[which(!is.na(metadata[,test_metadata])),]
  tab1=taxa_table[,intersect(colnames(taxa_table),rownames(metadata))]
  metadata=metadata[intersect(colnames(tab1),rownames(metadata)),]
  metadata[,test_metadata]=factor(as.character(metadata[,test_metadata]))
  metadata[,test_metadata]=droplevels(metadata[,test_metadata])

  tax_l=sapply(strsplit(rownames(tab1),"--"),function(i){length(i)})
  level1=c("kingdom","phylum","class","order","family","genus","species","strain")
  level_n=c(1:8)
  tab1n=tab1[which(tax_l==level_n[match(taxa_level,level1)]),]

  ado_p1=as.numeric(unlist(vegan::adonis(t(tab1n)~metadata[,test_metadata])$"aov.tab"[1,5:6]))

  par(mfrow=c(1,1),mar=c(5,5,5,5))
  if(method_mds=="pcoa"){
    gen_mds=vegan::capscale(t(tab1n)~1,distance=distance_type)
    var_per=round((gen_mds$CA$eig/sum(gen_mds$CA$eig))[1:6]*100,2)
    mds_p=paste("PCoA",c(1:6)," (",var_per,"%)",sep="")
    main1=paste(test_metadata,"R2 =",round(ado_p1[1],5)," P=",ado_p1[2])
    pcoa12=vegan::ordiplot(gen_mds,choices=c(xaxis,yaxis),type="none",cex.lab=1.5,xlab=mds_p[1],ylab=mds_p[2],main=main1,xlim=c(min(summary(gen_mds)$sites[,1])-0.3,max(summary(gen_mds)$sites[,1])+0.3))
  }else if (method_mds=="nmds"){
    gen_mds=vegan::metaMDS(t(tab1n),distance=distance_type)
    mds_p=paste("NMDS",c(1:6),sep="")
    main1=paste(test_metadata,"R2 =",round(ado_p1[1],5)," P=",ado_p1[2],"\nstress =",round(gen_mds$grstress,5))
    pcoa12=vegan::ordiplot(gen_mds,type="none",cex.lab=1.5,xlab=mds_p[1],ylab=mds_p[2],main=main1,xlim=c(min(gen_mds$points[,1])-0.3,max(gen_mds$points[,1])+0.3))
  }else{
    stop("Please use pcoa or nmds for method_mds")
  }

  col3=palette_group[match(levels(metadata[,test_metadata]),levels(metadata[,test_metadata]))][factor(metadata[,test_metadata])]
  pch1=16
  points(pcoa12,"sites",col=adjustcolor(col3, alpha.f = 0.3),pch=pch1,cex=1.5)
  for (j in 1:length(levels(metadata[,test_metadata]))){
    ordiellipse(pcoa12, metadata[,test_metadata], kind="se", conf=0.95, lwd=1, draw = "lines", col=palette_group[match(levels(metadata[,test_metadata]),levels(metadata[,test_metadata]))][j],show.groups=levels(metadata[,test_metadata])[j],label=T,font=2,cex=1.3)
  }
  legend("topright",levels(factor(metadata[,test_metadata])), cex=1.2, bty="n", col=palette_group[match(levels(metadata[,test_metadata]),levels(metadata[,test_metadata]))], pch=16)
}
