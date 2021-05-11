library(ggplot2)
library(ggpubr)
library(ggtree)
library(vegan)

#mutliple taxonomic tables, for example, with one for each level
taxa_table1=format_asv(taxa_file = "./data-raw/multiple_tsv")

#mutliple taxonomic tables in biom formats, for example, with one for each level
taxa_table2=format_asv(taxa_file = "./data-raw/multiple_biom",biom=T)

#ASV table (biom) from DADA2 with the taxonomy listed as the last column
taxa_table3=format_asv(taxa_file = "./data-raw/biom_taxonomy.biom",biom=T,onefile = T,ASV=T)

#ASV table (text) from DADA2 with the taxonomy listed as the last column
taxa_table4=format_asv(taxa_file = "./data-raw/table_taxonomy.txt",biom=F,onefile = T,ASV=T)

#One taxonomic table (test) with all levels
taxa_table5=format_asv(taxa_file = "./data-raw/taxa_all.csv",sep=",",biom=F,onefile = T,ASV=F)
taxa_table6=format_asv(taxa_file = "./data-raw/taxa_all.txt",sep="\t",biom=F,onefile = T,ASV=F)

#One taxonomic table (biom) with all levels
taxa_table7=format_asv(taxa_file = "./data-raw/table.from_txt_hdf5.biom",biom=T,onefile = T,ASV=F)
taxa_table8=format_asv(taxa_file = "./data-raw/table.from_txt_json.biom",biom=T,onefile = T,ASV=F)

taxa_table="./data-raw/biom_taxonomy.biom"
metadata_dir="./data-raw/metadata_cafe.csv"
taxa_table="./data-raw/multiple_biom"
#format the raw taxonomic abudance table
taxa_tab1=format_asv(taxa_file = taxa_table,biom=T,onefile = T,ASV=T)
#format metadata
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=2)
#subset the abundance table to only include samples for test
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="Study",stratify_by_value="Sugar")
#perform statistical test
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="Timepoint",method="wilcoxon")
#tree plot
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="Timepoint",fdr_cutoff=0.1)
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="Timepoint",fdr_cutoff=0.05)
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="Timepoint",fdr_cutoff=0.01)

#cor_plot
cor_plot1=meta_corplot(taxa_table =tab_s, metadata=metadata1,test_metadata="test_score",col_metadata="Timepoint",fdr_cutoff=0.3)

#Metaphlan2 and Kraken2 results
taxa_table9=format_wgs(taxa_file = "./data-raw/mali_kraken2.txt")
taxa_table10=format_wgs(taxa_file = "./data-raw/mali_phlan2.txt")

metadata_dir="./data-raw/metadata_mali.csv"
taxa_table="./data-raw/mali_phlan2.txt"
taxa_table="./data-raw/mali_kraken2.txt"
#format the raw taxonomic abudance table
taxa_tab1=format_wgs(taxa_file = taxa_table)
#format metadata
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=1)
#subset the abundance table to only include samples for test
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="time",stratify_by_value="1")
#perform statistical test
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="group",method="wilcoxon")
#tree plot
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="group",fdr_cutoff=0.1)
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="group",fdr_cutoff=0.1,taxa_removal="Candidatus_Saccharibacteria")
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="group",fdr_cutoff=0.05)
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="group",fdr_cutoff=0.01)

#PCoA plot
mds_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")

#alpha diversity boxplot
alpha_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",palette_group=c("red","blue","orange","green"))

#taxa boxplot
taxa_boxplot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",fdrs=fdrs1,log_norm=T,cutoff=0.01,palette_group=c("red","blue","orange","green"))



metadata_dir="/Users/shansun/Google\ Drive/mc_set1/metadata_parsed.csv"
taxa_dir1="/Users/shansun/Google\ Drive/mc_set1/set1/biom-with-taxonomy.biom"
taxa_dir2="/Users/shansun/Google\ Drive/mc_set1/set2/dada2/biom-with-taxonomy.biom"

#format the raw taxonomic abudance table
taxa_tab1=format_asv(taxa_file = taxa_dir1,biom=T,onefile = T,ASV=T)
taxa_tab2=format_asv(taxa_file = taxa_dir2,biom=T,onefile = T,ASV=T)
#overlap "M14705_DES" "M14912_DES"
taxa_tab2=taxa_tab2[,!grepl("M14705_DES",colnames(taxa_tab2))]
taxa_tab2=taxa_tab2[,!grepl("M14912_DES",colnames(taxa_tab2))]

colnames(taxa_tab1)=gsub("\\.","",colnames(taxa_tab1))
colnames(taxa_tab1)=paste(substr(colnames(taxa_tab1),1,6),substr(colnames(taxa_tab1),7,9),sep="_")

taxa_tab3=merge(taxa_tab1,taxa_tab2,by=0,all=T)
rownames(taxa_tab3)=taxa_tab3[,1]
taxa_tab3=taxa_tab3[,-1]
taxa_tab3[is.na(taxa_tab3)]=0
taxa_tab3=t(t(taxa_tab3)/colSums(taxa_tab3))*mean(colSums(taxa_tab3))

#format metadata
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=1)
metadata2=metadata1

metadata1$location="ASC"
rownames(metadata1)=paste(rownames(metadata1),"ASC",sep="_")
metadata2$location="DES"
rownames(metadata2)=paste(rownames(metadata2),"DES",sep="_")
metadata3=rbind(metadata1,metadata2)

metadata3$batch=rep(NA,nrow(metadata3))
metadata3$batch[which(rownames(metadata3)%in%colnames(taxa_tab1))]=1
metadata3$batch[which(rownames(metadata3)%in%colnames(taxa_tab2))]=2

write.csv(metadata3,"/Users/shansun/Google\ Drive/mc_set1/metadata_combined.csv")
write.csv(taxa_tab3,"/Users/shansun/Google\ Drive/mc_set1/taxa_combined.csv")


metadata_dir="/Users/shansun/Google\ Drive/mc_set1/metadata_combined.csv"
taxa_dir="/Users/shansun/Google\ Drive/mc_set1/taxa_combined.csv"

#format data
taxa_tab1=format_asv(taxa_file = taxa_dir,biom=F,onefile = T,ASV=F,sep=",")
#format metadata
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=1)

#subset the abundance table to only include samples for test
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="",stratify_by_value="",prevalence_cutoff=0.1, abundance_cutoff=0)

#perform statistical test
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="Case_Ctrl",method="wilcoxon")

#tree plot
plot1=tree_view(taxa_table =taxa_tab1, metadata=metadata1,fdrs=fdrs1,test_metadata="Case_Ctrl",prevalence_cutoff=0.1, abundance_cutoff=0)

#PCoA plot
mds_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="Case_Ctrl",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")

#alpha diversity boxplot
alpha_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="Case_Ctrl",palette_group=c("red","blue","orange","green"))

#taxa boxplot
taxa_boxplot(taxa_table = tab_s, metadata=metadata1,test_metadata="Case_Ctrl",fdrs=fdrs1,log_norm=T,cutoff=0.01,page=1,palette_group=c("red","blue","orange","green"))

#cor_plot
cor_plot1=meta_corplot(taxa_table =tab_s, metadata=metadata1,test_metadata="patient_age",col_metadata="Case_Ctrl",fdr_cutoff=0.1)

map1=read.table(file="/Users/shansun/Google\ Drive/mc_set1/set1/metadata.txt",sep="\t",quote="", header=T,row.names=1)
a1=metadata1[intersect(rownames(map1),rownames(metadata1)),]
a2=map1[intersect(rownames(map1),rownames(metadata1)),]
a3=cbind(a1[,3],a2[,6])

tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="batch",stratify_by_value="1",prevalence_cutoff=0.25, abundance_cutoff=0)
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="Case_Ctrl",method="wilcoxon")

tab_s2=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="batch",stratify_by_value="2",prevalence_cutoff=0.25, abundance_cutoff=0)
fdrs2=stat_test(taxa_table =tab_s2,metadata=metadata1,test_metadata="Case_Ctrl",method="wilcoxon")

t1=data.frame(sapply(by(t(tab_s),metadata1$Case_Ctrl[which(metadata1$batch==1)],colMeans),identity))
t2=data.frame(sapply(by(t(tab_s2),metadata1$Case_Ctrl[which(metadata1$batch==2)],colMeans),identity))

a1=sign(t1[,1]-t1[,2])
a2=sign(t2[,1]-t2[,2])
pdf("/Users/shansun/Google\ Drive/mc_set1/batch_comparison.pdf")
par(mfrow=c(1,1))
plot(-log10(fdrs1[,1])*a1,-log10(fdrs1[,2])*a1,xlab="batch1 log10(P)*direction",ylab="batch2 log10(P)*direction",cex.lab=1.2,cex.axis=1.2)
dev.off()



rsconnect::deployApp('/Users/shansun/git/plotmicrobiome')
