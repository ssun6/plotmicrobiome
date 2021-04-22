library(ggplot2)
library(ggpubr)
library(ggtree)
library(vegan)

#mutliple taxonomic tables, for example, with one for each level
taxa_table1=format_asv(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/multiple_tsv")

#mutliple taxonomic tables in biom formats, for example, with one for each level
taxa_table2=format_asv(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/multiple_biom",biom=T)

#ASV table (biom) from DADA2 with the taxonomy listed as the last column
taxa_table3=format_asv(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/biom_taxonomy.biom",biom=T,onefile = T,ASV=T)

#ASV table (text) from DADA2 with the taxonomy listed as the last column
taxa_table4=format_asv(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/table_taxonomy.txt",biom=F,onefile = T,ASV=T)

#One taxonomic table (test) with all levels
taxa_table5=format_asv(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/taxa_all.csv",sep=",",biom=F,onefile = T,ASV=F)
taxa_table6=format_asv(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/taxa_all.txt",sep="\t",biom=F,onefile = T,ASV=F)

#One taxonomic table (biom) with all levels
taxa_table7=format_asv(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/table.from_txt_hdf5.biom",biom=T,onefile = T,ASV=F)
taxa_table8=format_asv(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/table.from_txt_json.biom",biom=T,onefile = T,ASV=F)

taxa_table="/Users/shansun/git/plotmicrobiome/data-raw/biom_taxonomy.biom"
metadata_dir="/Users/shansun/git/plotmicrobiome/data-raw/metadata_cafe.csv"
taxa_table="/Users/shansun/git/plotmicrobiome/data-raw/multiple_biom"
#format the raw taxonomic abudance table
taxa_tab1=format_asv(taxa_file = taxa_table,biom=T)
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
taxa_table9=format_wgs(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/mali_kraken2.txt")
taxa_table10=format_wgs(taxa_file = "/Users/shansun/git/plotmicrobiome/data-raw/mali_phlan2.txt")

metadata_dir="/Users/shansun/git/plotmicrobiome/data-raw/metadata_mali.csv"
taxa_table="/Users/shansun/git/plotmicrobiome/data-raw/mali_phlan2.txt"
taxa_table="/Users/shansun/git/plotmicrobiome/data-raw/mali_kraken2.txt"
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


#subset the abundance table to only include samples for test
tab_s=table_subset(taxa_table = taxa_tab3,metadata=metadata3,stratify_by_metadata=NULL,prevalence_cutoff=0.25, abundance_cutoff=0)

#perform statistical test
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata3,test_metadata="Sex",method="wilcoxon")

#tree plot
plot1=tree_view(taxa_table =tab_s, metadata=metadata3,fdrs=fdrs1[,2],test_metadata="Sex",fdr_cutoff=0.1)

#PCoA plot
pcoa_plot(taxa_table = tab_s, metadata=metadata3,test_metadata="Sex",palette_group=c("red","blue","orange","green"),distance_type="bray")

#alpha diversity boxplot
alpha_plot(taxa_table = tab_s, metadata=metadata3,test_metadata="Sex",palette_group=c("red","blue","orange","green"))

#taxa boxplot
taxa_boxplot(taxa_table = tab_s, metadata=metadata3,test_metadata="Sex",fdrs=fdrs1[,2],log_norm=T,cutoff=0.01,palette_group=c("red","blue","orange","green"))
