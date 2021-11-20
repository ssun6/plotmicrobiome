#test
devtools::install_github("ssun6/plotmicrobiome")
setwd("/Users/shansun/git/plotmicrobiome")

library(plotmicrobiome)
#mutliple taxonomic tables, for example, with one for each level
taxa_table1=format_asv(taxa_file = "./data-raw/multiple_tsv",biom=F,onefile = F,ASV=F)

#mutliple taxonomic tables in biom formats, for example, with one for each level
taxa_table2=format_asv(taxa_file = "./data-raw/multiple_biom",biom=T,onefile = F,ASV=F)

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
#format the raw taxonomic abudance table
taxa_tab1=format_asv(taxa_file = taxa_table,biom=T,onefile = T,ASV=T)
#format metadata
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=2)

#subset the abundance table to only include samples for test
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="Study",stratify_by_value="Sugar")
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="Study",stratify_by_value="Sugar",exclude_ASV=T)

#PCoA plot
mds_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="Treatment",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")

#alpha diversity boxplot
alpha_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="Treatment",test_metadata_order=c("CTL","SUG"),palette_group=c("red","blue","orange","green"))

#taxa_barplot
taxa_barplot(taxa_table = tab_s, metadata=metadata1,test_metadata="Treatment",num_taxa=10,taxa_level="phylum",xlab_direction=1,legend_size=1)
taxa_barplot(taxa_table = tab_s, metadata=metadata1,test_metadata="Treatment",num_taxa=10,taxa_level="genus",xlab_direction=1,legend_size=1)

#taxa boxplot
taxa_boxplot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",fdrs=fdrs1,log_norm=T,cutoff=0.01,palette_group=c("red","blue","orange","green"))

#subset the abundance table to only include samples for test
tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,stratify_by_metadata="Study",stratify_by_value="Sugar")

#perform statistical test for Timepoint
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="Timepoint",method="wilcoxon")
#tree plot
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="Timepoint",fdr_cutoff=0.1)
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="Timepoint",fdr_cutoff=0.05)
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="Timepoint",fdr_cutoff=0.01)

#perform correlation test for test_score
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="test_score",method="spearman")
#tree plot
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="test_score",test_metadata_continuous=T,fdr_cutoff=0.3)
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="test_score",test_metadata_continuous=T,fdr_cutoff=0.4)

#cor_plot
cor_plot1=meta_corplot(taxa_table =tab_s, metadata=metadata1,test_metadata="test_score",col_metadata="Timepoint",fdr_cutoff=0.3)

#Metaphlan2 and Kraken2 results
taxa_table9=format_wgs(taxa_file = "./data-raw/mali_kraken2.txt")
taxa_table10=format_wgs(taxa_file = "./data-raw/mali_phlan2.txt")

metadata_dir="./data-raw/metadata_mali.csv"
#format the raw taxonomic abudance table
taxa_table="./data-raw/mali_phlan2.txt"
taxa_tab1=format_wgs(taxa_file = taxa_table,method="metaphlan")
taxa_table="./data-raw/mali_kraken2.txt"
taxa_tab1=format_wgs(taxa_file = taxa_table,method="kraken")

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
plot1=tree_view(taxa_table =tab_s, metadata=metadata1,fdrs=fdrs1,test_metadata="group",fdr_cutoff=0.000001)

#PCoA plot
mds_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",taxa_level="ASV_or_strain",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")

#alpha diversity boxplot
alpha_plot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",palette_group=c("red","blue","orange","green"))

#taxa_barplot
taxa_barplot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",num_taxa=10,taxa_level="phylum",xlab_direction=1,legend_size=1)
taxa_barplot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",num_taxa=10,taxa_level="genus",xlab_direction=1,legend_size=1)

#taxa boxplot
taxa_boxplot(taxa_table = tab_s, metadata=metadata1,test_metadata="group",fdrs=fdrs1,log_norm=T,cutoff=0.01,palette_group=c("red","blue","orange","green"))


#pathway data
path_table="./data-raw/humann2_mali_unstratified.txt"
path_tab1=format_pathway(taxa_file = path_table,sep="\t")
colnames(path_tab1)=sapply(strsplit(colnames(path_tab1),"_S"),"[[",1)
write.table(path_tab1,file="./data-raw/humann2_pathway.txt",sep="\t",quote=F)

path_table="./data-raw/humann2_pathway.txt"
path_tab1=format_pathway(taxa_file = path_table,sep="\t")

metadata_dir="./data-raw/metadata_mali.csv"
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=1)

tab_s=table_subset(taxa_table = path_tab1,metadata=metadata1,stratify_by_metadata="time",stratify_by_value="1",taxa_file=F)

#PCoA plot
mds_plot(taxa_table = tab_s, metadata=metadata1,one_level=T,test_metadata="group",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")
mds_plot(taxa_table = tab_s, metadata=metadata1,one_level=T,log_norm=T,test_metadata="group",method_mds = "pcoa",palette_group=c("red","blue","orange","green"),distance_type="bray")

#alpha diversity boxplot
alpha_plot(taxa_table = tab_s, metadata=metadata1,one_level=T,test_metadata="group",palette_group=c("red","blue","orange","green"))

#taxa boxplot
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="group",method="wilcoxon")
taxa_boxplot(taxa_table = tab_s, metadata=metadata1,one_level=T,test_metadata="group",fdrs=fdrs1,log_norm=T,cutoff=0.01,palette_group=c("red","blue","orange","green"))




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

fdrs1=read.table(file="/Users/shansun/Google\ Drive/bartelt/new/Treatment_Day_wilcoxon2.csv",sep=",",row.names=1,header=T)
fdrs2=read.table(file="/Users/shansun/Google\ Drive/bartelt/new/Treatment_Day_wilcoxon.csv",sep=",",row.names=1,header=T)
p_compare(fdrs1,fdrs2,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=3,cor_method="spearman")
p_compare(fdrs1,fdrs2,p_col1=2,p_col2=2,indicator1=4,indicator2=4,point_color="black",lab_cutoff=3,cor_method="kendall")

metadata_dir="/Users/shansun/Google\ Drive/mc_set1/test/metadata_combined.csv"
taxa_dir="/Users/shansun/Google\ Drive/mc_set1/test/taxa_combined.csv"
taxa_tab1=format_asv(taxa_file = taxa_dir,biom=F,onefile = T,ASV=F,sep=",")
metadata1=meta_format(metadata=metadata_dir,metadata_sep=",",meta_sample_name_col=1)

tab_s=table_subset(taxa_table = taxa_tab1,metadata=metadata1,prevalence_cutoff=0.25, abundance_cutoff=0)
fdrs1=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="Case_Ctrl",method="lm",test_metadata_continuous=F,lm_anova=F,model_lm="factor(grade_school)+factor(ppi)+factor(batch)",random_effect_var=NULL)
fdrs2=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="Case_Ctrl",method="lme",test_metadata_continuous=F,lm_anova=F,model_lm="factor(grade_school)+factor(ppi)+factor(batch)",random_effect_var="external_subject_id")
fdrs3=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="Case_Ctrl",method="t.test",test_metadata_continuous=F,lm_anova=F,model_lm="factor(grade_school)+factor(ppi)+factor(batch)",random_effect_var="external_subject_id")

fdrs4=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="patient_age",method="lm",test_metadata_continuous=T,lm_anova=F,model_lm="factor(grade_school)+factor(ppi)+factor(batch)",random_effect_var="external_subject_id")
fdrs5=stat_test(taxa_table =tab_s,metadata=metadata1,test_metadata="patient_age",method="spearman",test_metadata_continuous=T,lm_anova=F,model_lm="factor(grade_school)+factor(ppi)+factor(batch)",random_effect_var="external_subject_id")



library(shiny)
runApp('/Users/shansun/git/plotmicrobiome')
#pathway bug, add taxa_file parameter in table subset
#kraken bug, not strain level in alpha diversity

library(BiocManager)
options(repos = BiocManager::repositories())
rsconnect::deployApp('/Users/shansun/git/plotmicrobiome')

#tree plot has three levels
#be able to change order of metadata in boxplots
#Variable preview: for data filter
