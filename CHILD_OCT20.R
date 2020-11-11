######################################################################################
######################################################################################
######################################################################################
# TITLE: Maternal consumption of artificially sweetened beverages during pregnancy is 
# associated with infant gut microbiota modifications and increased infant body mass index.
# CODE BY: ISABELLE LAFOREST-LAPOINTE

######################################################################################
######################################################################################
######################################################################################
# INSTALL AND LOAD LIBRARIES
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

to_install<-c("bestNormalize","biomformat","car","dada2","decontam","DECIPHER",
              "DESeq2","devtools","DirichletMultinomial","dplyr","exactRankTests",
              "FSA","gapminder","GenomicRanges","gganimate","ggplot2",
              "ggpubr","gifski","grid","gridExtra","Hmisc","httr","igraph",
              "leaflet","lme4","lmerTest","magrittr","MASS","Matrix","mediation",
              "microbiome","multcomp","nlme","openxlsx","patchwork","pbkrtest",
              "phangorn","pheatmap","phyloseq","picante","plotly","plyr","png",
              "randomForest","RColorBrewer","rcompanion","Rcpp","RcppArmadillo",
              "RCurl","reshape","reshape2","scales","slam","tibble","tidyr","vegan")

#BiocManager::install(to_install)

unlist(lapply(to_install, require, character.only = TRUE))

# Set seed for reproducible results
set.seed(1)

######################################################################################
######################################################################################
######################################################################################
# LOAD asv TABLE
asv<-read.csv("CHILDsequences.csv",row.names=1,header=T)
dim(asv) # 371 1268
colnames(asv)[1:5]
row.names(asv)

######################################################################################
######################################################################################
######################################################################################
# LOAD TAXONOMY TABLE (why is taxonomy table bigger -- to fix)
taxa<-read.csv("CHILDtaxa.csv",header=T,row.names=1)
dim(taxa) # 1702 7
colnames(taxa)
row.names(taxa)[1:5]

######################################################################################
######################################################################################
######################################################################################
# LOAD AND FORMAT METADATA IMPLEMENTED IN PREVIOUS MKDWN SCRIPT
CHILD<-read.csv("CHILD_implemented_MAR20.csv", row.names=1)
row.names(CHILD)
dim(CHILD) #200 38
CHILD$AS_bev4<-as.factor(CHILD$AS_bev4)
CHILD$sex<-as.factor(CHILD$sex)
CHILD$csec<-as.factor(CHILD$csec)
CHILD$CHILD_ID1<-as.factor(CHILD$CHILD_ID1)
CHILD$BF_3m_status<-as.factor(CHILD$BF_3m_status)
CHILD$FF_3m<-as.factor(CHILD$FF_3m)
CHILD$mom_sec<-as.factor(CHILD$mom_sec)
CHILD$infant_sec<-as.factor(CHILD$infant_sec)
CHILD$Dog_PN<-as.factor(CHILD$Dog_PN)
CHILD$Cat_PN<-as.factor(CHILD$Cat_PN)
CHILD$GDM_CHILD<-as.factor(CHILD$GDM_CHILD)
CHILD$solids_3m<-as.factor(CHILD$solids_3m)
CHILD$solids_6m<-as.factor(CHILD$solids_6m)
CHILD$Prebirth_abs_oralIV_yn<-as.factor(CHILD$Prebirth_abs_oralIV_yn)
CHILD$Mother_abs_birth_yn<-as.factor(CHILD$Mother_abs_birth_yn)
CHILD$Sample_time<-as.factor(CHILD$Sample_time)
summary(CHILD)
colnames(CHILD)

######################################################################################
######################################################################################
######################################################################################
# CONFIRM THAT ALL FILES HAVE THE ROW.NAMES AND ORDER
combi<-function(data1,data2) {
  combi<-intersect(row.names(data1),row.names(data2))
  return(combi)
}
both<-combi(asv,CHILD)  
asv<-asv[both,];dim(asv) #198 1268
CHILD<-CHILD[both,];dim(CHILD) #198  38
asv<-asv[,-c(which(apply(asv,2,sum)==0))];dim(asv) #198 954

both<-intersect(colnames(asv),row.names(taxa));length(both) #954
taxa<-taxa[both,]; dim(taxa) #954 7
taxa$joint<-as.factor(paste(taxa$Phylum,taxa$Class,taxa$Order,taxa$Family,taxa$Genus,taxa$Species,1:dim(taxa)[1],sep=""))
summary(taxa)

######################################################################################
######################################################################################
######################################################################################
# Obtain summary statistics from ASV dataset
# Total number of sequences included
sum(apply(asv,1,sum)) # 4553000
# Number of samples
dim(asv)[1] # 198
# Minimum, maximum, and mean number of sequences per sample
summary(apply(asv,1,sum))
# Number of ASVs
dim(asv)[2] # 954
# Minimum, maximum, and mean number of ASVs per sample
summary(apply(decostand(asv,method="pa"),1,sum))
# Minimum, maximum, and mean number of samples containing each ASVs
summary(apply(decostand(asv,method="pa"),2,sum))
# Who is the most ubiquitous ASV? Akkermansia municiphila
taxa[which(apply(decostand(asv,method="pa"),2,sum)==154),]$joint

######################################################################################
######################################################################################
######################################################################################
# TRANSFER TO PHYLOSEQ OBJECT
data_CHILD <- phyloseq(otu_table(asv, taxa_are_rows=FALSE), 
                        sample_data(CHILD), tax_table(as.matrix(taxa)))
data_CHILD

######################################################################################
######################################################################################
######################################################################################
# PHYLOSEQ VARIANCE STABILIZING TRANSFORMATION WITH BLIND TRANSFORMATION
test.phyloseq.dds<-phyloseq_to_deseq2(data_CHILD,~csec)

# Create function for computation (taken from online website of phyloseq)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds), 1, gm_mean))

# Blind version for all treatments
vst.blind <- DESeq2::varianceStabilizingTransformation(test.phyloseq.dds, blind=T)
vst.blind.Mat <- SummarizedExperiment::assay(vst.blind) # Extract transformed asv table
vst.blind.Mat<-t(vst.blind.Mat)
vst.blind.Mat[which(vst.blind.Mat<0)]<-0

# Prepare for ordination for all treatments
comm.vst.blind.Mat <- vegdist(vst.blind.Mat, "bray")
PCoA.comm.vst.blind.Mat<-capscale(vst.blind.Mat~1,distance="bray")
PCoA.comm.vst.blind.Mat$CA$eig[1:3]/sum(PCoA.comm.vst.blind.Mat$CA$eig)
PCoA.comm.vst.blind.Mat
eig <- PCoA.comm.vst.blind.Mat$CA$eig
eig[1]/sum(abs(eig)) #11.4%
eig[2]/sum(abs(eig)) #5.7%
eig[3]/sum(abs(eig)) #5.2%

######################################################################################
######################################################################################
######################################################################################
# TAXONOMICAL COMMUNITY ANALYSIS
comm.taxoGenera <- aggregate(t(otu_table(data_CHILD)), by=list(class=taxa$Genus), sum)
rownames(comm.taxoGenera) <- comm.taxoGenera[,1]
comm.taxoGenera <- comm.taxoGenera[,-1]
t(comm.taxoGenera)->comm.taxoGenera
colnames(comm.taxoGenera)
dim(comm.taxoGenera) #198 91

######################################################################################
######################################################################################
######################################################################################
# DIRICHLET MULTINOMIAL MIXTURES (DMM) MODELLING
count<-as.matrix(comm.taxoGenera)
cnts<-log10(colSums(count))
densityplot(cnts,xlim=range(cnts),xlab="Taxon")
fit<-mclapply(1:15,dmn,count=count,verbose=TRUE)
lplc<-sapply(fit,laplace)
plot(lplc, type="b", xlab="Number of Dirichlet Components",ylab="Model Fit")
(best <- fit[[which.min(lplc)]])
splom(log(fitted(best)))
p0 <- fitted(fit[[1]], scale=TRUE) # scale by theta
p3 <- fitted(best, scale=TRUE)
colnames(p3) <- paste("m", 1:4, sep="")
(meandiff <- colSums(abs(p3 - as.vector(p0))))
sum(meandiff)
diff <- rowSums(abs(p3 - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
df <- head(cbind(Mean=p0[o], p3[o,], diff=diff[o], cdiff), 10);df
heatmapdmn(count, fit[[1]], best, 30)

cluster<-numeric(length=dim(mixture(best))[1])
for (i in 1:dim(mixture(best))[1]){
  cluster[i]<-which(mixture(best)[i,]==max(mixture(best)[i,]))
}

######################################################################################
######################################################################################
######################################################################################
# DMM
md <- CHILD
md$CLUSTER_COLUMN<-cluster

#relevel cluster order so that it goes from young to mixed to older profile
summary(md$CLUSTER_COLUMN)
as.factor(md$CLUSTER_COLUMN)->md$CLUSTER_COLUMN
levels(md$CLUSTER_COLUMN)[levels(md$CLUSTER_COLUMN)=="1"]<-"2m"
levels(md$CLUSTER_COLUMN)[levels(md$CLUSTER_COLUMN)=="2"]<-"1"
levels(md$CLUSTER_COLUMN)[levels(md$CLUSTER_COLUMN)=="2m"]<-"2"
relevel(md$CLUSTER_COLUMN,"1")->md$CLUSTER_COLUMN
summary(md$CLUSTER_COLUMN)

genera <- t(comm.taxoGenera / rowSums(comm.taxoGenera))
apply(genera, 2, sum)
top15 <- head(names(rev(sort(rowSums(genera)))), 15)
top15[which(top15=="[Ruminococcus]")]<-"Ruminococcus"

which(colnames(md)=="CLUSTER_COLUMN")->t;t
order(md[colnames(genera), t])->order_cluster
md[order_cluster,]->md
anno<-as.data.frame(md$CLUSTER_COLUMN)
row.names(anno)<-row.names(md)
colnames(anno)<-"Bacterial genus"

colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=5, name="Spectral")), bias=3)(100)  #use to set colour (name =) and scale (bias =)

mat <- genera[top15, rownames(md)];dim(mat)
mat <- t(apply(mat, 1L, scales::rescale))   #uncomment to normalize each taxa (row)

######################################################################################
######################################################################################
######################################################################################
# FIGURE S1 ##########################################################################
######################################################################################
######################################################################################
######################################################################################
pheatmap::pheatmap(
  mat            = mat, 
  color          = colors,   #uncheck if setting the colour scale manual
  annotation_col = anno, 
  show_colnames  = FALSE,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  gaps_col       = cumsum(unname(table(md$CLUSTER_COLUMN))),
  labels_row     = sub("^.*__", "", top15)
)

both<-intersect(row.names(CHILD),row.names(md))
CHILD$cluster<-md[both,]$CLUSTER_COLUMN
CHILD$cluster<-as.factor(CHILD$cluster)
cbind(CHILD,scores(PCoA.comm.vst.blind.Mat)$sites)->CHILD

######################################################################################
######################################################################################
######################################################################################
# ALPHA-DIVERSITY AND SPECIES RICHNESS

alpha_diversity <- estimate_richness(otu_table(asv, taxa_are_rows=FALSE),
                                     measure = c("Shannon", "Observed","Chao1"))

#Chao1
CHILD$chao1<-alpha_diversity$Chao1
hist(log(CHILD$chao1)) #not normal
# Non parametric test
kruskal.test(chao1~cluster, data=CHILD)
dunnTest(data=CHILD,chao1~cluster, method="bh")

#Shannon
CHILD$alphadiv<-alpha_diversity$Shannon
hist(CHILD$alphadiv)
shapiro.test(CHILD$alphadiv) #not normal
# Non parametric test
kruskal.test(alphadiv~cluster, data=CHILD)
dunnTest(data=CHILD,alphadiv~cluster, method="bh")

# Evenness
H <- alpha_diversity$Shannon
S1 <- alpha_diversity$Observed
S <- log(S1)
evenness <- H/S
alpha_diversity$Evenness = evenness
CHILD$evenness<-alpha_diversity$Evenness
hist(CHILD$evenness)
shapiro.test(CHILD$evenness) #not normal
# Non parametric test
kruskal.test(evenness~cluster, data=CHILD)
dunnTest(data=CHILD,evenness~cluster, method="bh")

dev.off()
d1<-ggplot(CHILD, aes(cluster, alphadiv,fill=cluster, color=cluster))+theme_bw()+
  geom_boxplot(alpha=0.9)+
  xlab("")+ylab("ALPHA-DIVERSITY")+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"))+
  scale_fill_manual(name="Treatment",values=c("dodgerblue","forestgreen","red","deeppink"))+
  scale_color_manual(name="Treatment",values=c("black","black","black","black"))+
  annotate("text", x = 1, y = 4, label = "a", size=6) +
  annotate("text", x = 2, y = 4, label = "b", size=6) +
  annotate("text", x = 3, y = 4, label = "c", size=6) +
  annotate("text", x = 4, y = 4, label = "c", size=6) +
  guides(fill=F,color=F);d1

d2<-ggplot(CHILD, aes(cluster, chao1,fill=cluster, color=cluster))+theme_bw()+
  geom_boxplot(alpha=0.9)+
  xlab("")+ylab("SPECIES RICHNESS")+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"))+
  scale_fill_manual(name="Treatment",values=c("dodgerblue","forestgreen","red","deeppink"))+
  scale_color_manual(name="Treatment",values=c("black","black","black","black"))+
  annotate("text", x = 1, y = 98, label = "a", size=6) +
  annotate("text", x = 2, y = 98, label = "b", size=6) +
  annotate("text", x = 3, y = 98, label = "c", size=6) +
  annotate("text", x = 4, y = 98, label = "d", size=6) +
  guides(fill=F,color=F);d2

d3<-ggplot(CHILD, aes(cluster, evenness,fill=cluster, color=cluster))+theme_bw()+
  geom_boxplot(alpha=0.9)+
  xlab("CLUSTER")+ylab("SPECIES EVENNESS")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(name="Treatment",values=c("dodgerblue","forestgreen","red","deeppink"))+
  scale_color_manual(name="Treatment",values=c("black","black","black","black"))+
  annotate("text", x = 1, y = 1, label = "a", size=6) +
  annotate("text", x = 2, y = 1, label = "b", size=6) +
  annotate("text", x = 3, y = 1, label = "c", size=6) +
  annotate("text", x = 4, y = 1, label = "c", size=6) +
  guides(fill=F,color=F);d3

######################################################################################
######################################################################################
######################################################################################
# CREATE A SEPARATE FILE FOR EACH CLUSTER
Cluster1<-CHILD[which(CHILD$cluster=="1"),];dim(Cluster1) #48 44
Cluster2<-CHILD[which(CHILD$cluster=="2"),];dim(Cluster2) #59 44
Cluster3<-CHILD[which(CHILD$cluster=="3"),];dim(Cluster3) #47 44
Cluster4<-CHILD[which(CHILD$cluster=="4"),];dim(Cluster4) #44 44

table(summary(as.factor(Cluster1$CHILD_ID1)))
table(summary(as.factor(Cluster2$CHILD_ID1)))
table(summary(as.factor(Cluster3$CHILD_ID1)))
table(summary(as.factor(Cluster4$CHILD_ID1)))

######################################################################################
######################################################################################
######################################################################################
# ORDINATIONS PER CLUSTER

# CLUSTER 1
data_clust1 <- phyloseq(otu_table(asv, taxa_are_rows=FALSE), 
                        sample_data(Cluster1), tax_table(as.matrix(taxa)))

# PHYLOSEQ VARIANCE STABILIZING TRANSFORMATION WITH BLIND TRANSFORMATION
test.phyloseq.dds<-phyloseq_to_deseq2(data_clust1,~csec)

# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds), 1, gm_mean))

# Blind version for all treatments
vst.blind <- DESeq2::varianceStabilizingTransformation(test.phyloseq.dds, blind=T)
vst.blind.Mat1 <- SummarizedExperiment::assay(vst.blind) # Extract transformed asv table
vst.blind.Mat1<-t(vst.blind.Mat1)
vst.blind.Mat1[which(vst.blind.Mat1<0)]<-0

# Prepare for ordination for all treatments
comm.vst.blind.Mat1 <- vegdist(vst.blind.Mat1, "bray")
PCoA.comm.vst.blind.Mat1<-capscale(vst.blind.Mat1~1,distance="bray")
PCoA.comm.vst.blind.Mat1$CA$eig[1:3]/sum(PCoA.comm.vst.blind.Mat1$CA$eig)
eig1 <- PCoA.comm.vst.blind.Mat1$CA$eig
eig[1]/sum(abs(eig)) #14.7%
eig[2]/sum(abs(eig)) #8.7%
eig[3]/sum(abs(eig)) #7.2%
PCoA.comm.vst.blind.Mat1->PCoA_clust1

# CLUSTER 2
data_clust2 <- phyloseq(otu_table(asv, taxa_are_rows=FALSE), 
                        sample_data(Cluster2), tax_table(as.matrix(taxa)))

# PHYLOSEQ VARIANCE STABILIZING TRANSFORMATION WITH BLIND TRANSFORMATION
test.phyloseq.dds<-phyloseq_to_deseq2(data_clust2,~csec)

# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds), 1, gm_mean))

# Blind version for all treatments
vst.blind <- DESeq2::varianceStabilizingTransformation(test.phyloseq.dds, blind=T)
vst.blind.Mat2 <- SummarizedExperiment::assay(vst.blind) # Extract transformed asv table
vst.blind.Mat2<-t(vst.blind.Mat2)
vst.blind.Mat2[which(vst.blind.Mat2<0)]<-0

# Prepare for ordination for all treatments
comm.vst.blind.Mat2 <- vegdist(vst.blind.Mat2, "bray")
PCoA.comm.vst.blind.Mat2<-capscale(vst.blind.Mat2~1,distance="bray")
PCoA.comm.vst.blind.Mat2$CA$eig[1:3]/sum(PCoA.comm.vst.blind.Mat2$CA$eig)
eig <- PCoA.comm.vst.blind.Mat2$CA$eig
eig1[1]/sum(abs(eig1)) #8.8%
eig1[2]/sum(abs(eig1)) #7.7%
eig1[3]/sum(abs(eig1)) #6.4%
PCoA.comm.vst.blind.Mat2->PCoA_clust2

# CLUSTER 3
data_clust3 <- phyloseq(otu_table(asv, taxa_are_rows=FALSE), 
                        sample_data(Cluster3), tax_table(as.matrix(taxa)))

# PHYLOSEQ VARIANCE STABILIZING TRANSFORMATION WITH BLIND TRANSFORMATION
test.phyloseq.dds<-phyloseq_to_deseq2(data_clust3,~csec)

# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds,
                                        geoMeans = apply(counts(test.phyloseq.dds), 1, gm_mean))

# Blind version for all treatments
vst.blind <- DESeq2::varianceStabilizingTransformation(test.phyloseq.dds, blind=T)
vst.blind.Mat3 <- SummarizedExperiment::assay(vst.blind) # Extract transformed asv table
vst.blind.Mat3<-t(vst.blind.Mat3)
vst.blind.Mat3[which(vst.blind.Mat3<0)]<-0

# Prepare for ordination for all treatments
comm.vst.blind.Mat3 <- vegdist(vst.blind.Mat3, "bray")
PCoA.comm.vst.blind.Mat3<-capscale(vst.blind.Mat3~1,distance="bray")
PCoA.comm.vst.blind.Mat3$CA$eig[1:3]/sum(PCoA.comm.vst.blind.Mat3$CA$eig)
eig <- PCoA.comm.vst.blind.Mat3$CA$eig
eig[1]/sum(abs(eig)) #14.1%
eig[2]/sum(abs(eig)) #8.7%
eig[3]/sum(abs(eig)) #6.8%
PCoA.comm.vst.blind.Mat3->PCoA_clust3

# CLUSTER 4
data_clust4 <- phyloseq(otu_table(asv, taxa_are_rows=FALSE), 
                        sample_data(Cluster4), tax_table(as.matrix(taxa)))

# PHYLOSEQ VARIANCE STABILIZING TRANSFORMATION WITH BLIND TRANSFORMATION
test.phyloseq.dds<-phyloseq_to_deseq2(data_clust4,~csec)

# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds,
                                        geoMeans = apply(counts(test.phyloseq.dds), 1, gm_mean))

# Blind version for all treatments
vst.blind <- DESeq2::varianceStabilizingTransformation(test.phyloseq.dds, blind=T)
vst.blind.Mat4 <- SummarizedExperiment::assay(vst.blind) # Extract transformed asv table
vst.blind.Mat4<-t(vst.blind.Mat4)
vst.blind.Mat4[which(vst.blind.Mat4<0)]<-0

# Prepare for ordination for all treatments
comm.vst.blind.Mat4 <- vegdist(vst.blind.Mat4, "bray")
PCoA.comm.vst.blind.Mat4<-capscale(vst.blind.Mat4~1,distance="bray")
PCoA.comm.vst.blind.Mat4$CA$eig[1:3]/sum(PCoA.comm.vst.blind.Mat4$CA$eig)
eig <- PCoA.comm.vst.blind.Mat4$CA$eig
eig[1]/sum(abs(eig)) #9.7%
eig[2]/sum(abs(eig)) #8.6%
eig[3]/sum(abs(eig)) #7.3%
PCoA.comm.vst.blind.Mat4->PCoA_clust4

asv1<-as.matrix(otu_table(data_clust1, taxa_are_rows = F));dim(asv1)
asv2<-as.matrix(otu_table(data_clust2, taxa_are_rows = F));dim(asv2)
asv3<-as.matrix(otu_table(data_clust3, taxa_are_rows = F));dim(asv3)
asv4<-as.matrix(otu_table(data_clust4, taxa_are_rows = F));dim(asv4)
colnames(asv1)->asv_sequences
colnames(asv1)<-paste("asv",c(1:length(colnames(asv1))),sep="")
colnames(asv2)<-paste("asv",c(1:length(colnames(asv2))),sep="")
colnames(asv3)<-paste("asv",c(1:length(colnames(asv3))),sep="")
colnames(asv4)<-paste("asv",c(1:length(colnames(asv4))),sep="")

x1<-vegan::diversity(otu_table(asv1, taxa_are_rows=FALSE), index="shannon",
                     MARGIN=1, base=exp(1))
Cluster1$alphadiv<-x1
hist(Cluster1$alphadiv)
# Non parametric
kruskal.test(alphadiv~AS_bev4, data=Cluster1) #NS
kruskal.test(alphadiv~csec, data=Cluster1) #NS
kruskal.test(alphadiv~BF_3m_status, data=Cluster1) #NS
kruskal.test(alphadiv~diet_3m, data=Cluster1) #NS

# Compute Chao1 richness index with phyloseq
x1<-estimate_richness(otu_table(asv1, taxa_are_rows=FALSE),measures="Chao1")
Cluster1$chao1<-x1[,1]
hist(log(Cluster1$chao1))
kruskal.test(chao1~AS_bev4, data=Cluster1) #NS
kruskal.test(chao1~csec, data=Cluster1) #NS
kruskal.test(chao1~BF_3m_status, data=Cluster1) #NS
kruskal.test(chao1~diet_3m, data=Cluster1) #NS

######################################################################################
######################################################################################
######################################################################################
# FIGURE 1K ##########################################################################
######################################################################################
######################################################################################
######################################################################################

test <- qplot(CHILD$MDS1, CHILD$MDS2, xlab="PCoA1",
               ylab="PCoA2", shape=cluster, fill=cluster, color=cluster, data=(CHILD))
t1<-(test +
    scale_shape_manual(name="CLUSTER",values=c(21,21,21,21),
                       labels=c("1 : 3 months only", "2 : 3 months and 12 months", "3 : 3 months and 12 months", "4 : mostly 12 months"))+
    stat_ellipse(level=0.95, geom = "polygon", alpha = 1/6, linetype=2, aes(fill = cluster,color=cluster))+
    theme_bw() + theme(legend.title = element_text(colour="black", size=18, face="bold"),
                       legend.text = element_text(colour="black", size = 18),
                       axis.title=element_text(face="bold",size="18", color="black"),
                       legend.background = element_blank(),
                       plot.title=element_text(face="bold",size=18))+
    scale_color_manual(name="CLUSTER",values=c("dodgerblue","forestgreen","red","deeppink","black","black","black","black"),
                       labels=c("1 : 3 months only", "2 : 3 months and 12 months", "3 : 3 months and 12 months", "4 : mostly 12 months"),
                       guide=FALSE)+
    scale_fill_manual(name="CLUSTER",values=c("dodgerblue","forestgreen","red","deeppink"),
                      labels=c("1 : 3 months only", "2 : 3 months and 12 months", "3 : 3 months and 12 months", "4 : mostly 12 months"))+
    geom_point(aes(size=0.5,color=site),show.legend=FALSE)+
    guides(size=FALSE, shape = guide_legend(override.aes = list(size=4)))+
    #guides(size=FALSE, shape=FALSE, color=FALSE, fill=FALSE)+
  annotate("text", x = -1.1, y = 2, label = "p = 0.001", size=8)+
  annotate("text", x = -1.1, y = 2.2, label = "R2 = 0.04", size=8))

######################################################################################
######################################################################################
######################################################################################
# Test of statistical differences in all variables between clusters
kruskal.test(AS_bev4~cluster, data=CHILD) #NS
kruskal.test(AS_soda4~cluster, data=CHILD) #NS
kruskal.test(SS_bev4~cluster, data=CHILD) #NS
kruskal.test(SS_soda4~cluster, data=CHILD) #NS
kruskal.test(sex~cluster, data=CHILD) #NS
kruskal.test(mom_bmi_best~cluster, data=CHILD) #NS
kruskal.test(race_mom4~cluster, data=CHILD) #NS
kruskal.test(site~cluster, data=CHILD) #NS
kruskal.test(older_sibs3~cluster, data=CHILD) #NS
kruskal.test(mom_age~cluster, data=CHILD) #NS
kruskal.test(ses_mom_edu~cluster, data=CHILD) #NS
kruskal.test(GDM_CHILD~cluster, data=CHILD) #NS
kruskal.test(age_stool_3m_mos~cluster, data=CHILD) #NS
kruskal.test(age_stool_1y_mos~cluster, data=CHILD) #NS
kruskal.test(solids_3m~cluster, data=CHILD) #NS
kruskal.test(solids_6m~cluster, data=CHILD) #NS
kruskal.test(calories~cluster, data=CHILD) #NS
kruskal.test(addsugar~cluster, data=CHILD) #NS
kruskal.test(totsugar~cluster, data=CHILD) #NS
kruskal.test(Prebirth_abs_oralIV_yn~cluster, data=CHILD) #NS
kruskal.test(Child_6moto1Y_abs_total~cluster, data=CHILD) #NS
kruskal.test(Dog_PN~cluster, data=CHILD) #NS
kruskal.test(Cat_PN~cluster, data=CHILD) #NS
kruskal.test(mom_sec~cluster, data=CHILD) #NS
kruskal.test(infant_sec~cluster, data=CHILD) #NS
kruskal.test(diet_3m~cluster, data=CHILD) #NS

contrasts(CHILD$cluster) <- contr.treatment
# SAMPLE TIME
par(mfrow=c(1,1))
plot(Sample_time~cluster, data=CHILD)
kruskal.test(Sample_time~cluster, data=CHILD) # < 2.2e-16
model <- glm(Sample_time ~ cluster, data=CHILD, family=binomial(logit))
summary(model)
Anova(model)

tb<-table(CHILD$Sample_time,CHILD$cluster)
lbs<-matrix(NA, 8,3)
lbs[,1]<-rep(levels(CHILD$cluster),2)
lbs[,2]<-c(rep(levels(CHILD$Sample_time)[1],4),rep(levels(CHILD$Sample_time)[2],4))
lbs[,3]<-c(tb[1,],tb[2,])
colnames(lbs)<-c("cluster","Sample_time","freq")
lbs<-as.data.frame(lbs)
lbs$freq<-as.numeric(as.character(lbs$freq))
lbs$label_ypos<-as.numeric(as.character(c(sum(tb[,1]),sum(tb[,2]),sum(tb[,3]),
                                          sum(tb[,4])+1.5,"NA",tb[2,2],tb[2,3],tb[2,4])))
lbs$Sample_time<-as.factor(lbs$Sample_time)
lbs$Sample_time<-factor(lbs$Sample_time,levels=c(levels(lbs$Sample_time)[2:1]))
summary(lbs)

p1<-ggplot(lbs, aes(cluster,freq,fill=Sample_time))+theme_bw()+
  geom_bar(stat="identity",color="black",width = 0.7)+
  xlab("")+ylab("COUNT")+
  geom_text(aes(y=label_ypos, label=freq), vjust=1.6, 
            color="black", size=5)+
  theme(axis.text.x  = element_text(size=20, color="black"),
        axis.text.y  = element_text(size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position=c(0.82,0.89),
        legend.background = element_blank())+
  annotate("text", x = 1.2, y = 62, label = "p < 2.2e-16", size=6)+
  scale_fill_brewer(name="INFANT AGE",palette="Blues",labels=c("3 months","12 months"))+
  guides(color=F);p1

# BIRTH MODE
plot(csec~cluster, data=CHILD)
kruskal.test(csec~cluster, data=CHILD) # 0.0515
model <- glm(csec ~ cluster, data=CHILD, family=binomial(logit))
summary(model)
Anova(model)
summary(glht(model, mcp(cluster="Tukey")))

tb<-table(CHILD$csec,CHILD$cluster)
lbs<-matrix(NA, length(levels(CHILD$csec))*length(levels(CHILD$cluster)),3)
lbs[,1]<-rep(levels(CHILD$cluster),2)
lbs[,2]<-c(rep(levels(CHILD$csec)[2],4),rep(levels(CHILD$csec)[1],4))
lbs[,3]<-c(tb[2,],tb[1,])
colnames(lbs)<-c("cluster","csec","freq")
lbs<-as.data.frame(lbs)
lbs$freq<-as.numeric(as.character(lbs$freq))
lbs$label_ypos<-as.numeric(as.character(c(tb[2,1],tb[2,2],tb[2,3],tb[2,4],
                                          sum(tb[,1]),sum(tb[,2]),sum(tb[,3]),sum(tb[,4]))))
summary(lbs)
lbs$csec1<-NA
lbs$csec1[which(lbs$csec==0)]<-"Vaginal"
lbs$csec1[which(lbs$csec==1)]<-"Csec"
lbs$csec1<-as.factor(lbs$csec1)
lbs$csec1<-factor(lbs$csec1,levels=c(levels(lbs$csec1)[2:1]))

p2<-ggplot(lbs, aes(cluster,freq,fill=csec1))+theme_bw()+
  geom_bar(stat="identity",color="black",width = 0.7)+
  xlab("")+ylab("")+
  geom_text(aes(y=label_ypos, label=freq), vjust=1.6, 
            color="black", size=5)+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position=c(0.82,0.89),
        legend.background = element_blank())+
  scale_fill_brewer(name="BIRTH MODE",palette="Blues",labels=c("Vaginal","Cesarean"))+
  annotate("text", x = 1, y = 62, label = "p = 0.05", size=6)+
  guides(color=F);p2

# AS_BEV4
plot(AS_bev4~cluster, data=CHILD)
kruskal.test(AS_bev4~cluster, data=CHILD) # 0.46
model <- glm(AS_bev4 ~ cluster, data=CHILD, family=binomial(logit))
summary(model)
Anova(model)
summary(glht(model, mcp(cluster="Tukey")))

tb<-table(CHILD$AS_bev4,CHILD$cluster)
lbs<-matrix(NA, length(levels(CHILD$AS_bev4))*length(levels(CHILD$cluster)),3)
lbs[,1]<-rep(levels(CHILD$cluster),2)
lbs[,2]<-c(rep(levels(CHILD$AS_bev4)[1],4),rep(levels(CHILD$AS_bev4)[2],4))
lbs[,3]<-c(tb[1,],tb[2,])
colnames(lbs)<-c("cluster","AS_bev4","freq")
lbs<-as.data.frame(lbs)
lbs$freq<-as.numeric(as.character(lbs$freq))
lbs$label_ypos<-as.numeric(as.character(c(sum(tb[,1]),sum(tb[,2]),sum(tb[,3]),
                                          sum(tb[,4]),tb[2,1],tb[2,2],tb[2,3],tb[2,4])))
summary(lbs)
lbs$AS_bev41[which(lbs$AS_bev4==0)]<-"Low"
lbs$AS_bev41[which(lbs$AS_bev4==3)]<-"High"
lbs$AS_bev41<-as.factor(lbs$AS_bev41)
lbs$AS_bev41<-factor(lbs$AS_bev41,levels=levels(lbs$AS_bev41)[2:1])

p3<-ggplot(lbs, aes(cluster,freq,fill=AS_bev41))+theme_bw()+
  geom_bar(stat="identity",color="black",width = 0.7)+
  xlab("")+ylab("COUNT")+
  geom_text(aes(y=label_ypos, label=freq), vjust=1.6, 
            color="black", size=5)+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position=c(0.82,0.89),
        legend.background = element_blank())+
  annotate("text", x = 1, y = 62, label = "p = 0.46", size=6)+
  scale_fill_brewer(name="ASB",palette="Blues",labels=c("Zero","High"))+
  guides(color=F);p3

# BF_3m
plot(BF_3m_status~cluster, data=CHILD)
kruskal.test(BF_3m_status~cluster, data=CHILD) # 0.0002
model <- glm(BF_3m_status ~ cluster, data=CHILD, family=binomial(logit))
summary(model)
Anova(model)
summary(glht(model, mcp(cluster="Tukey")))

tb<-table(CHILD$BF_3m_status,CHILD$cluster)
lbs<-matrix(NA, length(levels(CHILD$BF_3m_status))*length(levels(CHILD$cluster)),3)
lbs[,1]<-rep(levels(CHILD$cluster),3)
lbs[,2]<-c(rep(levels(CHILD$BF_3m_status)[1],4),rep(levels(CHILD$BF_3m_status)[2],4),
           rep(levels(CHILD$BF_3m_status)[3],4))
lbs[,3]<-c(tb[1,],tb[2,],tb[3,])
colnames(lbs)<-c("cluster","BF_3m_status","freq")
lbs<-as.data.frame(lbs)
lbs$label_ypos<-as.numeric(as.character(c(sum(tb[,1]),sum(tb[,2]),sum(tb[,3]),
                                          sum(tb[,4]),sum(tb[3:2,1]),sum(tb[3:2,2]),
                                          sum(tb[3:2,3]),sum(tb[3:2,4]),tb[3,1],
                                          tb[3,2],tb[3,3],tb[3,4])))
lbs$freq<-as.numeric(as.character(lbs$freq))
summary(lbs)

p4<-ggplot(lbs, aes(cluster,freq,fill=BF_3m_status))+theme_bw()+
  geom_bar(stat="identity",color="black",width = 0.7)+
  xlab("")+ylab("")+
  geom_text(aes(y=label_ypos, label=freq), vjust=1.6, 
            color="black", size=5)+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position=c(0.82,0.87),
        legend.background = element_blank())+
  annotate("text", x = 1, y = 62, label = "p < 0.001", size=6)+
  scale_fill_brewer(name="BF at 3M",palette="Blues", labels=c("Exclusive","Partial","None"))+
  guides(color=F);p4

# FF_3m
plot(FF_3m~cluster, data=CHILD)
kruskal.test(FF_3m~cluster, data=CHILD) # 0.0002
model <- glm(FF_3m ~ cluster, data=CHILD, family=binomial(logit))
summary(model)
Anova(model)
summary(glht(model, mcp(cluster="Tukey")))

tb<-table(CHILD$FF_3m,CHILD$cluster)
lbs<-matrix(NA, length(levels(CHILD$FF_3m))*length(levels(CHILD$cluster)),3)
lbs[,1]<-rep(levels(CHILD$cluster),2)
lbs[,2]<-c(rep(levels(CHILD$FF_3m)[1],4),rep(levels(CHILD$FF_3m)[2],4))
lbs[,3]<-c(tb[1,],tb[2,])
colnames(lbs)<-c("cluster","FF_3m","freq")
lbs<-as.data.frame(lbs)
lbs$freq<-as.numeric(as.character(lbs$freq))
lbs$label_ypos<-as.numeric(as.character(c(sum(tb[,1]),sum(tb[,2]),sum(tb[,3]),sum(tb[,4]),tb[2,1],tb[2,2],tb[2,3],tb[2,4])))
summary(lbs)
lbs$FF_3m1<-NA
lbs$FF_3m1[which(lbs$FF_3m==0)]<-"No"
lbs$FF_3m1[which(lbs$FF_3m==1)]<-"Yes"

p5<-ggplot(lbs, aes(cluster,freq,fill=FF_3m1))+theme_bw()+
  geom_bar(stat="identity",color="black",width = 0.7)+
  xlab("CLUSTER")+ylab("COUNT")+
  geom_text(aes(y=label_ypos, label=freq), vjust=1.6, 
            color="black", size=5)+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position=c(0.82,0.89),
        legend.background = element_blank())+
  annotate("text", x = 1, y = 62, label = "p < 0.001", size=6)+
  scale_fill_brewer(name="FF at 3M",palette="Blues")+
  guides(color=F);p5

# Mother_abs_birth_yn
plot(Mother_abs_birth_yn~cluster, data=CHILD)
kruskal.test(Mother_abs_birth_yn~cluster, data=CHILD) # 0.006
model <- glm(Mother_abs_birth_yn ~ cluster, data=CHILD, family=binomial(logit))
summary(model)
Anova(model)
summary(glht(model, mcp(cluster="Tukey")))

tb<-table(CHILD$Mother_abs_birth_yn,CHILD$cluster)
lbs<-matrix(NA, length(levels(CHILD$Mother_abs_birth_yn))*length(levels(CHILD$cluster)),3)
lbs[,1]<-rep(levels(CHILD$cluster),2)
lbs[,2]<-c(rep(levels(CHILD$Mother_abs_birth_yn)[1],4),rep(levels(CHILD$Mother_abs_birth_yn)[2],4))
lbs[,3]<-c(tb[1,],tb[2,])
colnames(lbs)<-c("cluster","Mother_abs_birth_yn","freq")
lbs<-as.data.frame(lbs)
lbs$freq<-as.numeric(as.character(lbs$freq))
lbs$label_ypos<-as.numeric(as.character(c(sum(tb[,1]),sum(tb[,2]),sum(tb[,3]),sum(tb[,4]),tb[2,1],tb[2,2],tb[2,3],tb[2,4])))
summary(lbs)
lbs$ABX<-NA
lbs$ABX[which(lbs$Mother_abs_birth_yn==0)]<-"No"
lbs$ABX[which(lbs$Mother_abs_birth_yn==1)]<-"Yes"

p6<-ggplot(lbs, aes(cluster,freq,fill=ABX))+theme_bw()+
  geom_bar(stat="identity",color="black",width = 0.7)+
  xlab("CLUSTER")+ylab("")+
  geom_text(aes(y=label_ypos, label=freq), vjust=1.6, 
            color="black", size=5)+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 16),
        legend.position=c(0.82,0.89),
        legend.background = element_blank())+
  annotate("text", x = 1, y = 62, label = "p = 0.005", size=6)+
  scale_fill_brewer(name="ABX",palette="Blues")+
  guides(color=F);p6

# BF duration
plot(BF_duration_imp~cluster, data=CHILD)
kruskal.test(BF_duration_imp~cluster, data=CHILD) # <0.0001
dunnTest(data=CHILD,BF_duration_imp~cluster, method="bh")

p7<-ggplot(CHILD, aes(cluster, BF_duration_imp))+theme_bw()+
  geom_boxplot(aes(fill="main"))+
  xlab("CLUSTER")+ylab("BF DURATION (mths)")+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"))+
    scale_fill_brewer(palette="Blues")+
  annotate("text", x = 1, y = 26, label = "b", size=5, col="blue") +
  annotate("text", x = 2, y = 26, label = "a", size=5, col="blue") +
  annotate("text", x = 3, y = 26, label = "a", size=5, col="blue") +
  annotate("text", x = 4, y = 26, label = "a", size=5, col="blue") +
  annotate("text", x = 1.2, y = 28, label = "p < 0.0001", size=6)+
  guides(fill=F);p7

# Hei2010

plot(hei2010~cluster, data=CHILD)
kruskal.test(hei2010~cluster, data=CHILD) # 0.01
dunnTest(data=CHILD,hei2010~cluster, method="bh")

p8<-ggplot(CHILD, aes(cluster, hei2010))+theme_bw()+
  geom_boxplot(aes(fill="main"))+
  xlab("CLUSTER")+ylab("HEI2010")+
  theme(axis.text.x  = element_text (size=20, color="black"),
        axis.text.y  = element_text (size=20, color="black"),
        axis.title.x  = element_text(size=24, color="black"),
        axis.title.y  = element_text(size=24, color="black"))+
  scale_fill_brewer(palette="Blues")+
  annotate("text", x = 1, y = 95, label = "a", size=5, col="blue") +
  annotate("text", x = 2, y = 95, label = "ab", size=5, col="blue") +
  annotate("text", x = 3, y = 95, label = "b", size=5, col="blue") +
  annotate("text", x = 4, y = 95, label = "ab", size=5, col="blue") +
  annotate("text", x = 1, y = 99, label = "p = 0.01", size=6)+
  guides(fill=F);p8

######################################################################################
######################################################################################
######################################################################################
# FIGURE 1 ###########################################################################
######################################################################################
######################################################################################
######################################################################################

# Figure to be edited externally
( d1 | d2 | d1 | d2) / ( p1 | p2 |  p1 | p2 ) / ( p3 | p4 | p3 | p4 ) / (  p5 | p6 | p7 | p8 )

# Test for differences in age at stool sample

kruskal.test(age_stool_3m_mos~cluster, data=CHILD) # 0.46
dunnTest(data=CHILD,age_stool_3m_mos~cluster, method="bh")
p9<-ggplot(CHILD, aes(cluster, age_stool_3m_mos))+theme_bw()+
  geom_boxplot(aes(fill="main"))+
  xlab("CLUSTER")+ylab("Age Sample 3M")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_brewer(palette="Blues")+
  guides(fill=F);p9

kruskal.test(age_stool_1y_mos~cluster, data=CHILD) # 0.53
p10<-ggplot(CHILD, aes(cluster, age_stool_1y_mos))+theme_bw()+
  geom_boxplot(aes(fill="main"))+
  xlab("CLUSTER")+ylab("Age Sample 1Y")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_brewer(palette="Blues")+
  guides(fill=F);p10

aggregate(cbind(CHILD$age_stool_3m_mos,CHILD$age_stool_1y_mos),by=list(CHILD$cluster),mean)
aggregate(cbind(CHILD$age_stool_3m_mos,CHILD$age_stool_1y_mos),by=list(CHILD$cluster),sd)

ggarrange(p9,p10)

######################################################################################
######################################################################################
######################################################################################
# BARCHART
rel_abun_genus<-sweep(t(comm.taxoGenera), 2, apply(t(comm.taxoGenera),2,sum), `/`)
apply(rel_abun_genus,2,sum) # Here you confirm that all samples have a total relative abundance of 1
t<-length(apply(rel_abun_genus,1,sum))
order(apply(rel_abun_genus,1,sum))[(t-14):t]
top20_genus<-names(apply(rel_abun_genus,1,sum)[order(apply(rel_abun_genus,1,sum))[(t-14):t]]);top20_genus
rel_abun_genus[top20_genus,]->rel_abun_genus
rel_abun_genus<-t(rel_abun_genus)
dim(CHILD);dim(rel_abun_genus)
cbind(CHILD,rel_abun_genus)->metadataCHILD
summary(metadataCHILD)

# Create aggregated dataset for barchart figure
data_long <- gather(metadataCHILD, class, relative_abundance, "Prevotella":"Bacteroides", factor_key=TRUE)
summary(data_long)
agg <- aggregate(data_long$relative_abundance, by=list(class=data_long$class,cluster=data_long$cluster), mean)
1-aggregate(agg$x,by=list(agg$cluster),sum)$x
others<-matrix(NA,4,3);others
others[,1]<-"other"
others[,2]<-1:4;others
others[,3]<-1-aggregate(agg$x,by=list(agg$cluster),sum)$x
others<-as.data.frame(others)
colnames(others)<-colnames(agg)
agg<-rbind(agg,others)
agg$x<-as.numeric(agg$x)
agg$class<-factor(agg$class, levels=levels(agg$class)[c(16,1:15)])

# Test for statistical differences in top 10 genera between clusters

comm.taxoGenera <- aggregate(t(otu_table(data_CHILD)), by=list(class=taxa$Genus), sum)
rownames(comm.taxoGenera) <- comm.taxoGenera[,1]
comm.taxoGenera <- comm.taxoGenera[,-1]
t(comm.taxoGenera)->comm.taxoGenera
colnames(comm.taxoGenera)
dim(comm.taxoGenera) #198 91

rel_abun_genus<-sweep(t(comm.taxoGenera), 2, apply(t(comm.taxoGenera),2,sum), `/`)
apply(rel_abun_genus,2,sum) # Here you confirm that all samples have a total relative abundance of 1
t<-length(apply(rel_abun_genus,1,sum))
order(apply(rel_abun_genus,1,sum))[(t-14):t]
top15_genus<-names(apply(rel_abun_genus,1,sum)[order(apply(rel_abun_genus,1,sum))[(t-14):t]]);top15_genus
rel_abun_genus[top15_genus,]->rel_abun_genus
rel_abun_genus<-t(rel_abun_genus)
dim(CHILD);dim(rel_abun_genus)
cbind(CHILD,rel_abun_genus)->metadataCHILD
summary(metadataCHILD)

kruskal.test(Bacteroides~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Bacteroides~cluster, method="bh")
plot(Bacteroides~cluster, data=metadataCHILD)
s1<-ggplot(metadataCHILD, aes(cluster, Bacteroides*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  xlab("")+ylab("RELATIVE ABUNDANCE")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values="#006837")+
  ggtitle("Bacteroides")+
  annotate("text", x = 1, y = 100, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 100, label = "b", size=5, col="black") +
  annotate("text", x = 3, y = 100, label = "a", size=5, col="black") +
  annotate("text", x = 4, y = 100, label = "b", size=5, col="black") +
  guides(fill=F);s1

kruskal.test(Escherichia~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Escherichia~cluster, method="bh")
s2<-ggplot(metadataCHILD, aes(cluster, Escherichia*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  xlab("")+ylab("")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#1a9850"))+
  ggtitle("Escherichia")+
  annotate("text", x = 1, y = 100, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 100, label = "b", size=5, col="black") +
  annotate("text", x = 3, y = 100, label = "a", size=5, col="black") +
  annotate("text", x = 4, y = 100, label = "b", size=5, col="black") +
  guides(fill=F);s2

kruskal.test(Bifidobacterium~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Bifidobacterium~cluster, method="bh")
s3<-ggplot(metadataCHILD, aes(cluster, Bifidobacterium*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  xlab("CLUSTER")+ylab("")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#66bd63"))+
  ggtitle("Bifidobacterium")+
  annotate("text", x = 1, y = 63, label = "ab", size=5, col="black") +
  annotate("text", x = 2, y = 63, label = "a", size=5, col="black") +
  annotate("text", x = 3, y = 63, label = "b", size=5, col="black") +
  annotate("text", x = 4, y = 63, label = "a", size=5, col="black") +
  ylim(c(0,63))+
  guides(fill=F);s3

kruskal.test(`[Ruminococcus]`~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,`[Ruminococcus]`~cluster, method="bh")
s4<-ggplot(metadataCHILD, aes(cluster, `[Ruminococcus]`*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  xlab("CLUSTER")+ylab("")+
  ggtitle("Ruminococcus")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#a6d96a"))+
  annotate("text", x = 1, y = 82, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 82, label = "b", size=5, col="black") +
  annotate("text", x = 3, y = 82, label = "c", size=5, col="black") +
  annotate("text", x = 4, y = 82, label = "b", size=5, col="black") +
  ylim(c(0,82))+
  guides(fill=F);s4

kruskal.test(Akkermansia~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Akkermansia~cluster, method="bh")
s5<-ggplot(metadataCHILD, aes(cluster, Akkermansia*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  ggtitle("Akkermansia")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#313695"))+
  xlab("")+ylab("RELATIVE ABUNDANCE")+
  annotate("text", x = 1, y = 90, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 90, label = "b", size=5, col="black") +
  annotate("text", x = 3, y = 90, label = "b", size=5, col="black") +
  annotate("text", x = 4, y = 90, label = "b", size=5, col="black") +
  ylim(c(0,90))+
  guides(fill=F);s5

kruskal.test(Clostridium~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Clostridium~cluster, method="bh")
s6<-ggplot(metadataCHILD, aes(cluster, Clostridium*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  xlab("")+ylab("")+
  ggtitle("Clostridium")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#4575b4"))+
  annotate("text", x = 1, y = 63, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 63, label = "a", size=5, col="black") +
  annotate("text", x = 3, y = 63, label = "b", size=5, col="black") +
  annotate("text", x = 4, y = 63, label = "a", size=5, col="black") +
  ylim(c(0,63))+
  guides(fill=F);s6

kruskal.test(Veillonella~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Veillonella~cluster, method="bh")
s7<-ggplot(metadataCHILD, aes(cluster, Veillonella*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  xlab("")+ylab("RELATIVE ABUNDANCE")+ggtitle("Veillonella")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#74add1"))+
  annotate("text", x = 1, y = 68, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 68, label = "a", size=5, col="black") +
  annotate("text", x = 3, y = 68, label = "b", size=5, col="black") +
  annotate("text", x = 4, y = 68, label = "c", size=5, col="black") +
  ylim(c(0,68))+
  guides(fill=F);s7

kruskal.test(Faecalibacterium~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Faecalibacterium~cluster, method="bh")
s8<-ggplot(metadataCHILD, aes(cluster, Faecalibacterium*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  xlab("")+ylab("")+
  ggtitle("Faecalibacterium")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#abd9e9"))+
  annotate("text", x = 1, y = 55, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 55, label = "b", size=5, col="black") +
  annotate("text", x = 3, y = 55, label = "b", size=5, col="black") +
  annotate("text", x = 4, y = 55, label = "c", size=5, col="black") +
  ylim(c(0,55))+
  guides(fill=F);s8

kruskal.test(Klebsiella~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Klebsiella~cluster, method="bh")
s9<-ggplot(metadataCHILD, aes(cluster, Klebsiella*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  xlab("CLUSTER")+ylab("RELATIVE ABUNDANCE")+
  ggtitle("Klebsiella")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#e0f3f8"))+
  annotate("text", x = 1, y = 75, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 75, label = "a", size=5, col="black") +
  annotate("text", x = 3, y = 75, label = "a", size=5, col="black") +
  annotate("text", x = 4, y = 75, label = "b", size=5, col="black") +
  ylim(c(0,75))+
  guides(fill=F);s9

kruskal.test(Parabacteroides~cluster, data=metadataCHILD)
dunnTest(data=metadataCHILD,Parabacteroides~cluster, method="bh")
s10<-ggplot(metadataCHILD, aes(cluster, Parabacteroides*100))+theme_bw()+
  geom_boxplot(outlier.colour="black",aes(fill="main"))+
  ggtitle("Parabacteroides")+
  xlab("CLUSTER")+ylab("")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("#ffffbf"))+
  annotate("text", x = 1, y = 77, label = "a", size=5, col="black") +
  annotate("text", x = 2, y = 77, label = "b", size=5, col="black") +
  annotate("text", x = 3, y = 77, label = "a", size=5, col="black") +
  annotate("text", x = 4, y = 77, label = "b", size=5, col="black") +
  ylim(c(0,77))+
  guides(fill=F);s10

######################################################################################
######################################################################################
######################################################################################
# FIGURE 2 ###########################################################################
######################################################################################
######################################################################################
######################################################################################

ggarrange(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, ncol=4, nrow=3,
          labels = c("A","B","C","D","E","F","G","H","I","J"),
          common.legend = F)

ggarrange(ggarrange(s1, s2, s3, s4,ncol=4,labels = c("A","B","C","D")),
          ggarrange(s5, s6,ncol=4,labels = c("E","F")),
          ggarrange(s7, s8,ncol=4,labels = c("G","H")),
          ggarrange(s9, s10, ncol=4,labels = c("I","J")),
          nrow=4,
          common.legend = F)

######################################################################################
######################################################################################
######################################################################################
# FIGURE 2K ##########################################################################
######################################################################################
######################################################################################
######################################################################################

b1<-ggplot(agg,aes(cluster,x*100,fill=class))+
  geom_bar(stat="identity",alpha=0.9)+
  xlab("")+
  ylab("Relative abundance (%)")+
  ggtitle("")+
  theme_bw() + theme(legend.title = element_text(colour="black", size=12, face="bold"),
                     legend.text = element_text(colour="black", size = 12),
                     axis.text.x  = element_text(size=12, color="black"),
                     axis.text.y  = element_text(size=12, color="black"),
                     axis.title.x  = element_text(size=14, color="black"),
                     axis.title.y  = element_text(size=14, color="black"),
                     legend.position="bottom")+
  scale_fill_manual(name="Class",values=c("black",'#a50026','#d73027','#f46d43','#fdae61',
                                          '#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1',
                                          '#4575b4','#313695','#a6d96a','#66bd63','#1a9850',
                                          '#006837'),
                    labels=c("Other","Prevotella","Lachnospira","Roseburia","Haemophilus",
                             "Blautia","Parabacteroides","Klebsiella",
                             "Faecalibacterium","Veillonella","Clostridium",
                             "Akkermansia","Ruminococcus","Bifidobacterium",
                             "Escherichia","Bacteroides"))+
  guides(fill=guide_legend(ncol=3));b1

######################################################################################
######################################################################################
######################################################################################
# ENVFIT MODELS ON ALL DATA AND PER CLUSTER
colnames(CHILD)
envfit_all<-envfit(scores(PCoA.comm.vst.blind.Mat)$sites, CHILD[,c(4:16,19:29,31:38)]);envfit_all
envfit1<-envfit(scores(PCoA_clust1)$sites, Cluster1[,c(4:16,19:29,31:38)]);envfit1
envfit2<-envfit(scores(PCoA_clust2)$sites, Cluster2[,c(4:16,19:29,31:38)]);envfit2
envfit3<-envfit(scores(PCoA_clust3)$sites, Cluster3[,c(4:16,19:29,31:38)]);envfit3
envfit4<-envfit(scores(PCoA_clust4)$sites, Cluster4[,c(4:16,19:29,31:38)]);envfit4

c(c(envfit_all$vectors$pvals[which(envfit_all$vectors$pvals<0.05)],
  envfit_all$factors$pvals[which(envfit_all$factors$pvals<0.05)]),
c(envfit1$vectors$pvals[which(envfit1$vectors$pvals<0.05)],
  envfit1$factors$pvals[which(envfit1$factors$pvals<0.05)]),
c(envfit2$vectors$pvals[which(envfit2$vectors$pvals<0.05)],
  envfit2$factors$pvals[which(envfit2$factors$pvals<0.05)]),
c(envfit3$vectors$pvals[which(envfit3$vectors$pvals<0.05)],
  envfit3$factors$pvals[which(envfit3$factors$pvals<0.05)]),
c(envfit4$vectors$pvals[which(envfit4$vectors$pvals<0.05)],
  envfit4$factors$pvals[which(envfit4$factors$pvals<0.05)]))->to_plot

colnames(CHILD)[c(unique(sort(match(names(to_plot),colnames(CHILD)))))]
c(unique(sort(match(names(to_plot),colnames(CHILD)))))
tp<-colnames(CHILD)[c(4:16,19:29,31:38)];length(tp)

# ALL
nv<-length(envfit_all$vectors$r);nv
nf<-length(envfit_all$factors$r);nf
all<-matrix(nrow=nv+nf,ncol=2);all
all[1:nv,1]<-envfit_all$vectors$r;all
all[1:nv,2]<-envfit_all$vectors$pvals;all
all[(nv+1):(nv+nf),1]<-envfit_all$factors$r;all
all[(nv+1):(nv+nf),2]<-envfit_all$factors$pvals;all
all<-as.data.frame(all)
row.names(all)[1:nv]<-names(envfit_all$vectors$r)
row.names(all)[(nv+1):(nv+nf)]<-names(envfit_all$factors$r)
colnames(all)<-c("r2","pval")
all;dim(all)

# CLUSTER 1
nv<-length(envfit1$vectors$r)
nf<-length(envfit1$factors$r)
clust1<-matrix(nrow=nv+nf,ncol=2);clust1
clust1[1:nv,1]<-envfit1$vectors$r;clust1
clust1[1:nv,2]<-envfit1$vectors$pvals;clust1
clust1[(nv+1):(nv+nf),1]<-envfit1$factors$r;clust1
clust1[(nv+1):(nv+nf),2]<-envfit1$factors$pvals;clust1
clust1<-as.data.frame(clust1)
row.names(clust1)[1:nv]<-names(envfit1$vectors$r)
row.names(clust1)[(nv+1):(nv+nf)]<-names(envfit1$factors$r)
colnames(clust1)<-c("r2","pval")

# CLUSTER 2
nv<-length(envfit2$vectors$r)
nf<-length(envfit2$factors$r)
clust2<-matrix(nrow=nv+nf,ncol=2);clust2
clust2[1:nv,1]<-envfit2$vectors$r;clust2
clust2[1:nv,2]<-envfit2$vectors$pvals;clust2
clust2[(nv+1):(nv+nf),1]<-envfit2$factors$r;clust2
clust2[(nv+1):(nv+nf),2]<-envfit2$factors$pvals;clust2
clust2<-as.data.frame(clust2)
row.names(clust2)[1:nv]<-names(envfit2$vectors$r)
row.names(clust2)[(nv+1):(nv+nf)]<-names(envfit2$factors$r)
colnames(clust2)<-c("r2","pval")

# CLUSTER 3
nv<-length(envfit3$vectors$r)
nf<-length(envfit3$factors$r)
clust3<-matrix(nrow=nv+nf,ncol=2);clust3
clust3[1:nv,1]<-envfit3$vectors$r;clust3
clust3[1:nv,2]<-envfit3$vectors$pvals;clust3
clust3[(nv+1):(nv+nf),1]<-envfit3$factors$r;clust3
clust3[(nv+1):(nv+nf),2]<-envfit3$factors$pvals;clust3
clust3<-as.data.frame(clust3)
row.names(clust3)[1:nv]<-names(envfit3$vectors$r)
row.names(clust3)[(nv+1):(nv+nf)]<-names(envfit3$factors$r)
colnames(clust3)<-c("r2","pval")

# CLUSTER 4
nv<-length(envfit4$vectors$r)
nf<-length(envfit4$factors$r)
clust4<-matrix(nrow=nv+nf,ncol=2);clust4
clust4[1:nv,1]<-envfit4$vectors$r;clust4
clust4[1:nv,2]<-envfit4$vectors$pvals;clust4
clust4[(nv+1):(nv+nf),1]<-envfit4$factors$r;clust4
clust4[(nv+1):(nv+nf),2]<-envfit4$factors$pvals;clust4
clust4<-as.data.frame(clust4)
row.names(clust4)[1:nv]<-names(envfit4$vectors$r)
row.names(clust4)[(nv+1):(nv+nf)]<-names(envfit4$factors$r)
colnames(clust4)<-c("r2","pval")

all;dim(all);all$var[which(all$pval<0.05)]
clust1;dim(clust1);clust1$var[which(clust1$pval<0.05)]
clust2;dim(clust2);clust2$var[which(clust2$pval<0.05)]
clust3;dim(clust3);clust3$var[which(clust3$pval<0.05)]
clust4;dim(clust4);clust4$var[which(clust4$pval<0.05)]

all$group<-"All"
clust1$group<-"Cluster 1"
clust2$group<-"Cluster 2"
clust3$group<-"Cluster 3"
clust4$group<-"Cluster 4"

all$var<-row.names(all)
clust1$var<-row.names(clust1)
clust2$var<-row.names(clust2)
clust3$var<-row.names(clust3)
clust4$var<-row.names(clust4)

for_barplot<-rbind(all,clust1,clust2,clust3,clust4)
for_barplot$var<-as.factor(for_barplot$var)
levels(for_barplot$var)
for_barplot$var<-factor(for_barplot$var,
                        levels=levels(for_barplot$var)[rev(c(1,2,3,6,15,30,31,32,
                                                             4,5,10,11,13,28,29,
                                                             19,14,17,18,23,25,
                                                             9,16,24,26,
                                                             7,12,21,27,
                                                             8,20,22))])

titles<-c(All="ALL","Cluster 1"="CLUSTER 1","Cluster 2"="CLUSTER 2","Cluster 3"="CLUSTER 3","Cluster 4"="CLUSTER 4")

######################################################################################
######################################################################################
######################################################################################
# FIGURE S3 ##########################################################################
######################################################################################
######################################################################################
######################################################################################

(p<-ggplot(data=for_barplot, aes(x=var, y=r2,fill=var)) +
    geom_bar(stat="identity")+
    facet_wrap(~group,strip.position="top",nrow=1,ncol=5,labeller=labeller(group=titles))+
    xlab("")+
    #geom_text(data=dat_text,mapping=aes(x=var, y=y, label=label),size=5)+
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(size=12, face="bold"),
          axis.text.y  = element_text(size=10, color="black"),
          axis.text.x  = element_text(size=8, color="black",angle = 45))+
    scale_fill_manual(values=c("#980043","#dd1c77","#df65b0","#006d2c","#31a354",
                               "#74c476","#a1d99b","yellow",'#ffeda0','#edf8b1','#c7e9b4',
                               '#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494',
                               '#081d58','#fed976','#feb24c','#fd8d3c',
                               '#fc4e2a','#e31a1c','#bd0026','#67000d','#8c6bb1',
                               '#88419d','#810f7c','#4d004b','#737373','#525252',
                               '#252525','#000000'))+
    scale_x_discrete(label=rev(c("Added sugar","ASB in beverages", "ASB in sodas","Calories",
                                 "HEI2010","Sugar in beverages","Sugar in sodas","Total sugar",
                                 "BF at 3M","BF duration","Diet at 3M", "Diet at 6M","FF at 3M",
                                 "Solids at 3M","Solids at 6M",
                                 "Secretor status","Gestational diab.","Age","BMI","Ethnicity","Education",
                                 "Birth mode","Secretor status","Age","Sex",
                                 "Cat","Dog","Older siblings","Site",
                                 "Infant 6-12M ABX","Mother Intrapartum ABX","Mother Prebirth ABX")))+
    ylab(lab="UNIVARIATE R2")+
    coord_flip()+
    guides(fill=F))

# Selection for Figure 3
tp<-colnames(CHILD)[c(4,9,10,11,13,19,28,36,37,38)]
envfit_all<-envfit(scores(PCoA.comm.vst.blind.Mat)$sites, CHILD[,tp]);envfit_all
envfit1<-envfit(scores(PCoA_clust1)$sites, Cluster1[,tp]);envfit1
envfit2<-envfit(scores(PCoA_clust2)$sites, Cluster2[,tp]);envfit2
envfit3<-envfit(scores(PCoA_clust3)$sites, Cluster3[,tp]);envfit3
envfit4<-envfit(scores(PCoA_clust4)$sites, Cluster4[,tp]);envfit4

c(c(envfit_all$vectors$pvals[which(envfit_all$vectors$pvals<0.05)],
    envfit_all$factors$pvals[which(envfit_all$factors$pvals<0.05)]),
  c(envfit1$vectors$pvals[which(envfit1$vectors$pvals<0.05)],
    envfit1$factors$pvals[which(envfit1$factors$pvals<0.05)]),
  c(envfit2$vectors$pvals[which(envfit2$vectors$pvals<0.05)],
    envfit2$factors$pvals[which(envfit2$factors$pvals<0.05)]),
  c(envfit3$vectors$pvals[which(envfit3$vectors$pvals<0.05)],
    envfit3$factors$pvals[which(envfit3$factors$pvals<0.05)]),
  c(envfit4$vectors$pvals[which(envfit4$vectors$pvals<0.05)],
    envfit4$factors$pvals[which(envfit4$factors$pvals<0.05)]))->to_plot

# ALL
nv<-length(envfit_all$vectors$r);nv
nf<-length(envfit_all$factors$r);nf
all<-matrix(nrow=nv+nf,ncol=2);all
all[1:nv,1]<-envfit_all$vectors$r;all
all[1:nv,2]<-envfit_all$vectors$pvals;all
all[(nv+1):(nv+nf),1]<-envfit_all$factors$r;all
all[(nv+1):(nv+nf),2]<-envfit_all$factors$pvals;all
all<-as.data.frame(all)
row.names(all)[1:nv]<-names(envfit_all$vectors$r)
row.names(all)[(nv+1):(nv+nf)]<-names(envfit_all$factors$r)
colnames(all)<-c("r2","pval")
all;dim(all)

# CLUSTER 1
nv<-length(envfit1$vectors$r)
nf<-length(envfit1$factors$r)
clust1<-matrix(nrow=nv+nf,ncol=2);clust1
clust1[1:nv,1]<-envfit1$vectors$r;clust1
clust1[1:nv,2]<-envfit1$vectors$pvals;clust1
clust1[(nv+1):(nv+nf),1]<-envfit1$factors$r;clust1
clust1[(nv+1):(nv+nf),2]<-envfit1$factors$pvals;clust1
clust1<-as.data.frame(clust1)
row.names(clust1)[1:nv]<-names(envfit1$vectors$r)
row.names(clust1)[(nv+1):(nv+nf)]<-names(envfit1$factors$r)
colnames(clust1)<-c("r2","pval")

# CLUSTER 2
nv<-length(envfit2$vectors$r)
nf<-length(envfit2$factors$r)
clust2<-matrix(nrow=nv+nf,ncol=2);clust2
clust2[1:nv,1]<-envfit2$vectors$r;clust2
clust2[1:nv,2]<-envfit2$vectors$pvals;clust2
clust2[(nv+1):(nv+nf),1]<-envfit2$factors$r;clust2
clust2[(nv+1):(nv+nf),2]<-envfit2$factors$pvals;clust2
clust2<-as.data.frame(clust2)
row.names(clust2)[1:nv]<-names(envfit2$vectors$r)
row.names(clust2)[(nv+1):(nv+nf)]<-names(envfit2$factors$r)
colnames(clust2)<-c("r2","pval")

# CLUSTER 3
nv<-length(envfit3$vectors$r)
nf<-length(envfit3$factors$r)
clust3<-matrix(nrow=nv+nf,ncol=2);clust3
clust3[1:nv,1]<-envfit3$vectors$r;clust3
clust3[1:nv,2]<-envfit3$vectors$pvals;clust3
clust3[(nv+1):(nv+nf),1]<-envfit3$factors$r;clust3
clust3[(nv+1):(nv+nf),2]<-envfit3$factors$pvals;clust3
clust3<-as.data.frame(clust3)
row.names(clust3)[1:nv]<-names(envfit3$vectors$r)
row.names(clust3)[(nv+1):(nv+nf)]<-names(envfit3$factors$r)
colnames(clust3)<-c("r2","pval")

# CLUSTER 4
nv<-length(envfit4$vectors$r)
nf<-length(envfit4$factors$r)
clust4<-matrix(nrow=nv+nf,ncol=2);clust4
clust4[1:nv,1]<-envfit4$vectors$r;clust4
clust4[1:nv,2]<-envfit4$vectors$pvals;clust4
clust4[(nv+1):(nv+nf),1]<-envfit4$factors$r;clust4
clust4[(nv+1):(nv+nf),2]<-envfit4$factors$pvals;clust4
clust4<-as.data.frame(clust4)
row.names(clust4)[1:nv]<-names(envfit4$vectors$r)
row.names(clust4)[(nv+1):(nv+nf)]<-names(envfit4$factors$r)
colnames(clust4)<-c("r2","pval")

all;dim(all);all$var[which(all$pval<0.05)]
clust1;dim(clust1);clust1$var[which(clust1$pval<0.05)]
clust2;dim(clust2);clust2$var[which(clust2$pval<0.05)]
clust3;dim(clust3);clust3$var[which(clust3$pval<0.05)]
clust4;dim(clust4);clust4$var[which(clust4$pval<0.05)]

all$group<-"All"
clust1$group<-"Cluster 1"
clust2$group<-"Cluster 2"
clust3$group<-"Cluster 3"
clust4$group<-"Cluster 4"

all$var<-row.names(all)
clust1$var<-row.names(clust1)
clust2$var<-row.names(clust2)
clust3$var<-row.names(clust3)
clust4$var<-row.names(clust4)

for_barplot<-rbind(all,clust1,clust2,clust3,clust4)
for_barplot$var<-as.factor(for_barplot$var)
levels(for_barplot$var)
for_barplot$var<-factor(for_barplot$var,
                        levels=levels(for_barplot$var)[rev(c(1,5,6,7,9,3,2,4,8,10))])
levels(for_barplot$var)

######################################################################################
######################################################################################
######################################################################################
# FIGURE 3A ###########################################################################
######################################################################################
######################################################################################
######################################################################################

# To decide colors from main palette
a <- c("#980043","#dd1c77","#df65b0","#006d2c","#31a354",
       "#74c476","#a1d99b","yellow",'#edf8b1','#c7e9b4',
       '#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494',
       '#081d58','#ffeda0','#fed976','#feb24c','#fd8d3c',
       '#fc4e2a','#e31a1c','#bd0026','#67000d','#8c6bb1',
       '#88419d','#810f7c','#4d004b','#737373','#525252',
       '#252525','#000000')
cols(a)

(p<-ggplot(data=for_barplot, aes(x=var, y=r2,fill=var)) +
  geom_bar(stat="identity")+
  facet_wrap(~group,strip.position="top",nrow=1,ncol=5,labeller=labeller(group=titles))+
  xlab("")+
    #geom_text(data=dat_text,mapping=aes(x=var, y=y, label=label),size=5)+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size=12, face="bold"),
        axis.text.y  = element_text(size=10, color="black"),
        axis.text.x  = element_text(size=8, color="black",angle = 45))+
  scale_fill_manual(values=rev(a[c(27,22,20,18,8,6,4,12,14,15)]))+
  scale_x_discrete(label=rev(c("ASB","BMI", "SECRETOR STATUS", "INTRAPARTUM ABX",
                           "ETHNICITY","BIRTH MODE","BF at 3M","SECRETOR STATUS","OLDER SIBLINGS",
                           "INFANT AGE")))+
  ylab(lab="UNIVARIATE R2")+
  coord_flip()+
  guides(fill=F))


t0<-ggplot(CHILD, aes(Sample_time, fill=Sample_time))+theme_bw()+
  geom_bar(stat="count",color="black",width = 0.7)+
  xlab("")+ylab("COUNT")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("grey","grey"))+
  ylim(c(0,100))+
  guides(fill=F);t0
t1<-ggplot(Cluster1, aes(Sample_time, fill=Sample_time))+theme_bw()+
  geom_bar(stat="count",color="black",width = 0.7)+
  xlab("")+ylab("")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("grey","grey"))+
  ylim(c(0,50))+
  xlim("3","12")+
  guides(fill=F);t1
t2<-ggplot(Cluster2, aes(Sample_time, fill=Sample_time))+theme_bw()+
  geom_bar(stat="count",color="black",width = 0.7)+
  xlab("")+ylab("")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("grey","grey"))+
  xlim("3","12")+
  ylim(c(0,50))+
  guides(fill=F);t2
t3<-ggplot(Cluster3, aes(Sample_time, fill=Sample_time))+theme_bw()+
  geom_bar(stat="count",color="black",width = 0.7)+
  xlab("")+ylab("")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("grey","grey"))+
  ylim(c(0,50))+
  guides(fill=F);t3
t4<-ggplot(Cluster4, aes(Sample_time, fill=Sample_time))+theme_bw()+
  geom_bar(stat="count",color="black",width = 0.7)+
  xlab("")+ylab("")+
  theme(axis.text.x  = element_text(size=12, color="black"),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_text(size=14, color="black"),
        axis.title.y  = element_text(size=14, color="black"))+
  scale_fill_manual(values=c("grey","grey"))+
  ylim(c(0,50))+
  guides(fill=F);t4

ggarrange(t0, t1, t2, t3, t4, ncol=5, nrow=1,
          labels = c("","","","",""),
          legend = F,
          font.label = list(size = 20))

######################################################################################
######################################################################################
######################################################################################
# TABLE 1 ############################################################################
######################################################################################
######################################################################################
######################################################################################

# TABLE 1

model_test<-adonis(vst.blind.Mat ~ Sample_time+race_mom4+BF_3m_status+Mother_abs_birth_yn+csec+
                     older_sibs3+infant_sec+AS_bev4, data=CHILD);model_test
model_test<-adonis(vst.blind.Mat ~ Sample_time+race_mom4+BF_3m_status+Mother_abs_birth_yn+
                     older_sibs3+infant_sec+csec*AS_bev4, data=CHILD);model_test
model_test1<-adonis(vst.blind.Mat1 ~ Sample_time+race_mom4+BF_3m_status+Mother_abs_birth_yn+csec+
                     older_sibs3+infant_sec+AS_bev4, data=Cluster1);model_test1
model_test2<-adonis(vst.blind.Mat2 ~ race_mom4+BF_3m_status+Mother_abs_birth_yn+
                     older_sibs3+infant_sec+csec*AS_bev4, data=Cluster2);model_test2
model_test3<-adonis(vst.blind.Mat3 ~ Sample_time+race_mom4+BF_3m_status+Mother_abs_birth_yn+csec+
                     older_sibs3+infant_sec+AS_bev4, data=Cluster3);model_test3
model_test4<-adonis(vst.blind.Mat4 ~ Sample_time+race_mom4+BF_3m_status+Mother_abs_birth_yn+csec+
                     older_sibs3+infant_sec+AS_bev4, data=Cluster4);model_test4

######################################################################################
# LINEAR MODELS
# Here I also tested the influence of the first two axes of the ordination on
# microbial community structure but nothing showed up

contrasts(CHILD$AS_bev4) <- contr.treatment
contrasts(CHILD$csec) <- contr.treatment
contrasts(CHILD$BF_3m_status) <- contr.treatment

######################################################################################
######################################################################################
######################################################################################
# TABLE 2 ############################################################################
######################################################################################
######################################################################################
######################################################################################

# TABLE 2
CHILD3M<-CHILD[which(CHILD$Sample_time=="3"),];dim(CHILD3M)
CHILD12M<-CHILD[which(CHILD$Sample_time=="12"),];dim(CHILD12M)

# Test for normality of microbiome PCoA axes and transform accordingly
shapiro.test(CHILD3M$MDS1) #NS
shapiro.test(CHILD3M$MDS2) #p=0.025
bestNormalize(CHILD3M$MDS2)

bestNormalize(CHILD12M$MDS1)
shapiro.test((CHILD12M$MDS1+1.188886)^0.5)
CHILD12M$MDS1_trans<-(CHILD12M$MDS1+1.188886)^0.5

bestNormalize(CHILD12M$MDS2)
orderNorm_obj<-orderNorm(CHILD12M$MDS2)
hist(orderNorm_obj$x.t)
shapiro.test(orderNorm_obj$x.t)
CHILD12M$MDS2_trans<-orderNorm_obj$x.t

#3M
full.model3M <- lm(BMIz_1y~AS_bev4,data=CHILD3M);summary(full.model3M)
Anova(full.model3M,type="III")
confint(full.model3M)
full.model3M <- lm(BMIz_1y~MDS1+MDS2,data=CHILD3M);summary(full.model3M)
Anova(full.model3M,type="III")
full.model3M <- lm(BMIz_1y~AS_bev4+MDS1,data=CHILD3M);summary(full.model3M)
Anova(full.model3M,type="III")
full.model3M <- lm(BMIz_1y~AS_bev4+MDS1+MDS2,data=CHILD3M);summary(full.model3M)
Anova(full.model3M,type="III")

#12M
full.model12M <- lm(BMIz_1y~AS_bev4+Mother_abs_birth_yn+csec+CHILD3metabo$Succinic.acid,data=CHILD12M);summary(full.model12M)
Anova(full.model12M,type="III")
confint(full.model12M)
full.model12M <- lm(BMIz_1y~MDS1_trans+MDS2_trans+CHILD3metabo$Succinic.acid,data=CHILD12M);summary(full.model12M)
Anova(full.model12M,type="III")
confint(full.model12M)
full.model12M <- lm(BMIz_1y~AS_bev4+MDS1_trans+CHILD3metabo$Succinic.acid,data=CHILD12M);summary(full.model12M)
Anova(full.model12M,type="III")
confint(full.model12M)
CHILD12M$resids <- residuals(full.model12M)

full.model12M <- lm(BMIz_1y~AS_bev4,data=CHILD12M);summary(full.model12M)
full.model12M <- lm(BMIz_1y~Succinic.acid,data=CHILD12M);summary(full.model12M)
full.model12M <- lm(BMIz_1y~AS_bev4*Succinic.acid,data=CHILD12M);summary(full.model12M)

med.fit<-lm(Succinic.acid~AS_bev4, data=CHILD12M)
out.fit<-lm(BMIz_1y~Succinic.acid +AS_bev4, data=CHILD12M)
med.out <- mediate(med.fit, out.fit, treat = "AS_bev4",
                   mediator = "Succinic.acid", robustSE = TRUE, sims = 9999)
summary(med.out)

ggplot(CHILD3M, aes(Succinic.acid, BMIz_1y))+geom_point()+geom_smooth(method="lm")+xlim(c(-1,1.5))+xlab("Succinate")+ylab("BMI at 12 months")

#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(CHILD12M$resids)

######################################################################################
######################################################################################
######################################################################################
# FIGURE S3 ##########################################################################
######################################################################################
######################################################################################
######################################################################################

#Make ordination of dataset showing main taxonomical groups of bacteria with envfit
rel_abun_genus<-sweep(t(comm.taxoGenera), 2, apply(t(comm.taxoGenera),2,sum), `/`)
apply(rel_abun_genus,2,sum) # Here you confirm that all samples have a total relative abundance of 1
total_rel_abun_gen<-apply(rel_abun_genus,1,sum)
total_rel_abun_gen[order(total_rel_abun_gen)]->total_rel_abun_gen;length(total_rel_abun_gen) #91
subset_taxo<-comm.taxoGenera[,names(total_rel_abun_gen)[c(80:91)]]

## Make vector of colors for values smaller than 0 (20 colors)
rc1 <- colorRampPalette(colors = c("red", "white"), space = "Lab")(10)

## Make vector of colors for values larger than 0 (180 colors)
rc2 <- colorRampPalette(colors = c("white", "green"), space = "Lab")(10)

## Combine the two color palettes
rampcols <- c(rc1, rc2)

mypal <- colorNumeric(palette = rampcols, domain = CHILD$BMIz_1y)

## If you want to preview the color range, run the following code
previewColors(colorNumeric(palette = rampcols, domain = NULL), values = -2:2)
colvec<-c("blue","yellow")
with(scores(PCoA.comm.vst.blind.Mat), points(mod, display = "sites", col = colvec[Use], pch = 19))
# Genus level PCoA on variance stabilized matrix
ordiplot(PCoA.comm.vst.blind.Mat, type="n", xlab="PCoA Axis 1", ylab="PCoA Axis 2", display="sites")
points(scores(PCoA.comm.vst.blind.Mat)$sites, pch=19, col=colvec[CHILD$Sample_time], cex=1.5)
points(scores(PCoA.comm.vst.blind.Mat)$sites, pch=1, cex=1.5)
# With labels
plot(envfit(scores(PCoA.comm.vst.blind.Mat)$sites, subset_taxo),p.max=0.01, col="black", lwd=1.5)
# No Labels
plot(envfit(scores(PCoA.comm.vst.blind.Mat)$sites, subset_taxo),p.max=0.01, col="black", lwd=1.5, labels="")

######################################################################################

genera <- t(asv / rowSums(asv))
apply(genera, 2, sum)
apply(genera, 1, sum)
top100 <- head(names(rev(sort(rowSums(genera)))), 100)
head(rev(sort(rowSums(genera))),100)/198*100
taxa[top100,]$joint
asv[,top100]->asv100;dim(asv100)

data_CHILD <- phyloseq(otu_table(asv100, taxa_are_rows=FALSE), 
                       sample_data(CHILD), tax_table(as.matrix(taxa)));data_CHILD
data_clust1 <- phyloseq(otu_table(asv100, taxa_are_rows=FALSE), 
                       sample_data(Cluster1), tax_table(as.matrix(taxa)));data_clust1
data_clust2 <- phyloseq(otu_table(asv100, taxa_are_rows=FALSE), 
                       sample_data(Cluster2), tax_table(as.matrix(taxa)));data_clust2
data_clust3 <- phyloseq(otu_table(asv100, taxa_are_rows=FALSE), 
                       sample_data(Cluster3), tax_table(as.matrix(taxa)));data_clust3
data_clust4 <- phyloseq(otu_table(asv100, taxa_are_rows=FALSE), 
                       sample_data(Cluster4), tax_table(as.matrix(taxa)));data_clust4


#DESEQ2

#on all dataset
test.phyloseq.dds<-phyloseq_to_deseq2(data_CHILD,~AS_bev4)
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds),
                                                                                       1, gm_mean))
diagdds = DESeq(test.phyloseq.dds)
res <- results(diagdds)
summary(res)
res <- res[order(res$padj),]

lall<-length(which(res$padj<0.05));lall
for (i in 1:lall){
  t<-which(row.names(taxa)==row.names(res)[i])
  print(taxa[t,]$joint)
}

deseq_dd<-matrix(ncol=5,nrow=20)
deseq_dd[1:lall,1]<-c("Bacteroides sp. ASV45",
                      "Prevotella copri ASV42")
deseq_dd[1:lall,2]<-res$log2FoldChange[1:lall]
deseq_dd[1:lall,3]<-res$padj[1:lall]
deseq_dd[1:lall,5]<-row.names(res)[1:lall]

#per cluster
#cluster 1
test.phyloseq.dds<-phyloseq_to_deseq2(data_clust1,~AS_bev4)
# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds),
                                                                            1, gm_mean))
diagdds = DESeq(test.phyloseq.dds)
res <- results(diagdds)
summary(res)
res <- res[order(res$padj),]

lc1<-length(which(res$padj<0.05));lc1
for (i in 1:lc1){
  t<-which(row.names(taxa)==row.names(res)[i])
  print(taxa[t,]$joint)
}

deseq_dd[(lall+1):(lall+lc1),1]<-c("Bacteroides uniformis ASV13",
                    "Bacteroides ovatus ASV6",
                    "Clostridium perfringens ASV90",
                    "Bacteroides ovatus ASV27",
                    "Parabacteroides sp. ASV18")

deseq_dd[(lall+1):(lall+lc1),2]<-res$log2FoldChange[1:lc1]
deseq_dd[(lall+1):(lall+lc1),3]<-res$padj[1:lc1]
deseq_dd[(lall+1):(lall+lc1),5]<-row.names(res)[1:lc1]

colnames(deseq_dd)<-c("Taxa","log2","padj","group","taxa")
deseq_dd[1:(lall+lc1),4]<-c("ALL","ALL","CLUSTER 1","CLUSTER 1","CLUSTER 1","CLUSTER 1","CLUSTER 1")

#cluster 2
test.phyloseq.dds<-phyloseq_to_deseq2(data_clust2,~AS_bev4)
# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds),
                                                                            1, gm_mean))
dds <- DESeq(test.phyloseq.dds)
res <- results(dds)
summary(res)
res <- res[order(res$padj),]

lc2<-length(which(res$padj<0.05));lc2
for (i in 1:lc2){
  t<-which(row.names(taxa)==row.names(res)[i])
  print(taxa[t,]$joint)
}

deseq_dd[(lall+lc1+1):(lall+lc1+lc2),1]<-c("Bacteroides ovatus ASV27",
                    "Parabacteroides sp. ASV83",
                    "Akkermansia muciniphila ASV19",
                    "Bacteroides sp. ASV45",
                   "Bacteroides sp. ASV25")
deseq_dd[(lall+lc1+1):(lall+lc1+lc2),2]<-res$log2FoldChange[1:lc2]
deseq_dd[(lall+lc1+1):(lall+lc1+lc2),3]<-res$padj[1:lc2]
deseq_dd[(lall+lc1+1):(lall+lc1+lc2),5]<-row.names(res)[1:lc2]
deseq_dd[(lall+lc1+1):(lall+lc1+lc2),4]<-c("CLUSTER 2","CLUSTER 2","CLUSTER 2","CLUSTER 2","CLUSTER 2")

#cluster 3
test.phyloseq.dds<-phyloseq_to_deseq2(data_clust3,~AS_bev4)
# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds),
                                                                            1, gm_mean))
dds <- DESeq(test.phyloseq.dds)
res <- results(dds)
summary(res)
res <- res[order(res$padj),]

lc3<-length(which(res$padj<0.05));lc3
for (i in 1:lc3){
  t<-which(row.names(taxa)==row.names(res)[i])
  print(taxa[t,]$joint)
}

deseq_dd[(lall+lc1+lc2+1):(lall+lc1+lc2+lc3),1]<-c("Rikenellaceae ASV24",
                                   "Enterobacteriaceae ASV14")
deseq_dd[(lall+lc1+lc2+1):(lall+lc1+lc2+lc3),2]<-res$log2FoldChange[1:lc3]
deseq_dd[(lall+lc1+lc2+1):(lall+lc1+lc2+lc3),3]<-res$padj[1:lc3]
deseq_dd[(lall+lc1+lc2+1):(lall+lc1+lc2+lc3),4]<-c("CLUSTER 3","CLUSTER 3")
deseq_dd[(lall+lc1+lc2+1):(lall+lc1+lc2+lc3),5]<-row.names(res)[1:lc3]

#cluster 4
test.phyloseq.dds<-phyloseq_to_deseq2(data_clust4,~AS_bev4)
# Estimate factor size
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds),
                                                                            1, gm_mean))
dds <- DESeq(test.phyloseq.dds)
res <- results(dds)
summary(res)
res <- res[order(res$padj),]

lc4<-length(which(res$padj<0.05));lc4
for (i in 1:lc4){
  t<-which(row.names(taxa)==row.names(res)[i])
  print(taxa[t,]$joint)
}

deseq_dd[(lall+lc1+lc2+lc3+1):(lall+lc1+lc2+lc3+lc4),1]<-c("Bacteroides caccae ASV38",
                                                           "Bacteroides uniformis ASV77",
                                                           "Bacteroides sp. ASV85",
                                                           "Bacteroides sp. ASV25",
                                                           "Bacteroides sp. ASV78",
                                                           "Escherichia coli ASV4"
)
deseq_dd[(lall+lc1+lc2+lc3+1):(lall+lc1+lc2+lc3+lc4),2]<-res$log2FoldChange[1:lc4]
deseq_dd[(lall+lc1+lc2+lc3+1):(lall+lc1+lc2+lc3+lc4),3]<-res$padj[1:lc4]
deseq_dd[(lall+lc1+lc2+lc3+1):(lall+lc1+lc2+lc3+lc4),4]<-c("CLUSTER 4","CLUSTER 4","CLUSTER 4","CLUSTER 4","CLUSTER 4","CLUSTER 4")
deseq_dd[(lall+lc1+lc2+lc3+1):(lall+lc1+lc2+lc3+lc4),5]<-row.names(res)[1:lc4]
deseq_dd<-as.data.frame(deseq_dd)
deseq_dd$group<-as.factor(deseq_dd$group)
deseq_dd$log2<-as.numeric(as.character(deseq_dd$log2))
deseq_dd$padj<-as.numeric(as.character(deseq_dd$padj))
summary(deseq_dd)
deseq_dd$Taxa

######################################################################################
######################################################################################
######################################################################################
# FIGURE 3B ##########################################################################
######################################################################################
######################################################################################
######################################################################################

h1<-ggplot(data = deseq_dd,
       aes(Taxa, log2,
           fill = log2 > 0))+
  geom_bar(stat = "identity")+
  theme_bw()+
  facet_wrap(~group,ncol=5)+
  xlab("")+
  ylab("Log2 fold change for no vs. daily maternal consumption of ASBs")+
  coord_flip()+
  theme(axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=16, color="black"),
        axis.text.x  = element_text(size=10, color="black",angle=90),
        strip.text.x = element_text(size=12, face="bold"))+
  guides(fill = FALSE);h1

######################################################################################
######################################################################################

rel_abun_asv<-sweep(t(asv), 2, apply(t(asv),2,sum), `/`)
apply(rel_abun_asv,2,sum)
rel_abun_asv<-t(rel_abun_asv)
aggregate(rel_abun_asv[,row.names(taxa)[which(row.names(taxa)==deseq_dd$taxa[1])]],
          by=list(CHILD$AS_bev4), mean)
aggregate(rel_abun_asv[,row.names(taxa)[which(row.names(taxa)==deseq_dd$taxa[2])]],
          by=list(CHILD$AS_bev4), mean)

asv1<-as.matrix(otu_table(data_clust1, taxa_are_rows = F));dim(asv1)
asv2<-as.matrix(otu_table(data_clust2, taxa_are_rows = F));dim(asv2)
asv3<-as.matrix(otu_table(data_clust3, taxa_are_rows = F));dim(asv3)
asv4<-as.matrix(otu_table(data_clust4, taxa_are_rows = F));dim(asv4)
rel_abun_asv1<-sweep(t(asv1), 2, apply(t(asv1),2,sum), `/`)
apply(rel_abun_asv1,2,sum)
rel_abun_asv1<-t(rel_abun_asv1);dim(rel_abun_asv1)
aggregate(rel_abun_asv1[,row.names(taxa)[which(row.names(taxa)==deseq_dd$taxa[3])]],
          by=list(Cluster1$AS_bev4), mean)
aggregate(rel_abun_asv[,row.names(taxa)[which(row.names(taxa)==deseq_dd$taxa[4])]],
          by=list(CHILD$AS_bev4), mean)

######################################################################################
######################################################################################

#on all dataset
test.phyloseq.dds<-phyloseq_to_deseq2(data_CHILD,~BMIz_1y)
test.phyloseq.dds = estimateSizeFactors(test.phyloseq.dds, geoMeans = apply(counts(test.phyloseq.dds),
                                                                            1, gm_mean))
diagdds = DESeq(test.phyloseq.dds)
res <- results(diagdds)
summary(res)
res <- res[order(res$padj),]

lall<-length(which(res$padj<0.05));lall
call<-rep(NA, lall);call
for (i in 1:lall){
  t<-which(row.names(taxa)==row.names(res)[i])
  call[i]<-as.character((taxa[t,]$joint))
}

deseq_dd<-matrix(ncol=5,nrow=lall)
deseq_dd[1:lall,2]<-res$log2FoldChange[1:lall]
deseq_dd[1:lall,3]<-res$padj[1:lall]
deseq_dd[1:lall,5]<-row.names(res)[1:lall]

colnames(deseq_dd)<-c("Taxa","log2","padj","group","taxa")
lall#10

deseq_dd[1:(lall),4]<-c("ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL","ALL")                                                        
deseq_dd[1:lall,1]<-c("Clostridium sp. ASV102",
                      "Pantoea agglomerans ASV115",
                      "Bacteroides eggerthii ASV61",
                      "Bacteroidesplebeius ASV124",
                      "Bacteroides sp. ASV78",
                      "Haemophiluspara influenzae ASV22",
                      "Paraprevotella sp. ASV94",
                      "Akkermansia muciniphila ASV1",
                      "Faecalibacterium prausnitzii ASV17",
                      "Clostridium perfringens ASV90")

deseq_dd<-as.data.frame(deseq_dd)
deseq_dd$group<-as.factor(deseq_dd$group)
deseq_dd$log2<-as.numeric(as.character(deseq_dd$log2))
deseq_dd$padj<-as.numeric(as.character(deseq_dd$padj))
summary(deseq_dd)
deseq_dd$Taxa

######################################################################################
######################################################################################
######################################################################################
# FIGURE 3B ##########################################################################
######################################################################################
######################################################################################
######################################################################################

h1<-ggplot(data = deseq_dd,
           aes(Taxa, log2,
               fill = log2 > 0))+
  geom_bar(stat = "identity")+
  theme_bw()+
  facet_wrap(~group,ncol=5)+
  xlab("")+
  ylab("Log2 fold change for no vs. daily maternal consumption of ASBs")+
  coord_flip()+
  theme(axis.text.y = element_text(size=12, color="black"),
        axis.title.x = element_text(size=16, color="black"),
        axis.text.x  = element_text(size=10, color="black",angle=90),
        strip.text.x = element_text(size=12, face="bold"))+
  guides(fill = FALSE);h1

