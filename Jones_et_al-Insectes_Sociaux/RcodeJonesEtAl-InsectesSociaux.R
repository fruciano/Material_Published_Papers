#### This script briefly outlines the analyses used in
#### Jones et al. 2018 - Insectes Sociaux 65(3): 419-429.
#### Each of these analyses was performed multiple times
#### (e.g., on multiple rarefied datasets, 
#### in multiple pairwise comparisons,
#### please, see the text).
#### Here, some of these repeated analyses are
#### not provided (see comments)
#### to keep the script short and easy to follow.

#### If you use this script in your own research
#### a citation to
#### Jones et al. 2018 - Insectes Sociaux 65(3): 419-429
#### would be kindly appreciated
#### However, more importantly, please cite the authors
#### whom have developed the methods we have used and
#### their R implementations (see the text of the paper
#### for citations)


library(vegan)
library(ape)
library(phyloseq)
library(phangorn)
library(ggplot2)
library(DESeq2)
library(sgof)



################################################
#### 				Prelims					####
#### (import data, rarefaction and similar) ####
################################################


Reads=read.table("Species.txt",row.names=1,header=T,sep="\t")
Reads=t(Reads)
# Read and transpose the matrix with the counts
# this is one of the outputs of the LotuS pipeline
# As it is often the case with this kind of software,
# this is a matrix of counts with OTUs in rows and
# individual bees in columns
# (this is why it gets transposed at the beginning)
# the individual names match the names in the SRA

factors=read.table("factorsNEW.txt",sep="\t",header=T,row.names=1)
# Import factors
# This is a tab-delimited file with individuals in rows
# and columns ID, Metadata (which is the behavioural group), and Apiary
# the ID matches the name of individuals in the file "Species.txt"

AllMerged=merge(factors,Reads,by="row.names")
rownames(AllMerged)=AllMerged[,1]
# Merge reads and factors

BehaviouralRareCurveByType=rarecurve(AllMerged[,4:ncol(AllMerged)],
	step=100, xlab = "Reads", ylab = "OTUs",
	label = TRUE, col=AllMerged[,3])
# Plot rarefaction curves by individual, color-code by group

BehaviouralMerged=AllMerged[-which(AllMerged[,1]=="APM59"),]
# Removes APM59 which has very few sequences


BehaviouralRarefy=rrarefy(BehaviouralMerged[,4:ncol(BehaviouralMerged)],
	min(apply(BehaviouralMerged[,4:ncol(BehaviouralMerged)],1,sum)) )
BehaviouralRarefy=cbind(BehaviouralMerged[,1:3],BehaviouralRarefy)
# Rarefy to the smallest number of reads
# Notice that for the main analyses only a single rarefied sample has been used
# but multiple rarefied samples have been produced to check for consistency
# (see text)


###############################################
#### Analyses on Bray-Curtis dissimilarity ####
###############################################

MDSbehav=metaMDS(BehaviouralRarefy[,4:ncol(BehaviouralRarefy)],k=11)
# non-metric multidimensional scaling
# k empirically adjusted to obtain sufficiently low stress

plotlims=range(c(MDSbehav$points[,1],MDSbehav$points[,2]))+c(-0.15,0.15)
qplot(MDSbehav$points[,1],MDSbehav$points[,2],asp=1,
color=BehaviouralRarefy[,2],shape=BehaviouralRarefy[,2],
xlim=plotlims,ylim=plotlims,
xlab="MDS axis 1",ylab="MDS axis 2")+
theme_classic()+ theme(aspect.ratio = 1)+ geom_point(size = 4)+
theme(legend.title=element_blank())+
stat_ellipse()
# Plot of scores along the first two
# MDS axes

# Notice that this set of analyses has also been repeated
# for multiple rarefied samples (see text)

MDSbehavD=metaMDSdist(BehaviouralRarefy[,4:ncol(BehaviouralRarefy)])
Betadispersion=betadisper(MDSbehavD, BehaviouralRarefy[,2])
	plot(Betadispersion)
	table(BehaviouralRarefy[,2])
	anova(Betadispersion)
	boxplot(Betadispersion)
TukeyHSD(Betadispersion)
# Analyses of dispersion

PERMANOVABrayBehararBehavAndColony1=adonis(
	MDSbehavD~BehaviouralMerged[,2]/BehaviouralMerged[,3])
# PERMANOVA
# Notice that this was also performed pairwise among groups
# (see text)



#################################################
#### Analyses on the Shannon diversity index ####
#################################################

BehaviouralShannon=cbind(BehaviouralRarefy[,1:3],
	vegan::diversity(BehaviouralRarefy[,4:ncol(BehaviouralRarefy)],index="shannon"))
colnames(BehaviouralShannon)=c("Individual","Behaviour","Apiary","Shannon")
# Perform analyses of Shannon diversity index
# this analysis has been repeated for multiple rarefied
# samples to check for consistency (see text)

ggplot(BehaviouralShannon,
		aes(x=BehaviouralShannon[,2],
		y=BehaviouralShannon[,4],
		fill=BehaviouralShannon[,2])) +
	geom_boxplot()+ theme_classic()
# Plot

anovaBehav=aov(Shannon~Behaviour, data=BehaviouralShannon)
summary(anovaBehav)
TukeyHSD(anovaBehav)
# Comparisons


###############################################
####  Analyses based on UniFrac distances  ####
###############################################

OTU=read.table("OTU.txt",row.names=1,header=T,sep="\t")
OTU=t(OTU)
# Read the OTU table and transpose it
# This is a table obtained from the LotuS pipeline
# with the same structure as "Species.txt"

BehaviouralOTUMerged=merge(factors,OTU,by="row.names")
BehaviouralOTUMerged=BehaviouralOTUMerged[-which(BehaviouralOTUMerged[,1]=="APM59"),]
# Removes APM59 which has very few sequences

BehaviouralOTURarefy=rrarefy(
	BehaviouralOTUMerged[,4:ncol(BehaviouralOTUMerged)],
	min(apply(BehaviouralOTUMerged[,4:ncol(BehaviouralOTUMerged)],1,sum))
	)
BehaviouralOTURarefy=cbind(BehaviouralOTUMerged[,1:3],BehaviouralOTURarefy)
# Rarefaction
# Please, notice this has been repeated multiple times to check for consistency
# (see text)

tree=read.tree('tree.nwk')
# This is the tree obtained from the LotuS pipeline
droplist<-setdiff(tree$tip, colnames(OTU))
tree<-drop.tip(tree,droplist)
# Import tree and prune it

source("midpoint.R")
# Source function to midpoint-root
# a tree, by Klaus Schliep
# This can be found at
# https://www.mail-archive.com/r-sig-phylo@r-project.org/msg03534.html


tree=midpoint(tree)
# Midpoint-root the tree

BehaviouralOTUrar=BehaviouralOTURarefy[,4:ncol(BehaviouralOTURarefy)]
BehaviouralOTUtablerar=otu_table(BehaviouralOTUrar,taxa_are_rows=F)
phyloseqBehaviourrar=phyloseq(BehaviouralOTUtablerar,tree)
# Prepare data for phyloseq and create objects of class phyloseq

UnifracBehaviourrar=distance(phyloseqBehaviourrar,method="unifrac")
UnifracWBehaviourrar=distance(phyloseqBehaviourrar,method="wunifrac")
# Compute UniFrac and weighted UniFrac distances among samples

PERMANOVAunifracBehaviourAndColonyrar=adonis(
	UnifracBehaviourrar~BehaviouralOTUMerged[,2]/BehaviouralOTUMerged[,3])
PERMANOVAwunifracBehaviourAndColonyrar=adonis(
	UnifracWBehaviourrar~BehaviouralOTUMerged[,2]/BehaviouralOTUMerged[,3])
# PERMANOVA
# Notice that this was also performed pairwise among groups
# (see text)

BetadispersionUnifrac=betadisper(UnifracBehaviourrar, BehaviouralOTUMerged[,2])
	plot(BetadispersionUnifrac)
	table(BehaviouralOTUMerged[,2])
	anova(BetadispersionUnifrac)
	boxplot(BetadispersionUnifrac)
TukeyHSD(BetadispersionUnifrac)

BetadispersionWUnifrac=betadisper(UnifracWBehaviourrar, BehaviouralOTUMerged[,2])
	plot(BetadispersionWUnifrac)
	table(BehaviouralOTUMerged[,2])
	anova(BetadispersionUnifrac)
	boxplot(BetadispersionUnifrac)
TukeyHSD(BetadispersionUnifrac)
# Analyses of dispersion


################################################
####  Analyses based on the ANCOM procedure ####
################################################

library("ancom.R")
### This assumes that the ANCOM R implementation
### has been installed in the system

BehaviouralRarefy4ANCOM=cbind(BehaviouralRarefy[,4:ncol(BehaviouralRarefy)],BehaviouralRarefy$Metadata)

BehaviouralRarefyANCOM=ANCOM(BehaviouralRarefy4ANCOM,multcorr=1)
plot_ancom(BehaviouralRarefyANCOM)
BehaviouralRarefyANCOM$detected
# Run the ANCOM analysis


################################################
####	  Analyses based on DESeq2 			####
################################################

# Notice that the analyses below are based on
# a single comparison (FOR vs NUR)
# they have been repeated for each comparison
# (see text)

Counts4Deseq=as.matrix(t(
	BehaviouralMerged[
		which(BehaviouralMerged$Metadata=="FOR" | BehaviouralMerged$Metadata=="NUR"),
			4:ncol(BehaviouralMerged)]
	))

Counts4Deseq=Counts4Deseq[which(apply(Counts4Deseq,1,sum)>500),]
# Remove taxa with less than 500 reads

Factors4Deseq=data.frame(row.names=colnames(Counts4Deseq), 
	BehaviouralMerged[
		which(BehaviouralMerged$Metadata=="FOR" | BehaviouralMerged$Metadata=="NUR"),
			"Metadata"])
colnames(Factors4Deseq)="BehaviouralMerged.Metadata"
BehaviourDESeq2Obj=DESeqDataSetFromMatrix(Counts4Deseq, Factors4Deseq, ~BehaviouralMerged.Metadata)
# Prepare data for DESeq2

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(BehaviourDESeq2Obj), 1, gm_mean)
diagdds = estimateSizeFactors(BehaviourDESeq2Obj, geoMeans = geoMeans)
DESeq2diagdds = DESeq(diagdds, fitType="local")


res = results(DESeq2diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(DESeq2diagdds )[rownames(sigtab), ], "matrix"))
head(sigtab)
# Perform DESeq2 analysis

resultsALL=cbind(as(res, "data.frame"))
resultsALL=resultsALL[order(resultsALL$pvalue), ]
BH=BH(resultsALL$pvalue)
Results=cbind(resultsALL,BH$Adjusted.pvalues)
# Get the results and adjust for multiple
# comparisons using the
# Benjamini-Hochberg approach
