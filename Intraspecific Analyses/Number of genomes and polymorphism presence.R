library(readxl)
library(ape)
library(phylolm)

#Read data sets
mam.data <- read_excel("Dataset S2 - Intraspecific Variation.xlsx", sheet = "Intraspecific_Analysis1") ##N.Genomes only

# Read tree block
mam.tree<-read.nexus("~/Downloads/tree_VertLife2-27-25.nex")
name.check(mam.tree[[1]], mam.data)

###making a 0/1 for variation present in N.genomes
mam.data$Variation.01<-0
mam.data$Variation.01[mam.data$Variation=="yes"]<-1

###removing tips without multiple genomes
missing<-name.check(mam.tree[[1]],mam.syn)$tree_not_data
mam.tree4<-drop.tip.multiPhylo(mam.tree,missing)

name.check(mam.tree4[[1]]$tip.labels,mam.data)

#Create data frame to store results of phyloglm on each tree in the block
phy.n.individuals<-data.frame(tree=1:1000,Intercept=0,Intercept.SD=0,Intercept.Z=0,Intercept.P=0,N.Genomes=0,N.Genomes.SD=0,N.Genomes.Z=0,N.Genomes.P=0,alpha=0)

#Phylogenetic glm
for (i in 1:length(mam.tree4)){
  polymorph<-phyloglm(Variation.01~N.Genomes, phy=mam.tree4[[i]],data=mam.data, method="logistic_MLE", btol=100)
  polymorph.sum<-summary(polymorph)
  ###Write results to data frame
  phy.n.individuals$Intercept[i]=polymorph.sum$coefficients[1,1]
  phy.n.individuals$Intercept.SD[i]=polymorph.sum$coefficients[1,2]
  phy.n.individuals$Intercept.Z[i]=polymorph.sum$coefficients[1,3]
  phy.n.individuals$Intercept.P[i]=polymorph.sum$coefficients[1,4]
  phy.n.individuals$N.Genomes[i]=polymorph.sum$coefficients[2,1]
  phy.n.individuals$N.Genomes.SD[i]=polymorph.sum$coefficients[2,2]
  phy.n.individuals$N.Genomes.Z[i]=polymorph.sum$coefficients[2,3]
  phy.n.individuals$N.Genomes.P[i]=polymorph.sum$coefficients[2,4]
  phy.n.individuals$alpha[i]=polymorph.sum$alpha
}

write.csv(phy.n.individuals,"N.genomes_phylo_results.csv")

###plotting results of running across the tree block
hist(phy.n.individuals$N.Genomes,breaks=seq(0,0.5,0.05))
abline(v=median(phy.n.individuals$N.Genomes))

hist(phy.n.individuals$N.Genomes, breaks = seq(0,0.5,0.05))
