library(readxl)
library(ape)
library(crayon)
library(geiger)
library(phylolm)
library(ggplot2)

#Read data set
mam.mutation <- read_excel("Dataset S2 - Intraspecific Variation.xlsx", sheet = "Intraspecific_Analysis2") ##Synonymous, nonsynonymous or indel

# Read tree block
mam.tree<-read.nexus("~/Downloads/tree_VertLife2-27-25.nex")

###removing tips without multiple genomes
missing<-name.check(mam.tree[[1]],subset(mam.mutation, Type == "1-synonymous"))$tree_not_data
mam.tree.multi<-drop.tip.multiPhylo(mam.tree,missing)

###removing tips without polymorphisms
no.variation<-rownames(subset(subset(mam.mutation, Total.Polymorphisms==0), Type == "1-synonymous"))
mam.tree.polymorph<-drop.tip.multiPhylo(mam.tree.multi,no.variation)
tips<-mam.tree.polymorph[[1]]$tip.label

##Adding tips for non-synonymus and indel mutations
###Adding non-synonymous tips
  treeBlock<-mam.tree.polymorph
  mam.tree5<-list()
  for (i in 1:1000) {
    #Identify the position (i.e., integer number) where insertNextT is located in the tree object
    tree.run<-treeBlock[[i]]
    for (j in 1: 43) {
      insertSpeciesNextTo = treeBlock[[i]]$tip.label[j]
      speciesToBeInserted = treeBlock[[i]]$tip.label[j]%+%"_nonsyn"
      branchLengthSplitAgo = 0.000001 ###
      treeNewStr = paste("(", speciesToBeInserted, ":", branchLengthSplitAgo, ");", sep="")
      treeNew <- read.tree(text = treeNewStr)
      pos = j
      #Merge both trees and write tree to new list of trees
      tree.run = bind.tree(tree.run, treeNew,where = pos, position = branchLengthSplitAgo)
    }
    mam.tree5[[i]]<-tree.run
  }
  #Write new nexus file    
    filename = "mam.tree5.nex"
    write.nexus(mam.tree5,file = filename)
    mam.tree5<-read.nexus("~/Downloads/mam.tree5.nex")
    
  ###Adding indel tips
  treeBlock2<-mam.tree5
  mam.tree5<-list()
  for (i in 1:1000) {
    #Identify the position (i.e., integer number) where insertNextT is located in the tree object
    tree.run<-treeBlock2[[i]]
    for (j in 1: 43) {
      insertSpeciesNextTo = treeBlock2[[i]]$tip.label[j]
      speciesToBeInserted = treeBlock2[[i]]$tip.label[j]%+%"_indel"
      branchLengthSplitAgo = 0.0000005 ###
      treeNewStr = paste("(", speciesToBeInserted, ":", branchLengthSplitAgo, ");", sep="")
      treeNew <- read.tree(text = treeNewStr)
      pos = j
      #Merge both trees and write tree to new list of trees
      tree.run = bind.tree(tree.run, treeNew,where = pos, position = branchLengthSplitAgo)
    }
    mam.tree5[[i]]<-tree.run
  }

#Write new nexus file
  write.nexus(mam.tree5,file = filename)
  mam.tree5<-read.nexus("~/Downloads/mam.tree5.nex")
  
  
###Making a dataset of only species with polymorphisms
mam.syn2<-subset(mam.mutation, Variation=="yes")

#Create data frame to store results of pglmm_compare on each tree in the block
phy.mutation.type.glmm<-data.frame(tree=1:1000,Intercept=0,Intercept.SE=0,Intercept.Z=0,Intercept.P=0,Nonsynonymous=0,Nonsynonymous.SE=0,Nonsynonymous.Z=0,Nonsynonymous.P=0,Indel=0,Indel.SE=0,Indel.Z=0,Indel.P=0,Random.SD=0)

name.check(mam.tree5[[1]],mam.syn2)

#Phylogenetic GLMM
for (i in 1:length(mam.tree5)){
  mutation_type_glmm<-pglmm_compare(Variants~Type, family = "poisson", phy=mam.tree5[[i]],data=mam.syn2)
  ###Write results to data frame
  phy.mutation.type.glmm$Intercept[i]=mutation_type_glmm$B[1,1]
  phy.mutation.type.glmm$Intercept.SE[i]=mutation_type_glmm$B.se[1,1]
  phy.mutation.type.glmm$Intercept.Z[i]=mutation_type_glmm$B.zscore[1,1]
  phy.mutation.type.glmm$Intercept.P[i]=mutation_type_glmm$B.pvalue[1,1]
  phy.mutation.type.glmm$Nonsynonymous[i]=mutation_type_glmm$B[2,1]
  phy.mutation.type.glmm$Nonsynonymous.SE[i]=mutation_type_glmm$B.se[2,1]
  phy.mutation.type.glmm$Nonsynonymous.Z[i]=mutation_type_glmm$B.zscore[2,1]
  phy.mutation.type.glmm$Nonsynonymous.P[i]=mutation_type_glmm$B.pvalue[2,1]
  phy.mutation.type.glmm$Indel[i]=mutation_type_glmm$B[3,1]
  phy.mutation.type.glmm$Indel.SE[i]=mutation_type_glmm$B.se[3,1]
  phy.mutation.type.glmm$Indel.Z[i]=mutation_type_glmm$B.zscore[3,1]
  phy.mutation.type.glmm$Indel.P[i]=mutation_type_glmm$B.zscore[3,1]
  phy.mutation.type.glmm$Random.SD[i]=mutation_type_glmm$ss
}

write.csv(phy.mutation.type.glmm,"Variant_type_PGLMM.csv")

#Phylogenetic GLM
phy.mutation.type<-data.frame(tree=1:1000,Intercept=0,Intercept.SE=0,Intercept.Z=0,Intercept.P=0,Nonsynonymous=0,Nonsynonymous.SE=0,Nonsynonymous.Z=0,Nonsynonymous.P=0,Indel=0,Indel.SE=0,Indel.Z=0,Indel.P=0,Scale=0)

for (i in 1:length(mam.tree5)){
  mutation_type<-phyloglm(Variants~Type, phy=mam.tree5[[i]],data=mam.syn2, method="poisson_GEE")
  mtype.sum<-summary(mutation_type)
  ###Write results to data frame
  phy.mutation.type$Intercept[i]=mtype.sum$coefficients[1,1]
  phy.mutation.type$Intercept.SE[i]=mtype.sum$coefficients[1,2]
  phy.mutation.type$Intercept.Z[i]=mtype.sum$coefficients[1,3]
  phy.mutation.type$Intercept.P[i]=mtype.sum$coefficients[1,4]
  phy.mutation.type$Nonsynonymous[i]=mtype.sum$coefficients[2,1]
  phy.mutation.type$Nonsynonymous.SE[i]=mtype.sum$coefficients[2,2]
  phy.mutation.type$Nonsynonymous.Z[i]=mtype.sum$coefficients[2,3]
  phy.mutation.type$Nonsynonymous.P[i]=mtype.sum$coefficients[2,4]
  phy.mutation.type$Indel[i]=mtype.sum$coefficients[3,1]
  phy.mutation.type$Indel.SE[i]=mtype.sum$coefficients[3,2]
  phy.mutation.type$Indel.Z[i]=mtype.sum$coefficients[3,3]
  phy.mutation.type$Indel.P[i]=mtype.sum$coefficients[3,4]
  phy.mutation.type$Scale[i]=mtype.sum$scale
}

#Add to the dataframe to get pairwise comparison of nonsynonymous and indels
phy.mutation.type$Non.Indel=0
phy.mutation.type$Non.Indel.SE=0
phy.mutation.type$Non.Indel.Z=0
phy.mutation.type$Non.Indel.P=0
phy.mutation.type$Non.Indel.Scale=0
phy.mutation.type$Non.Intercept=0
phy.mutation.type$Non.Intercept.SE=0
phy.mutation.type$Non.Intercept.Z=0
phy.mutation.type$Non.Intercept.P=0


mam.syn2$Type.2<-"1-Nonsynonymous"
mam.syn2$Type.2[mam.syn2$Type=="1-synonymous"]<-"2-Synonymous"
mam.syn2$Type.2[mam.syn2$Type=="3-indel"]<-"3-Indel"


#Phylogenetic glm
for (i in 1:length(mam.tree4)){
  mutation_type<-phyloglm(Variants~Type.2, phy=mam.tree5[[i]],data=mam.syn2, method="poisson_GEE")
  mtype.2.sum<-summary(mutation_type)
  ###Write results to data frame
  phy.mutation.type$Non.Intercept[i]=mtype.2.sum$coefficients[1,1]
  phy.mutation.type$Non.Intercept.SE[i]=mtype.2.sum$coefficients[1,2]
  phy.mutation.type$Non.Intercept.Z[i]=mtype.2.sum$coefficients[1,3]
  phy.mutation.type$Non.Intercept.P[i]=mtype.2.sum$coefficients[1,4]
  phy.mutation.type$Non.Indel[i]=mtype.2.sum$coefficients[3,1]
  phy.mutation.type$Non.Indel.SE[i]=mtype.2.sum$coefficients[3,2]
  phy.mutation.type$Non.Indel.Z[i]=mtype.2.sum$coefficients[3,3]
  phy.mutation.type$Non.Indel.P[i]=mtype.2.sum$coefficients[3,4]
  phy.mutation.type$Non.Indel.Scale[i]=mtype.2.sum$scale
}

write.csv(phy.mutation.type,"N.genomes_mutation_type.csv")

#Add to the dataframe to get incidence rates of variant types in phylogenetic GLM
phy.mutation.type$IR.Syn=exp(phy.mutation.type$Intercept)
phy.mutation.type$IR.Nonsyn=exp(phy.mutation.type$Intercept+phy.mutation.type$Nonsynonymous)
phy.mutation.type$IR.Indel=exp(phy.mutation.type$Intercept +phy.mutation.type$Indel)

hist(phy.mutation.type$IR.Syn, xlim = c(0,1.5), las = 1, breaks =seq(0,1.5, 0.01))
hist(phy.mutation.type$IR.Nonsyn, xlim = c(0,1.5), las = 1, breaks =seq(0,1.5, 0.01))
hist(phy.mutation.type$IR.Indel, xlim = c(0,1.5), las = 1, breaks =seq(0,1.5, 0.01))

#GLMM
non.phylo.model<-glm(Variants ~ Type + offset(log(N.Genomes))+ offset(log(N.SitesDNASP)), data = mam.syn2, family = "poisson")
summary(non.phylo.model)
IR.CI<-confint(non.phylo.model)

#Incidence rates
IR.Syn<-exp(-6.4149)
IR.Nonsyn<-exp(-6.4149+-2.3688)
IR.Indel<-exp(-6.4149+-3.6297)

#Plotting Incidence rates
IR.plot<-data.frame(IR = c(IR.Syn,IR.Nonsyn,IR.Indel),
                    IR.Low = c(exp(IR.CI[1,1]),exp(IR.CI[1,1]+IR.CI[2,1]),exp(IR.CI[1,1]+IR.CI[3,1])),
                    IR.High = c(exp(IR.CI[1,2]),exp(IR.CI[1,2]+IR.CI[2,2]),exp(IR.CI[1,2]+IR.CI[3,2]))
                    )

plot(x= c(0.25,0.5,0.75), y =IR.plot$IR, xlim = c(0,1), las = 1, ylim = c(0,0.0022))
arrows(x0=c(0.25,0.5,0.75),y0=IR.plot$IR,y1=IR.plot$IR.High,angle = 90)
arrows(x0=c(0.25,0.5,0.75),y0=IR.plot$IR,y1=IR.plot$IR.Low,angle = 90)

#Plotting Incidence rate ratios
IRR<-poissonirr(Variants ~ Type + offset(log(N.Genomes))+ offset(log(N.SitesDNASP)), data = mam.syn2)
IRR.plot<-data.frame(IRR = c(IRR$irr[1,1],IRR$irr[2,1]),IRR.High = c(IRR$irr[1,1]-IRR$irr[1,2]*1.96,IRR$irr[2,1]-IRR$irr[2,2]*1.96),IRR.Low = c(IRR$irr[1,1]+IRR$irr[1,2]*1.96,IRR$irr[2,1]+IRR$irr[2,2]*1.96))
plot(x= c(0.5,1), y =IRR.plot$IRR, xlim = c(0,1.5), ylim = c(0,1), las = 1)
abline(h=1)##Synonymous baseline
arrows(x0=c(0.5,1),y0=IRR.plot$IRR,y1=IRR.plot$IRR.Low,angle = 90)
arrows(x0=c(0.5,1),y0=IRR.plot$IRR,y1=IRR.plot$IRR.High,angle = 90)

##posthoc comparison of nonsynonymous and indels
non.phylo.model2<-glm(Variants ~ Type.2 + offset(log(N.Genomes))+ offset(log(N.SitesDNASP)), data = mam.syn2, family = "poisson")
IRR2<-poissonirr(Variants ~ Type.2 + offset(log(N.Genomes))+ offset(log(N.SitesDNASP)), data = mam.syn2)
IRR2.plot<-data.frame(IRR = c(IRR2$irr[1,1],IRR2$irr[2,1]),IRR2.High = c(IRR2$irr[1,1]-IRR2$irr[1,2]*1.96,IRR2$irr[2,1]-IRR2$irr[2,2]*1.96),IRR2.Low = c(IRR2$irr[1,1]+IRR2$irr[1,2]*1.96,IRR2$irr[2,1]+IRR2$irr[2,2]*1.96))
