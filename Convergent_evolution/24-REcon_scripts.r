
#library(devtools)
#install_github("nclark-lab/RERconverge")
library(RERconverge)

rerpath=find.package('RERconverge')
###continus traits
toytreefile = "subsetMammalGeneTrees.txt"
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

data("logAdultWeightcm")
mamRERw = getAllResiduals(toyTrees,useSpecies=names(logAdultWeightcm),
                          transform = "sqrt", weighted = T, scale = T)

noneutherians <- c("Platypus","Wallaby","Tasmanian_devil","Opossum")
par(mfrow=c(1,2))
avgtree=plotTreeHighlightBranches(toyTrees$masterTree, outgroup=noneutherians,
                                  hlspecies=c("Vole","Squirrel"), hlcols=c("blue","red"),
                                  main="Average tree") #plot average tree
bend3tree=plotTreeHighlightBranches(toyTrees$trees$BEND3, outgroup=noneutherians,hlspecies=c("Vole","Squirrel"), hlcols=c("blue","red"),
                                    main="BEND3 tree")
par(mfrow=c(1,1))
phenvExample <- foreground2Paths(c("Vole","Squirrel"),toyTrees,clade="terminal")
plotRers(mamRERw,"BEND3",phenv=phenvExample) #plot RERs

#plot RERs as tree
par(mfrow=c(1,1))
bend3rers = returnRersAsTree(toyTrees, mamRERw, "BEND3", plot = TRUE,
                             phenv=phenvExample) #plot RERs
##binary traits

marineb=read.tree(paste(rerpath,"/extdata/MarineTreeBinCommonNames_noCGM.txt",sep=""))
marinebrooted = root(marineb,outgroup=noneutherians)
plot(marinebrooted)
mb1 = marineb
mb1$edge.length = c(rep(1,length(mb1$edge.length)))
par(mfrow=c(1,2))
plot(marinebrooted, main="Trait tree from file (1)")
#alternative way of representing the tree
binplot1=plotTreeHighlightBranches(mb1, outgroup=noneutherians,
                                   hlspecies=which(marineb$edge.length==1), hlcols="blue",
                                   main="Foreground branches highlighted (1)")
marineextantforeground = c("Walrus","Seal","Killer_whale","Dolphin","Manatee")
marineb2a = foreground2Tree(marineextantforeground, toyTrees, clade="ancestral",
                            useSpecies=names(logAdultWeightcm))
phenvMarine=tree2Paths(marineb, toyTrees)
phenvMarine2=foreground2Paths(marineextantforeground, toyTrees, clade="all")
phenvMarine2b=tree2Paths(marineb2b, toyTrees,corMarine=correlateWithBinaryPhenotype(mamRERw, phenvMarine2, min.sp=10, min.pos=2,
                                      weighted="auto"))

####mydata
#mytre<- read.tree("mastertree.tre")
#mytre = root(mytre,outgroup="Osativa")
#plot(mytre, main="Trait tree from file")
genetrees=readTrees("allaa_pep3.tre")
RERw = getAllResiduals(genetrees,transform = "sqrt", weighted = T, scale = T)
#RERw2 = getAllResiduals(genetrees,transform = "none", weighted = T, scale = T)
#RERw3 = getAllResiduals(genetrees,transform = "sqrt", weighted = T, scale = F)
foreground= c("Crucihimalaya_himalaica","Eutrema_heterophyllum",
              "Salix_brachista","Prunus_mira","Saussurea_obvallata",
              "Rheum_alexandrae","Hordeum_vulgare")
phenv=foreground2Paths(foreground, genetrees, clade="all")
#phenv=tree2Paths(mytre, genetrees)

corMy=correlateWithBinaryPhenotype(RERw, phenv, min.sp=10, min.pos=2,
                                       weighted="auto")
write.csv(corMy,"results_20pep3_gctrees_5.5.csv")
#plot RERs
par(mfrow=c(1,1))
#phenvExample <- foreground2Paths(c("Vole","Squirrel"),toyTrees,clade="terminal")
plotRers(RERw,"OG0000765",phenv=phenv, sortrers = F) #plot RERs
#plot RERs as tree
par(mfrow=c(1,1))
retree = returnRersAsTree(genetrees, RERw, "OG0000765", plot = TRUE,
                             phenv=phenv) #plot RERs

newretree = treePlotRers(treesObj=genetrees, rermat=RERw, index="OG0003553",
                            type="c", nlevels=9, figwid=10)
par(mfrow=c(1,1))
hist(corMy$P)


###  plot 
plotRers(RERw,"OG0000000",
         phenv=phenv,
        # xlims = c(-3,3)
         )

treePlotRers(treesObj=genetrees, rermat=RERw, index="OG0000765",
             type="c",phenv=phenv, nlevels=9)
