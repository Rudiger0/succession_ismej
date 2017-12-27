####SCRIPT Ortiz-Álvarez 2017 Succession article####
setwd("")

library(ggplot2)
library(plyr)
library(reshape2)
library(vegan)
library(GUniFrac)
#mutiplot function: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

#########1. DATA FILES########
#genome_readtable
genometable <- read.table("succession_genometable.txt", sep="\t", header=T, row.names=1)
#metadata
meta <- read.table("meta_final.txt", sep="\t", header=T, check.names = F)
names(meta)
meta <- meta[,c(6,8,10:12,14)]
row.names(meta) <- meta$Sample_code
meta <- subset(meta, meta$Code != "A10") #optional: remove A10 from the analysis
#functional table
traits <- read.table(file = "succession_traits.txt", sep="\t", row.names=1, header=T, check.names = F)
#Seq count of sequences to match
seqs_tomatch <- read.table("seqcount_tomatch.txt", sep="\t", header=T, check.names=F)

####1b. Overview####
#Filter and quick overview of sequences and matches
genometable <- genometable[row.names(meta),]
melt(as.list(sort(rowSums(genometable)))) #counts reads per sample (no rarefied)
genometable.pa <- decostand(genometable, method="pa")
melt(as.list(sort(rowSums(genometable.pa)))) #counts nr matches per sample (no rarefied)

#Quick view of early-late match proportion
matches <- merge(as.data.frame(melt(as.list(sort(rowSums(genometable))))), meta, by.x="L1", by.y="Sample_code")
matches <- merge(matches, seqs_tomatch, by.x="L1", by.y="Sample_code")
matches <- mutate(matches, seq_proportion=value/Sequences*100)

ggplot(data=matches, aes(x = State, y=seq_proportion))+
  geom_jitter(alpha=0.2, size=0.5)+
  geom_boxplot(fill="grey80", outlier.shape = NULL, outlier.size = 1, alpha=0.5)+
  facet_wrap(~Habitat, scales="free", ncol=4)+
  theme_bw()

#Average % of matches per habitat
sum(matches$value)/sum(matches$Sequences)*100
a <- aggregate(matches$value~matches$Habitat, FUN=sum)
b <- aggregate(matches$Sequences~matches$Habitat, FUN=sum)
c <- a[,2]/b[,2]*100
names(c) <- (a[,1])

aggregate(matches$seq_proportion~matches$Habitat, FUN=median) #+matches$State

##Shared OTUs between habitat pairs
coral.hab <- merge(genometable, meta, by.x="row.names", by.y="Sample_code")
coral.hab <- split(coral.hab, coral.hab$State)
coral.hab.e <- coral.hab$E
coral.hab.l <- coral.hab$L

#enter early or late separatedly
coral.hab <- ddply(.data = coral.hab.e[,c(2:1845)], .(coral.hab.e$Habitat), colMeans) 
row.names(coral.hab) <- coral.hab[,1]
coral.hab[,1] <- NULL
coral.hab <- decostand(coral.hab, method="pa")
hab.jaccard <- as.data.frame(as.matrix((1-vegdist(coral.hab, method="jaccard"))))
hab.jaccard[hab.jaccard == 0] <- NA

mean(colMeans(hab.jaccard[,c(1:8)], na.rm=TRUE))
min(hab.jaccard[,c(1:8)], na.rm=TRUE)
max(hab.jaccard[,c(1:8)], na.rm=TRUE)

#To relative percentages
relative <- apply(t(genometable), 2, function(x){x/sum(x)})
colSums(relative)

#Prepare weighted.ko matrix using the relative matrix
relative <- as.data.frame(relative)
wtd.sim <- merge(relative, t(komerge), by="row.names") 
relative.sim <- wtd.sim[,2:(ncol(relative)+1)] #for samples
ko.sim <- wtd.sim[,(ncol(relative)+2):ncol(wtd.sim)] #for functions
simulated.ko <- apply(relative.sim, 2, function(x){apply(ko.sim,2,function(y){sum(x*y)})})

#equalize kegg profile per sample
simulated.ko <- decostand(t(simulated.ko), method="total")
rowSums(simulated.ko)
simulated.ko <- t(simulated.ko)
simulated.pa <- decostand(simulated.ko, method="pa")

#weighted genome ubiquity: taxonomy#
tax.ocu <- as.data.frame(colSums(genometable.pa))
tax.ocu <- as.data.frame(apply(relative, 2, function(x){sum(x*tax.ocu)}))
tax.ocu <- merge(tax.ocu, meta, by="row.names")
names(tax.ocu)[2] <- "weighted_taxocu"
aggregate(tax.ocu$weighted_taxocu~tax.ocu$State, FUN=mean)
tax.ocu.code <- aggregate(tax.ocu$weighted_taxocu~tax.ocu$State+tax.ocu$Code, FUN=mean)
names(tax.ocu.code) <- c("State", "Code", "weighted_taxocu")
tax.ocu.code <- split(tax.ocu.code, tax.ocu.code$State)
tax.ocu.change <- as.data.frame(tax.ocu.code$L$weighted_taxocu-tax.ocu.code$E$weighted_taxocu)
row.names(tax.ocu.change) <- tax.ocu.code$E$Code
tax.ocu.change <- cbind.data.frame(tax.ocu.change, as.data.frame(tax.ocu.change>0), tax.ocu.code$E$Code)
names(tax.ocu.change) <- c("tax.ocu.change", "magnitude", "Code")
tax.ocu.change <- unique(merge(tax.ocu.change, meta[,c(1, 3:4)], by="Code"))

#test for significant differences considering habitat effect: taxonomy
shapiro.test(tax.ocu$weighted_taxocu) #not normal
bartlett.test(tax.ocu$weighted_taxocu~tax.ocu$State) #not homocedastic
adonis(tax.ocu$weighted_taxocu~tax.ocu$Habitat+tax.ocu$State, permutations=999, method="euclidean")
adonis(tax.ocu$weighted_taxocu~tax.ocu$Code+tax.ocu$State, permutations=999, method="euclidean")

#weighted genome ubiquity: functions#
fun.ocu <- as.data.frame(colSums(simulated.pa))
fun.ocu <- as.data.frame(apply(relative, 2, function(x){sum(x*fun.ocu)}))
fun.ocu <- merge(fun.ocu, meta, by="row.names")
names(fun.ocu)[2] <- "weighted_funocu"
aggregate(fun.ocu$weighted_funocu~fun.ocu$State, FUN=mean)
fun.ocu.code <- aggregate(fun.ocu$weighted_funocu~fun.ocu$State+fun.ocu$Code, FUN=mean)
names(fun.ocu.code) <- c("State", "Code", "weighted_funocu")
fun.ocu.code <- split(fun.ocu.code, fun.ocu.code$State)
fun.ocu.change <- as.data.frame(fun.ocu.code$L$weighted_funocu-fun.ocu.code$E$weighted_funocu)
row.names(fun.ocu.change) <- fun.ocu.code$E$Code
fun.ocu.change <- cbind.data.frame(fun.ocu.change, as.data.frame(fun.ocu.change>0), fun.ocu.code$E$Code)
names(fun.ocu.change) <- c("fun.ocu.change", "magnitude", "Code")
fun.ocu.change <- unique(merge(fun.ocu.change, meta[,c(1, 3:4)], by="Code"))

#test for significant differences considering habitat effect: functions
shapiro.test(fun.ocu$weighted_funocu) #not normal
bartlett.test(fun.ocu$weighted_funocu~fun.ocu$State) #not homocedastic
adonis(fun.ocu$weighted_funocu~fun.ocu$Habitat*fun.ocu$State, permutations=999, method="euclidean")
adonis(fun.ocu$weighted_funocu~fun.ocu$Code*fun.ocu$State, permutations=999, method="euclidean")

#Plots
fig.ubiquity.sup <- ggplot(tax.ocu.change, aes(Code, tax.ocu.change))+
  geom_bar(aes(colour=magnitude, fill=magnitude),alpha=0.8, stat="identity")+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+#, breaks=c("TRUE", "FALSE"), labels=c("L-E>0", "L-E<0"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+#, breaks=c("TRUE", "FALSE"), labels=c("L-E>0", "L-E<0"))+
  geom_hline(aes(yintercept=mean(functional.beta, na.rm=T)), linetype=2)+
  #ylab(NULL)+
  xlab("Weighted occurrence change")+
  facet_wrap(~For_figures, scales="free", nrow=1)+
  theme_bw()#+

fig.ubiquity.fun <- ggplot(fun.ocu.change, aes(Code, fun.ocu.change))+
  geom_bar(aes(colour=magnitude, fill=magnitude),alpha=0.8, stat="identity")+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+#, breaks=c("TRUE", "FALSE"), labels=c("L-E>0", "L-E<0"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+#, breaks=c("TRUE", "FALSE"), labels=c("L-E>0", "L-E<0"))+
  geom_hline(aes(yintercept=mean(functional.beta, na.rm=T)), linetype=2)+
  xlab("Weighted occurrence change")+
  facet_wrap(~For_figures, scales="free", nrow=1)+
  theme_bw()#+

ggsave(filename = "Fig_supple_taxubicuity.pdf", plot= fig.ubiquity.sup, width = 25, height = 6.25, units = "cm")
ggsave(filename = "Fig_supple_funubicuity.pdf", plot= fig.ubiquity.fun, width = 25, height = 6.25, units = "cm")

######2. Dissimilarities For functions and structures######  
#taxonomic and functional matrices transformation
tax.hel<-decostand(t(relative), method="hellinger")
fun.hel<-decostand(t(simulated.ko), method="hellinger")
tax.dist <-vegdist(tax.hel, method="bray")
fun.dist <-vegdist(fun.hel, method="bray")

#nMDS (all samples together)
tax.nmds<-metaMDS(tax.dist)
fun.nmds<-metaMDS(fun.dist)
tax.nmds<-data.frame(tax.nmds$points)
fun.nmds<-data.frame(fun.nmds$points)

##Separated matrices for E/L: Taxonomy##
tax.hel.merge <- merge(tax.hel, meta, by="row.names")
tax.hel.merge <- split(tax.hel.merge, tax.hel.merge$State)
tax.early <- tax.hel.merge$E
tax.late <- tax.hel.merge$L
row.names(tax.early) <- tax.early$Row.names
row.names(tax.late) <- tax.late$Row.names

#Calculate distance matrices
early.dist<-vegdist(tax.early[,2:(ncol(tax.hel)+1)], method="bray")
late.dist<-vegdist(tax.late[,2:(ncol(tax.hel)+1)], method="bray")

#dist matrices by early or late state, for independant nmds's.
tax.nmds<-metaMDS(early.dist)
succ.nmds.points<-data.frame(tax.nmds$points)
succ.nmds.points.e <- merge(succ.nmds.points, meta, by="row.names")

tax.nmds<-metaMDS(late.dist)
succ.nmds.points <-data.frame(tax.nmds$points)
succ.nmds.points.l <- merge(succ.nmds.points, meta, by="row.names")

find_hull <- function(df) df[chull(df$MDS1, df$MDS2), ]
hulls.e <- ddply(succ.nmds.points.e, "Habitat", find_hull)
hulls.l <- ddply(succ.nmds.points.l, "Habitat", find_hull)

seg.e <- merge(succ.nmds.points.e,aggregate(cbind(mean.MDS1=MDS1,mean.MDS2=MDS2)~Habitat,succ.nmds.points.e,mean),by="Habitat")
seg.l <- merge(succ.nmds.points.l,aggregate(cbind(mean.MDS1=MDS1,mean.MDS2=MDS2)~Habitat,succ.nmds.points.l,mean),by="Habitat")

#Early
qpe <- ggplot(seg.e, aes(MDS1, MDS2))+
  geom_polygon(data=hulls.e, aes(fill=Habitat), alpha=0.5)+
  #geom_point(data=succ.nmds.points.e, aes(colour=Specific, shape=State), size=5, alpha=0.7)+
  #geom_text(data=seg.e,aes(label=Code), position=position_jitter(h=0.01,w=0.01), size=3, color="black")+
  geom_text(data=seg.e,aes(label=Code), size=3, color="black")+
  #facet_wrap(~State)
  scale_fill_manual(values = c("brown1","yellow", "darkolivegreen4", "darkseagreen3", "blueviolet", "dodgerblue", "lightseagreen", "sandybrown"))+
  scale_linetype_discrete() +
  theme(legend.position="none")+
  geom_segment(aes(x=mean.MDS1,y=mean.MDS2, xend=MDS1, yend=MDS2), linetype="dotted")+
  xlim(c(-0.50, 0.60))+ #structural
  ylim(c(-0.4, 0.45))+ #structural
  #xlim(c(-5, 10))+ #functional ko
  #ylim(c(-6, 4))+ #functional ko
  theme_bw()

#Late
qpl <- ggplot(seg.l, aes(MDS1, MDS2))+
  geom_polygon(data=hulls.l, aes(fill=Habitat), alpha=0.5)+
  #geom_point(data=succ.nmds.points.l, aes(colour=Specific, shape=State), size=5, alpha=0.7)+
  #geom_text(data=seg.l,aes(label=Code), position=position_jitter(h=0.01,w=0.01), size=3, color="black")+
  geom_text(data=seg.l,aes(label=Code), size=3, color="black")+
  #facet_wrap(~State)
  theme(legend.position="none")+
  scale_fill_manual(values = c("brown1","yellow", "darkolivegreen4", "darkseagreen3", "blueviolet", "dodgerblue", "lightseagreen", "sandybrown"))+
  geom_segment(aes(x=mean.MDS1,y=mean.MDS2, xend=MDS1, yend=MDS2), linetype="dotted")+
  xlim(c(-0.50, 0.55))+ #structural
  ylim(c(-0.4, 0.40))+ #structural
  #xlim(c(-5, 10))+ #functional ko
  #ylim(c(-6, 4))+ #functional ko
  theme_bw()

qpe
qpl

####2b. ANOSIM and distance matrices by state/habitat####
#Anosim (habitat overlap evaluation)
anosim(early.dist, tax.early$Habitat, permutations=999)
anosim(late.dist, tax.late$Habitat, permutations=999)

#It is possible to examine habitats separatedly too.
#Anosim: only between human and primate
tax.hel.merge <- merge(tax.hel, meta, by="row.names")
tax.hel.merge <- split(tax.hel.merge, tax.hel.merge$For_figures)

gut.otus <- tax.hel.merge$"Infant Gut microbiome"
row.names(gut.otus) <- gut.otus$Row.names
gut.otus$Row.names <- NULL
gut.otus <- split(gut.otus, gut.otus$State)
e.gut.otus <- gut.otus$E[,c(1:1844, 1848)]
e.gut.dist <- vegdist(e.gut.otus[,1:1844], method="bray")
l.gut.otus <- gut.otus$L[,c(1:1844, 1848)]
l.gut.dist <- vegdist(l.gut.otus[,1:1844], method="bray")

anosim(e.gut.dist, e.gut.otus$Habitat, permutations=999)
anosim(l.gut.dist, l.gut.otus$Habitat, permutations=999)

#####2c. Analysis of habitat and Stage in sample dissimilarities####
labels(tax.dist) == row.names(meta)
adonis(tax.dist~meta$Habitat*meta$State)

labels(fun.dist) == row.names(meta)
adonis(fun.dist~meta$Habitat*meta$State)

#####2d. BETA DIVERSITY DIFFERENCES (based on bray curtis dissimilarities)####
#Taxonomy
d <- as.matrix(tax.dist)
oldnames <- row.names(d)
row.names(d) <- c(1:117)
colnames(d) <- c(1:117)

df <- melt(d, varnames = c("row", "col"))
df[df$row > df$col,]
df[as.numeric(df$row) > as.numeric(df$col), ]

melted <- (df[as.numeric(df$row) > as.numeric(df$col), ])

oldnames <- as.data.frame(oldnames)
oldnames <- cbind(oldnames, c(1:117))
names(oldnames) <- c("oldnames", "newnames")
oldnames <- merge(oldnames, meta, by.x="oldnames", by.y="Sample_code")

melted.1 <- merge(melted, oldnames, by.x="row", by.y="newnames")
melted.2 <- merge(melted.1, oldnames, by.x="col", by.y="newnames")
#the melted.2 table, has the first element of the comparison named as .x, and
#the second named as .y. Therefore, now we can subset by E/L and then mean by the factor
#we want to: Code, Specific and Habitat

#Subset by successional stage (out are those between states)
melted.2.early <- subset(melted.2, State.x=="E" & State.y=="E")
melted.2.late <- subset(melted.2, State.x=="L" & State.y=="L")

melted.2.early.code <- subset(melted.2.early, Code.x==Code.y)
melted.2.late.code <- subset(melted.2.late, Code.x==Code.y)

melted.2.early.Specific <- subset(melted.2.early, Specific.x==Specific.y)
melted.2.late.Specific <- subset(melted.2.late, Specific.x==Specific.y)

melted.2.early.Habitat <- subset(melted.2.early, Habitat.x==Habitat.y)
melted.2.late.Habitat <- subset(melted.2.late, Habitat.x==Habitat.y)

#Mean by Code
code.means.early <- aggregate(melted.2.early.code$value~melted.2.early.code$Code.x, FUN=mean)
code.means.late <- aggregate(melted.2.late.code$value~melted.2.late.code$Code.x, FUN=mean)
code.means.late[,1]==code.means.early[,1]
code.mean <- code.means.late[,2]-code.means.early[,2]
names(code.mean) <- code.means.early[,1]
structural.beta <- code.mean

#Mean by Habitat
Habitat.means.early <- ddply(.data = melted.2.early.Habitat[,c(2:3)], .(melted.2.early.Habitat$Habitat.x), colMeans)
Habitat.means.late <- ddply(.data = melted.2.late.Habitat[,c(2:3)], .(melted.2.late.Habitat$Habitat.x), colMeans)
Habitat.means.late[,1]==Habitat.means.early[,1]
Habitat.mean <- (Habitat.means.late[,3]-Habitat.means.early[,3])/Habitat.means.late[,3]
names(Habitat.mean) <- Habitat.means.late$"melted.2.late.Habitat$Habitat.x"
Habitat.mean

#Functional
d <- as.matrix(fun.dist)
oldnames <- row.names(d)
row.names(d) <- c(1:117)
colnames(d) <- c(1:117)

df <- melt(d, varnames = c("row", "col"))
df[df$row > df$col,]
df[as.numeric(df$row) > as.numeric(df$col), ]

melted <- (df[as.numeric(df$row) > as.numeric(df$col), ])

oldnames <- as.data.frame(oldnames)
oldnames <- cbind(oldnames, c(1:117))
names(oldnames) <- c("oldnames", "newnames")
oldnames <- merge(oldnames, meta, by.x="oldnames", by.y="Sample_code")

melted.1 <- merge(melted, oldnames, by.x="row", by.y="newnames")
melted.2 <- merge(melted.1, oldnames, by.x="col", by.y="newnames")

#Subset by successional stage (out are those between states)
melted.2.early <- subset(melted.2, State.x=="E" & State.y=="E")
melted.2.late <- subset(melted.2, State.x=="L" & State.y=="L")

melted.2.early.code <- subset(melted.2.early, Code.x==Code.y)
melted.2.late.code <- subset(melted.2.late, Code.x==Code.y)

melted.2.early.Specific <- subset(melted.2.early, Specific.x==Specific.y)
melted.2.late.Specific <- subset(melted.2.late, Specific.x==Specific.y)

melted.2.early.Habitat <- subset(melted.2.early, Habitat.x==Habitat.y)
melted.2.late.Habitat <- subset(melted.2.late, Habitat.x==Habitat.y)

#Mean by Code
code.means.early <- aggregate(melted.2.early.code$value~melted.2.early.code$Code.x, FUN=mean)
code.means.late <- aggregate(melted.2.late.code$value~melted.2.late.code$Code.x, FUN=mean)
code.means.late[,1]==code.means.early[,1]
code.mean <- code.means.late[,2]-code.means.early[,2]
names(code.mean) <- code.means.early[,1]
functional.beta <- code.mean

#Mean by Habitat
Habitat.means.early <- ddply(.data = melted.2.early.Habitat[,c(2:3)], .(melted.2.early.Habitat$Habitat.x), colMeans)
Habitat.means.late <- ddply(.data = melted.2.late.Habitat[,c(2:3)], .(melted.2.late.Habitat$Habitat.x), colMeans)
Habitat.means.late[,1]==Habitat.means.early[,1]
Habitat.mean <- (Habitat.means.late[,3]-Habitat.means.early[,3])/Habitat.means.late[,3]
names(Habitat.mean) <- Habitat.means.late$"melted.2.late.Habitat$Habitat.x"
Habitat.mean

betadiversity <- cbind.data.frame(structural.beta, functional.beta)

########3. ALPHA DIVERSITY CHANGES#######
########3a.Shannon diversity for functional samples and FUNquad######
sim.shan <- split(merge(t(simulated.ko), meta, by="row.names"), merge(t(simulated.ko), meta, by="row.names")$State)
col.sim.shan <- lapply(sim.shan, function(x){ddply(.data=x[,2:8192], .(x$Code), colMeans)})
sinchan <- lapply(col.sim.shan, function(x){diversity(x[,2:8192], index="shannon")})

funshandiff <- as.data.frame(sinchan$L - sinchan$E)
row.names(funshandiff) <- col.sim.shan$E[,1]
functional.alpha <- funshandiff

########3b.Structural shannon#####
#Code using the average of 100 rarefactions directly in R
res <- lapply(as.list(1:100), function(x) Rarefy(genometable, 50)$otu.tab.rff)
res.shannon <- lapply(res, function(x) {diversity(as.data.frame(x), index="shannon")})
res.shannon.1 <- as.data.frame(rowMeans(as.data.frame(res.shannon)))
names(res.shannon.1) <- "ave.shannon"
index.meta <- merge(res.shannon.1, meta, by="row.names")

structural.alpha <- split(index.meta, index.meta$State)
structural.alpha.1 <- aggregate(as.data.frame(structural.alpha$E)$ave.shannon~as.data.frame(structural.alpha$E)$Code, FUN=mean)
structural.alpha.2 <- aggregate(as.data.frame(structural.alpha$L)$ave.shannon~as.data.frame(structural.alpha$L)$Code, FUN=mean)
#structural.alpha <- as.data.frame(structural.alpha.2[,2] - structural.alpha.1[,2])
names(structural.alpha.1) <- c("Code", "ave.shannon.E")
names(structural.alpha.2) <- c("Code", "ave.shannon.L")
structural.alpha <- merge(structural.alpha.1, structural.alpha.2, by="Code")
#row.names(structural.alpha) <- structural.alpha.2[,1]
structural.alpha <- mutate(structural.alpha, structural.alpha=ave.shannon.L-ave.shannon.E)
names(structural.alpha) <- c("Code", "ave.shannon.E", "ave.shannon.L", "structural.alpha")
structural.alpha.prior <- unique(merge(structural.alpha, meta, by="Code"))

ggplot(structural.alpha.prior, aes(Code, structural.alpha))+
  geom_bar(alpha=0.8, stat="identity")+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black", guide=FALSE))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70", guide=FALSE))+
  geom_hline(aes(yintercept=mean(structural.alpha)), linetype=2)+
  #xlab(NULL)+
  ylab("Structural (Genomic-OTUs)")+
  #scale_y_continuous(labels = fmt)+
  #geom_text(aes(label=as.factor(Row.names)))+
  theme_bw()


#########4.Table with the 4 diversity metrics#######
beta.metrics <- cbind.data.frame(structural.beta, functional.beta)
alpha.metrics <- cbind.data.frame(structural.alpha, functional.alpha)
diversity.metrics <- merge(beta.metrics, alpha.metrics, by="row.names", all=T)[,c(2:4, 7:8)]
names(diversity.metrics) <- c("structural.beta", "functional.beta", "Code", "structural.alpha", "functional.alpha")
row.names(diversity.metrics) <- diversity.metrics[,3]
diversity.metrics[,3] <- NULL

########4b. Barplot per thinggie########
diversity.meta <- unique(merge(diversity.metrics, meta[,c(1, 3:4)], by.x="row.names", by.y="Code"))
diversity.cats <- as.data.frame(diversity.metrics>0)
names(diversity.cats) <- c("sb", "fb", "sa", "fa")
diversity.meta <- merge(diversity.meta, diversity.cats, by.x="Row.names", by.y="row.names")

fmt<-function(x){format(x,nsmall = 1,scientific = FALSE)}

sa <- ggplot(diversity.meta, aes(as.factor(Row.names), structural.alpha))+
  geom_bar(aes(colour=sa, fill=sa), alpha=0.8, stat="identity")+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black", guide=FALSE))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70", guide=FALSE))+
  geom_hline(aes(yintercept=mean(structural.alpha)), linetype=2)+
  ylab("Structural (Genomic-OTUs)")+
  scale_y_continuous(labels = fmt)+
  theme_classic()

sb <- ggplot(diversity.meta, aes(as.factor(Row.names), structural.beta))+
  geom_bar(aes(colour=sb, fill=sb), alpha=0.8,stat="identity")+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+#, breaks=c("TRUE", "FALSE"), labels=c("L-E>0", "L-E<0"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+#, breaks=c("TRUE", "FALSE"), labels=c("L-E>0", "L-E<0"))+
  geom_hline(aes(yintercept=mean(structural.beta, na.rm=T)), linetype=2)+
  scale_y_continuous(breaks=c(-0.35, -0.15, 0, 0.15))+
  theme_classic()#+

fa <- ggplot(diversity.meta, aes(as.factor(Row.names), functional.alpha))+
  geom_bar(aes(colour=fa, fill=fa), alpha=0.8,stat="identity")+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black", guide=FALSE))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70", guide=FALSE))+
  geom_hline(aes(yintercept=mean(functional.alpha, na.rm=T)), linetype=2)+
  ylab("Functional (Genomic-OTUs)")+
  xlab("Shannon diversity changes")+
  theme_classic()#+

fb <- ggplot(diversity.meta, aes(as.factor(Row.names), functional.beta))+
  geom_bar(aes(colour=fb, fill=fb),alpha=0.8, stat="identity")+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+#, breaks=c("TRUE", "FALSE"), labels=c("L-E>0", "L-E<0"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+#, breaks=c("TRUE", "FALSE"), labels=c("L-E>0", "L-E<0"))+
  geom_hline(aes(yintercept=mean(functional.beta, na.rm=T)), linetype=2)+
  scale_y_continuous(breaks=c(-0.2, -0.1, 0, 0.1), limits=c(-0.2, 0.1))+
  xlab("Betadispersion changes (BC)")+
  theme_classic()#+

multiplot(sa, fa, sb, fb, cols=2)

##########5.Weighted trait information##########
relative.f <- relative[rownames(traits),]
row.names(traits)==row.names(relative.f) 
wtd.traits <- apply(relative.f, 2, function(x){apply(traits,2,function(y){sum(x*y)})})
row.names(wtd.traits) <- c("genome.size", "gene.count", "GC", "CDSp", "RNAp", "16Scount", "RNAnum", "rRNAnum", "tRNAnum")

########5a. Metabolism########

##WEIGHTED METABOLIC KEGGS (ONLY) PER SAMPLE
metab_kos <- read.csv("metabolism_ko.txt", sep="\t", header=T, check.names = F, row.names=1)
metab.filter <-komerge[rownames(metab_kos),]
metab.filter.t <- as.data.frame(t(metab.filter))
metab.filter <-metab.filter.t[rownames(relative),]

wtd.sim.prep <- merge(relative, metab.filter, by="row.names") 
relative.metab <- wtd.sim.prep[,2:118] #para metab.filter
metab.filter <- wtd.sim.prep[,119:195]
simulated <- apply(relative.metab, 2, function(x){apply(metab.filter,2,function(y){sum(x*y)})})

#Prepare weighted metabolism table for downstream analyses
row.names(simulated) == row.names(metab_kos)
simulated <- as.data.frame(simulated)
metab.hel <- ddply(.data = simulated[,c(1:117)], .(metab_kos$Metabolismo), colMeans) #Group o Metabolismo
row.names(metab.hel) <- metab.hel$"metab_kos$Metabolismo"
metab.hel$"metab_kos$Metabolismo" <- NULL
#metab.hel<-as.data.frame(t(decostand(metab.hel, method="hellinger"))) #scale by sample, not by metabolism.

metab.and.traits <- merge(t(wtd.traits), t(metab.hel), by="row.names") #ojo el transponer
metab.and.traits.meta <- merge(metab.and.traits, meta, by.x="Row.names", by.y="row.names")

########5b. rRNAnum & Autotrophy aggregate and changes. Plus Genomesize and GC ######

metab.and.traits.split <- split(metab.and.traits.meta, metab.and.traits.meta$State)

traits.early <- ddply(as.data.frame(metab.and.traits.split$E)[,c(2:37)], .(as.data.frame(metab.and.traits.split$E)$Code), colMeans)
traits.late <- ddply(as.data.frame(metab.and.traits.split$L)[,c(2:37)], .(as.data.frame(metab.and.traits.split$L)$Code), colMeans)

traits.changes <- as.data.frame(traits.late[,c(2:37)]-traits.early[,c(2:37)])
row.names(traits.changes) <- traits.late[,1]

auto.early <- aggregate(as.data.frame(metab.and.traits.split$E)$Arnon_C_Fix+as.data.frame(metab.and.traits.split$E)$Aer_C_Fix+as.data.frame(metab.and.traits.split$E)$Wood_C_Fix~as.data.frame(metab.and.traits.split$E)$Code, FUN=mean)
auto.late <- aggregate(as.data.frame(metab.and.traits.split$L)$Arnon_C_Fix+as.data.frame(metab.and.traits.split$L)$Aer_C_Fix+as.data.frame(metab.and.traits.split$L)$Wood_C_Fix~as.data.frame(metab.and.traits.split$L)$Code, FUN=mean)

traits.changes <- cbind.data.frame(traits.changes, as.data.frame(auto.late[,2]-auto.early[,2]), as.data.frame(auto.early[,1]))
names(traits.changes)[c(37,38)] <- c("auto.changes", "Code")

traits.cat <- as.data.frame(traits.changes[,c(1,3,8,23,28,37)]>0)
names(traits.cat) <- c("gsizec", "gcc", "rc", "nfix", "phosc", "ac")
traits.changes <- cbind.data.frame(traits.changes, traits.cat)

fmt<-function(x){format(x,nsmall = 3,scientific = TRUE)}

rRNA.bar <- ggplot(traits.changes, aes(as.factor(Code), rRNAnum))+
  geom_bar(aes(colour=rc, fill=rc), alpha=0.8, stat="identity")+
  scale_y_continuous(labels = fmt)+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+
  geom_hline(aes(yintercept=mean(rRNAnum)), linetype=2)+
  theme_bw()#+

auto.bar <- ggplot(traits.changes, aes(as.factor(Code), auto.changes))+
  geom_bar(aes(colour=ac, fill=ac), alpha=0.8, stat="identity")+
  scale_y_continuous(labels = fmt)+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+
  geom_hline(aes(yintercept=mean(auto.changes)), linetype=2)+
  xlab("Sample code")+
  theme(panel.background = element_rect(fill = "white", colour="black"), panel.grid.minor = element_line(colour = NULL), panel.grid.major = element_line(colour = NULL))

gsize.bar <- ggplot(traits.changes, aes(as.factor(Code), genome.size))+
  geom_bar(aes(colour=gsizec, fill=gsizec), alpha=0.8,stat="identity")+
  ylab(label="Genome Size")+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+
  geom_hline(aes(yintercept=mean(genome.size)), linetype=2)+
  theme(panel.background = element_rect(fill = "white", colour="black"), panel.grid.minor = element_line(colour = NULL), panel.grid.major = element_line(colour = NULL))

GC.bar <- ggplot(traits.changes, aes(as.factor(Code), GC))+
  geom_bar(aes(colour=gcc, fill=gcc), alpha=0.8, stat="identity")+
  scale_y_continuous(labels = fmt)+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+
  geom_hline(aes(yintercept=mean(GC)), linetype=2)+
  theme(panel.background = element_rect(fill = "white", colour="black"), panel.grid.minor = element_line(colour = NULL), panel.grid.major = element_line(colour = NULL))

phosphate <- ggplot(traits.changes, aes(as.factor(Code), P_transport_high))+
  geom_bar(aes(colour=phosc, fill=phosc), alpha=0.8, stat="identity")+
  scale_y_continuous(labels = fmt)+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+
  geom_hline(aes(yintercept=mean(P_transport_high)), linetype=2)+
  xlab("Sample code")+
  theme_bw()

nitrogen <- ggplot(traits.changes, aes(as.factor(Code), N_Fix))+
  geom_bar(aes(colour=nfix, fill=nfix), alpha=0.8, stat="identity")+
  scale_y_continuous(labels = fmt, limits = c(-0.3, 0.06))+
  scale_colour_manual(values=c("FALSE" ="black", "TRUE" = "black"))+
  scale_fill_manual(values=c("FALSE" ="gray20", "TRUE" = "gray70"))+
  geom_hline(aes(yintercept=mean(N_Fix)), linetype=2)+
  xlab("Sample code")+
  theme(panel.background = element_rect(fill = "white", colour="black"), panel.grid.minor = element_line(colour = NULL), panel.grid.major = element_line(colour = NULL))
