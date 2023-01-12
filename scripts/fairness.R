
##########################################
##  fairness: R script for Fairness paper
##  Note_1: Code to reproduce the analyses in 
##    McDonnell, Stoltz, and Taylor (2022) in SER.
##  Note_2: The data necessary to replicate the analyses
##    are not publicly available. Please contact Erin
##    McDonnell (erin.mcdonnell@nd.edu) to inquire about 
##    access to the data.
##  Author: Marshall A. Taylor
##########################################

### BEGIN ###

##########################################

# Install necessary packages
#install.packages("pacman")
pacman::p_load(igraph,qgraph,ggplot2,grid,
               gridExtra,gridBase,plotrix,poLCA,pheatmap,
               imager,scales,psychometric,cluster,ggpubr,
               reshape2,Hmisc,corrplot,factoextra,psych,
               ggcorrplot,ggeffects,devtools,dplyr,lavaan,
               install = T)

# devtools::install_version("corclass", version = "0.1.1",
#                          repos = "http://cran.us.r-project.org")
library(corclass)

source("scripts/cca2.R")
source("scripts/clusGap2.R")

##########################################

##########################################

# Load in and prep data
data <- read.csv("data/fairness.csv", header = T) #NOTE: Contact Erin McDonnell
                                              #(erin.mcdonnell@nd.edu) to inquire 
                                              #about access to the data

data <- data[which(data$wave == 4),] #Only use wave 4

vars <- c("painter","autoa","redwageb","fixedinc","pb","grocery",
          "apples","dolls","hiring")
data <- data[complete.cases(data[vars]),] #Remove missing data

##########################################

##########################################

# Perform CCA
  # 5-class solution
cca.results <- cca(data[vars], filter.significance = T,
                   filter.value = 0.05, zero.action = "drop",
                   verbose = F)

  # 4-class solution (see gap statistic section below)
cca.results2 <- cca2(data[vars], filter.significance = T,
                   filter.value = 0.05, zero.action = "drop",
                   verbose = F) #cca2() stops after 6 steps, which corresponds
                                #to a 4-class partition. Could have also stopped
                                #it after 4 or 5 steps, which also correspond to
                                #k = 4

##########################################

# Merge cca group labels with the dataset
cca.group <- cbind(rownames(cca.results2$cormat), 
                   cca.results2$membership) %>% 
  as.data.frame() # Get groups and their row IDs

colnames(cca.group) <- c("rowid", "group")

data$rowid <- rownames(data) #Convert row IDs to a variable

data <- left_join(data, cca.group, by = "rowid") #Merge

##########################################

##########################################

# Descriptive statistics

  # Get percentage breakdown
prop.table(table(data$group))

  # Frequencies, percentages, and medians
data2 <- data[which(is.na(data$group) == F),]
describe(data2[vars])

for (i in vars) {
  print(table(data2[i]))
  print(table(data2[i]) %>% prop.table())
}

##########################################

##########################################

# Network visualizations
cor1 <- cca.results2$modules[[1]]$cormat
cor2 <- cca.results2$modules[[2]]$cormat
cor3 <- cca.results2$modules[[3]]$cormat
cor4 <- cca.results2$modules[[4]]$cormat

g1 <- c(2,5,6,7,8)
g2 <- c(1,3,9)
g3 <- c(4)
col.groups <- list(g1,g2,g3)

  # Dual-Entitlement Omnivores Class
qgraph(cor1,
       graph = "cor",
       minimum=.15, maximum=.75, threshold="sig", sampleSize=157,
        #Not plotting corrs < |.15|. Max non-diagonal abs(corr) = .723.
       groups = col.groups,
       layout = "spring", repulsion = 1.86,
       posCol="black", negCol="black", negDashed=T,
       borders=T, shape = "circle", label.prop = 0.75, 
       curveAll=F, edge.labels=F, edge.label.cex = 0.45, esize = 8,
       title="", color=c("#343d46","#65737e","#c0c5ce"),
       labels=c("pnt","aut","rdw","fxd","pb","grc","app","dll","hrn"),
       label.color=c("black","white","black","black","white","white",
                     "white","white","black"), legend=F
)

  # Procedural (Un)Fairness Class
qgraph(cor2,
       graph = "cor",
       minimum = 0.15, maximum=.75, threshold="sig", sampleSize=389,
       groups = col.groups,
       layout = "spring", repulsion = 1.86,
       posCol="black", negCol="black", negDashed=T,
       borders=T, shape = "circle", label.prop = 0.75, 
       curveAll=F, edge.labels=F, edge.label.cex = 0.45, esize = 8,
       title="", color=c("#343d46","#65737e","#c0c5ce"),
       labels=c("pnt","aut","rdw","fxd","pb","grc","app","dll","hrn"),
       label.color=c("black","white","black","black","white","white",
                     "white","white","black"), legend=F
)

  # Labor-Retail Separate Spheres
qgraph(cor3,
       graph = "cor",
       minimum = 0.15, maximum=.75, threshold="sig", sampleSize=191,
       groups = col.groups,
       layout = "spring", repulsion = 1.86,
       posCol="black", negCol="black", negDashed=T,
       borders=T, shape = "circle", label.prop = 0.75, 
       curveAll=F, edge.labels=F, edge.label.cex = 0.45, esize = 8,
       title="", color=c("#343d46","#65737e","#c0c5ce"),
       labels=c("pnt","aut","rdw","fxd","pb","grc","app","dll","hrn"),
       label.color=c("black","white","black","black","white","white",
                     "white","white","black"), legend=F
)

  # Free-Market Proponents (and Opponents)
qgraph(cor4,
       graph = "cor",
       minimum = 0.15, maximum=.75, threshold="sig", sampleSize=82,
       groups = col.groups,
       layout = "spring", repulsion = 1.25,
       posCol="black", negCol="black", negDashed=T,
       borders=T, shape = "circle", label.prop = 0.75, 
       curveAll=F, edge.labels=F, edge.label.cex = 0.45, esize = 8,
       title="", color=c("#343d46","#65737e","#c0c5ce"),
       labels=c("pnt","aut","rdw","fxd","pb","grc","app","dll","hrn"),
       label.color=c("black","white","black","black","white","white",
                     "white","white","black"), legend=F
)

##########################################

# Correlation matrices
corr.data <- data[which(is.na(data$group) == F),] #Get the 819 cases

colnames(corr.data)[colnames(corr.data)=="redwageb"] <- "redwage"
colnames(corr.data)[colnames(corr.data)=="autoa"] <- "auto" #Rename

vars2 <- c("fixedinc","grocery","pb","apples","auto","dolls","hiring",
           "redwage","painter") #New var list to reflect
                                              #change in names

corr.data1 <- corr.data[which(corr.data$group=="1"),] #Subsets
corr.data2 <- corr.data[which(corr.data$group=="2"),]
corr.data3 <- corr.data[which(corr.data$group=="3"),]
corr.data4 <- corr.data[which(corr.data$group=="4"),]

corr1 <- round(cor(corr.data1[vars2]), 2) #Get bivariate corrs
corr2 <- round(cor(corr.data2[vars2]), 2)
corr3 <- round(cor(corr.data3[vars2]), 2)
corr4 <- round(cor(corr.data4[vars2]), 2)

pmat1 <- cor_pmat(corr.data1[vars2]) #Get corr p-values
pmat2 <- cor_pmat(corr.data2[vars2])
pmat3 <- cor_pmat(corr.data3[vars2])
pmat4 <- cor_pmat(corr.data4[vars2])

  # Dual-Entitlement Omnivores Class
ggcorrplot(corr1, hc.order = F, type = "upper",
                       lab = TRUE, p.mat=pmat1, sig.level=0.05,
                       colors=c("black","white","black")) + 
  scale_y_discrete(limits = rev(unique(colnames(corr1)[colnames(corr1)!="fixedinc"]))) +
  scale_x_discrete(limits = unique(rownames(corr1[-1])))

  # Procedural (Un)Fairness Class
ggcorrplot(corr2, hc.order = F, type = "upper",
                       lab = TRUE, p.mat=pmat2, sig.level=0.05,
                       colors=c("black","white","black")) + 
  scale_y_discrete(limits = rev(unique(colnames(corr2)[colnames(corr2)!="fixedinc"]))) +
  scale_x_discrete(limits = unique(rownames(corr2[-1])))

  # Labor-Retail Separate Spheres
ggcorrplot(corr3, hc.order = F, type = "upper",
                       lab = TRUE, p.mat=pmat3, sig.level=0.05,
                       colors=c("black","white","black")) + 
  scale_y_discrete(limits = rev(unique(colnames(corr3)[colnames(corr3)!="fixedinc"]))) +
  scale_x_discrete(limits = unique(rownames(corr3[-1])))

  # Free-Market Proponents (and Opponents)
ggcorrplot(corr4, hc.order = F, type = "upper",
                       lab = TRUE, p.mat=pmat4, sig.level=0.05,
                       colors=c("black","white","black")) + 
  scale_y_discrete(limits = rev(unique(colnames(corr4)[colnames(corr4)!="fixedinc"]))) +
  scale_x_discrete(limits = unique(rownames(corr4[-1])))

##########################################

##########################################

# Average weighted degree centralities
cor1.net <- cca.results2$modules[[1]]$cormat
diag(cor1.net) <- 0
cor1.adj <- graph.adjacency(cor1.net, mode="undirected", diag=F, weighted=T)

cor2.net <- cca.results2$modules[[2]]$cormat
diag(cor2.net) <- 0
cor2.adj <- graph.adjacency(cor2.net, mode="undirected", diag=F, weighted=T)

cor3.net <- cca.results2$modules[[3]]$cormat
diag(cor3.net) <- 0
cor3.adj <- graph.adjacency(cor3.net, mode="undirected", diag=F, weighted=T)

cor4.net <- cca.results2$modules[[4]]$cormat
diag(cor4.net) <- 0
cor4.adj <- graph.adjacency(cor4.net, mode="undirected", diag=F, weighted=T)

degree1 <- (strength(cor1.adj))/8
degree2 <- (strength(cor2.adj))/8
degree3 <- (strength(cor3.adj))/8
degree4 <- (strength(cor4.adj))/8

degree1.mat <- as.matrix(degree1)
degree2.mat <- as.matrix(degree2)
degree3.mat <- as.matrix(degree3)
degree4.mat <- as.matrix(degree4)

degree <- cbind(degree1,degree2,degree3,degree4)

degree <- melt(degree, id=rownames(degree))

pd <- position_dodge(0.1)
ggplot(data=degree, aes(x=Var1, y=value, group=Var2, shape=Var2)) +
  geom_line(aes(linetype=Var2), show.legend=F, position=pd) +
  geom_point(aes(shape=Var2), position=pd) +
  theme_bw() +
  theme(axis.line=element_line(),
        axis.title.y=element_text(size=10), legend.title=element_text(face="bold"),
        legend.position="bottom") +
  xlab("Scenario") + ylab("Average Correlation with the Scenario") +
  scale_shape_discrete(name="",breaks=c("degree1","degree4","degree3","degree2"),
                       labels=c("Class A:\nDual-Entitlement Omnivores","Class B:\nAnything (But) Non-Competition",
                                "Class C:\nLabor-Retail Separate Spheres",
                                "Class D:\nSeparate Spheres with\nProcedural (Un)Fairness")) +
  scale_x_discrete(labels=c("painter"="pnt",
                            "autoa"="aut","redwageb"="rdw","fixedinc"="fxd",
                            "pb"="pb","apples"="app",
                            "dolls"="dll","hiring"="hrn","grocery"="grc")) +
  guides(shape=guide_legend(nrow=2,byrow=T))


  #Testing some of these within-scenario/across-class average 
    #weighted degree centralities (substitute in necessary numbers):
r.test(n = 157, n2 = 191, r12 = .4, r34 = .08, twotailed = T)

##########################################

##########################################

# Multiple Group Test of Schema Invariance - did this in Stata 
  # (see fairness_sem.do)
write.csv(data, file="sem_data.csv")

##########################################

##########################################

# Gap Statistic

  # Comparing modularity scores
graph.adjacency(cca.results2$cormat, 
                mode="undirected", weighted=T) %>%
  leading.eigenvector.community(., steps = 6) #Q = .33 when k = 4
                                              #Q = .36 when k = 5,
                                              #which occurs at 7 steps

gap.data <- corr.data[vars2] #Subset

set.seed(123)
myCluster <- function(x, k) list(cluster = leading.eigenvector.community(
  graph.adjacency(cca(x, filter.significance=T, filter.value=0.05, 
                      zero.action=c("drop"), verbose=F)$cormat, 
                  mode="undirected", weighted=T), steps=k)$membership)
gap <- clusGap2(gap.data, K.max=7, B=500, FUN=myCluster)

gap.stats <- print(gap, method="Tibs2001SEmax")$Tab[c(1,2,3,4,8),] %>% 
  as.data.frame()
gap.stats$cluster <- rownames(gap.stats)

  # Gap stat comparison
p1 <- ggplot(gap.stats, aes(x=cluster, y=gap, group=1)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=gap-SE.sim, ymax=gap+SE.sim), width=.1) +
  geom_vline(xintercept = 4, linetype="dashed", color = "black") +
  xlab("Number of Classes") + ylab("E*(log[W]) - log[W]") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), 
        axis.line=element_line(), panel.border=element_blank())

  # Gap stat differences 
gap.stats$diff <- NA

for (i in 1:4) {
  gap.stats[i,6] <- gap.stats[i,3] - 
    (gap.stats[(i+1),3] - gap.stats[(i+1),4])
}

p2 <- ggplot(gap.stats[-5,], aes(x=cluster, y=diff)) + 
  geom_bar(stat="identity", aes(fill = diff>0), position="dodge", fill="gray50", col="black") +
  xlab("Number of Classes") + ylab("Gap(k) - [Gap(k+1) - SE(k+1)]") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), 
        axis.line=element_line(), panel.border=element_blank()) +
  scale_fill_discrete(guide="none")

ggarrange(p1, p2, nrow=1, ncol=2) #Preferred k = 4

##########################################

### END ###
