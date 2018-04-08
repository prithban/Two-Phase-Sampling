#############################
### Analysis with Hierarchical Clustering
#############################
rm(list=ls())
library(sampling)
library(sm)
library(mclust)
library(xlsx)
library(reshape)
library(ggplot2)
library(ggdendro)
library(dendextend)

##################
## Screening Test: Wave-I
##################

### Read Data ###
furn = read.csv("demoa_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

### Clustering ###
clusters = furn1 %>%  dist %>% hclust %>% as.dendrogram

### Ordering Sample ###
clustercut = cutree(clusters, k=2)

Y = data.frame(cbind(EDEPI, clustercut))
Y = Y[order(Y[,2]),]

### Calculating the screening test sensitivity and specificity
sensitivity.screen = specifisity.screen = NULL
tab.screen = table(classify, Y$EDEPI)
sensitivity.screen = c(sensitivity.screen, tab.screen[1,1]/sum(tab.screen[,1]))
specifisity.screen = c(specifisity.screen, tab.screen[2,2]/sum(tab.screen[,2]))

### Choice of Sample/Sample Proportion ###
ppv = NULL
f1=1; ## 1 MDD, 2 Non-MDD
n1= sum(Y[,2]==1);  n2= sum(Y[,2]==2); s1=rep(1,n1);
pi1=n1/length(Y[,2]);
ppv = c(ppv, pi1)

### Simple Random Sampling ###
set.seed(1000)
preval = NULL

### Sampling Proportions ###
sense.test = spec.test = rep(0,3)
sense.GS = spec.GS = NULL
f2<-c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n<-round(n2*f2[j])
  preval.t = NULL
  s = srswor(n, n2)
  s = c(s1,s)
  k = which(s==1)
  
  ### Sample for GS ###
  new<-Y[k,]
  colnames(new) <- c("T","P")
  new<- as.data.frame(new)	
  
  ### Calculation of Lambda ###
  (tb<-table(P=new$P, T=new$T))
  
  sense.GS = c(sense.GS, tb[1,1]/(tb[1,1]+tb[2,1]))
  spec.GS = c(spec.GS, tb[2,2]/(tb[1,2]+tb[2,2]))
}

sense.test = rbind(sense.test, sense.GS)
spec.test = rbind(spec.test, spec.GS)


##################
## Screening Test: Phase-II
################## 

### Read Data ###
furn = read.csv("demob_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

### Clustering ###
clusters = furn1 %>%  dist %>% hclust %>% as.dendrogram

### Ordering Sample ###
clustercut = cutree(clusters, k=2)

Y = data.frame(cbind(EDEPI, clustercut))
Y = Y[order(Y[,2]),]

### Calculating the screening test sensitivity and specificity
sensitivity.screen = specifisity.screen = NULL
tab.screen = table(classify, Y$EDEPI)
sensitivity.screen = c(sensitivity.screen, tab.screen[1,1]/sum(tab.screen[,1]))
specifisity.screen = c(specifisity.screen, tab.screen[2,2]/sum(tab.screen[,2]))

### Choice of Sample/Sample Proportion ###
ppv = NULL
f1=1; ## 1 MDD, 2 Non-MDD
n1= sum(Y[,2]==1);  n2= sum(Y[,2]==2); s1=rep(1,n1);
pi1=n1/length(Y[,2]);
ppv = c(ppv, pi1)

### Sampling Proportions ###
sense.GS = spec.GS = NULL
f2<-c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n<-round(n2*f2[j])
  preval.t = NULL
  #for(i in 1:500){
  s = srswor(n, n2)
  s = c(s1,s)
  k = which(s==1)
  
  ### Sample for GS ###
  new<-Y[k,]
  colnames(new) <- c("T","P")
  new<- as.data.frame(new)	
  
  ### Calculation of Lambda ###
  (tb<-table(P=new$P, T=new$T))
  sense.GS = c(sense.GS, tb[1,1]/(tb[1,1]+tb[2,1]))
  spec.GS = c(spec.GS, tb[2,2]/(tb[1,2]+tb[2,2]))
  
}

sense.test = rbind(sense.test, sense.GS)
spec.test = rbind(spec.test, spec.GS)


##################
## Screening Test: Phase-III
################## 

### Read Data ###
furn = read.csv("democ_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

### Clustering ###
clusters = furn1 %>%  dist %>% hclust %>% as.dendrogram

### Ordering Sample ###
clustercut = cutree(clusters, k=2)

Y = data.frame(cbind(EDEPI, clustercut))
Y = Y[order(Y[,2]),]

### Calculating the screening test sensitivity and specificity
sensitivity.screen = specifisity.screen = NULL
tab.screen = table(classify, Y$EDEPI)
sensitivity.screen = c(sensitivity.screen, tab.screen[1,1]/sum(tab.screen[,1]))
specifisity.screen = c(specifisity.screen, tab.screen[2,2]/sum(tab.screen[,2]))

### Choice of Sample/Sample Proportion ###
ppv = NULL
f1=1; ## 1 MDD, 2 Non-MDD
n1= sum(Y[,2]==1);  n2= sum(Y[,2]==2); s1=rep(1,n1);
pi1=n1/length(Y[,2]);
ppv = c(ppv, pi1)

### Simple Random Sampling ###
set.seed(1000)
preval = NULL

### Sampling Proportions ###
sense.GS = spec.GS = NULL
f2<-c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n<-round(n2*f2[j])
  preval.t = NULL
  s = srswor(n, n2)
  s = c(s1,s)
  k = which(s==1)
  
  ### Sample for GS ###
  new<-Y[k,]
  colnames(new) <- c("X","T","P")
  new<- as.data.frame(new)	
  
  ### Calculation of Lambda ###
  (tb<-table(P=new$P, T=new$T))
  sense.GS = c(sense.GS, tb[1,1]/(tb[1,1]+tb[2,1]))
  spec.GS = c(spec.GS, tb[2,2]/(tb[1,2]+tb[2,2]))
}

sense.test = rbind(sense.test, sense.GS)
spec.test = rbind(spec.test, spec.GS)

### Final sensitivity and specificity

sense.test = sense.test[-1,]
colnames(sense.test) = c("5%", "10%", "20%")
rownames(sense.test) = c("Baseline", "3 Month", "12 Month")

spec.test = spec.test[-1,]
colnames(spec.test) = c("5%", "10%", "20%")
rownames(spec.test) = c("Baseline", "3 Month", "12 Month")
