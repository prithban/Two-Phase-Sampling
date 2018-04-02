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
furn <- read.xlsx("C:/Users/pbanerjee/Documents/Two_Phase_Sampling/Data/demoa_imputed.xls",sheetIndex=1);

### Variables for Principal Component ###
furn1 <- furn[,c("ASRH1", "Gender", "MARITAL", "EDUCATION",	"AGE", "ADL", "IADL", "Mobility", "AMMSE", "ACHARLSD", "POVERTY", 
                 "RACEP", "SMOKEYR", "ABMI")]

### Principle Components ###
pca.furn <- prcomp(furn1, center = T, scale = T)
x = pca.furn$x

### Clustering ###
Stest <- Mclust(x[,1:10] , G=2)

Y <- cbind(furn[,1:2],furn$ADEPI,Stest$classification)[,-2];
Y <- Y[order(Y[,3]),]

### Calculating the screening test sensitivity and specificity
sensitivity.screen = specifisity.screen = NULL
tab.screen = table(Y$`Stest$classification`, Y$`furn$ADEPI`)
sensitivity.screen = c(sensitivity.screen, tab.screen[1,1]/sum(tab.screen[,1]))
specifisity.screen = c(specifisity.screen, tab.screen[2,2]/sum(tab.screen[,2]))

### Choice of Sample/Sample Proportion ###
f1<-1; ## 1 MDD, 2 Non-MDD
n1<- sum(Y[,3]==1);  n2<- sum(Y[,3]==2); s1<-rep(1,n1);
pi1<-n1/length(Y[,3]);


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
  #for(i in 1:500){
  s = srswor(n, n2)
  s = c(s1,s)
  k = which(s==1)
  #length(k)
  #rm(list="s")	
  
  ### Sample for GS ###
  new<-Y[k,]
  #comp<- cbind(temp[,1][s==0], temp[,2][s==0], temp[,3][s==0])
  colnames(new) <- c("X","T","P")
  new<- as.data.frame(new)	
  
  ### Calculation of Lambda ###
  (tb<-table(P=new$P, T=new$T))
  #lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
  #lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])
  sense.GS = c(sense.GS, tb[1,1]/(tb[1,1]+tb[2,1]))
  spec.GS = c(spec.GS, tb[2,2]/(tb[1,2]+tb[2,2]))
  ### Prevalence and Variance ###
  #preval.t <- c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  #}
  #preval = rbind(preval, c(length(k), preval.t))
}

sense.test = rbind(sense.test, sense.GS)
spec.test = rbind(spec.test, spec.GS)


##################
## Screening Test: Phase-II
################## 

### Read Data ###
furn <- read.xlsx("C:/Users/pbanerjee/Documents/Two_Phase_Sampling/Data/democ_imputed.xls",sheetIndex=1);

### Variables for Principal Component ###
furn1 <- furn[,c("CSRH1", "Gender", "MARITAL", "EDUCATION",	"AGE", "ADL", "IADL", "Mobility", "AMMSE", "ACHARLSD", "POVERTY", 
                 "RACEP", "SMOKEYR", "ABMI")]

### Principle Components ###
pca.furn <- prcomp(furn1, center = T, scale = T)
x = pca.furn$x

### Clustering ###
Stest <- Mclust(x[,1:9] , G=2)

Y <- cbind(furn[,1:2],furn$CDEPI,Stest$classification)[,-2];
Y <- Y[order(Y[,3]),] 

### Calculating the screening test sensitivity and specificity
tab.screen = table(Y$`Stest$classification`, Y$`furn$CDEPI`)
sensitivity.screen = c(sensitivity.screen, tab.screen[1,1]/sum(tab.screen[,1]))
specifisity.screen = c(specifisity.screen, tab.screen[2,2]/sum(tab.screen[,2]))

### Choice of Sample/Sample Proportion ###
f1<-1; ## 1MDD, 2 Non-MDD
n1<- sum(Y[,3]==1);  n2<- sum(Y[,3]==2); s1<-rep(1,n1);
pi1<-n1/length(Y[,2]);


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
  #length(k)
  #rm(list="s")	
  
  ### Sample for GS ###
  new<-Y[k,]
  #comp<- cbind(temp[,1][s==0], temp[,2][s==0], temp[,3][s==0])
  colnames(new) <- c("X","T","P")
  new<- as.data.frame(new)	
  
  ### Calculation of Lambda ###
  (tb<-table(P=new$P, T=new$T))
  #lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
  #lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])
  sense.GS = c(sense.GS, tb[1,1]/(tb[1,1]+tb[2,1]))
  spec.GS = c(spec.GS, tb[2,2]/(tb[1,2]+tb[2,2]))
  ### Prevalence and Variance ###
  #preval.t <- c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  #}
  #preval = rbind(preval, c(length(k), preval.t))
}

sense.test = rbind(sense.test, sense.GS)
spec.test = rbind(spec.test, spec.GS)


##################
## Screening Test: Phase-III
################## 

### Read Data ###
furn <- read.xlsx("C:/Users/pbanerjee/Documents/Two_Phase_Sampling/Data/demob_imputed.xls",sheetIndex=1);

### Variables for Principal Component ###
furn1 <- furn[,c("BSRH1", "Gender", "MARITAL", "EDUCATION",	"AGE", "ADL", "IADL", "Mobility", "AMMSE", "ACHARLSD", "POVERTY", 
                 "RACEP", "SMOKEYR", "ABMI")]

### Principle Components ###
pca.furn <- prcomp(furn1, center = T, scale = T)
x = pca.furn$x

### Clustering ###
Stest <- Mclust(x[,1:9] , G=2)

Y <- cbind(furn[,1:2],furn$BDEPI,Stest$classification)[,-2];
Y <- Y[order(Y[,3]),] 

### Calculating the screening test sensitivity and specificity
tab.screen = table(Y$`Stest$classification`, Y$`furn$BDEPI`)
sensitivity.screen = c(sensitivity.screen, tab.screen[1,1]/sum(tab.screen[,1]))
specifisity.screen = c(specifisity.screen, tab.screen[2,2]/sum(tab.screen[,2]))

### Choice of Sample/Sample Proportion ###
f1<-1; ## 1MDD, 2 Non-MDD
n1<- sum(Y[,3]==1);  n2<- sum(Y[,3]==2); s1<-rep(1,n1);
pi1<-n2/length(Y[,2]);

### Simple Random Sampling ###
set.seed(1000)
preval = NULL

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
  #length(k)
  #rm(list="s")	
  
  ### Sample for GS ###
  new<-Y[k,]
  #comp<- cbind(temp[,1][s==0], temp[,2][s==0], temp[,3][s==0])
  colnames(new) <- c("X","T","P")
  new<- as.data.frame(new)	
  
  ### Calculation of Lambda ###
  (tb<-table(P=new$P, T=new$T))
  #lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
  #lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])
  sense.GS = c(sense.GS, tb[1,1]/(tb[1,1]+tb[2,1]))
  spec.GS = c(spec.GS, tb[2,2]/(tb[1,2]+tb[2,2]))
  ### Prevalence and Variance ###
  #preval.t <- c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  #}
  #preval = rbind(preval, c(length(k), preval.t))
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
