#############################
### Analysis with Model-based Clustering
#############################
rm(list=ls())
library(sampling)
library(sm)
library(mclust)
library(xlsx)
library(reshape)
library(ggplot2)
#library(scatterplot3d)
#library(rgl)

######################################################

##################
## Screening Test: Phase-I
################## 

### Read Data ###
furn = read.csv("demoa_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

### Principle Components ###
pca.furn <- prcomp(furn1, center = T, scale = T)
x = pca.furn$x

### Scree Plot ###
library(ggplot2)
variances <- data.frame(variances=pca.furn$sdev**2, principal_component=1:length(pca.furn$sdev))
variances[,1] <- variances[,1]/sum(variances[,1])
varPlot <- ggplot(variances, aes(principal_component, variances)) + geom_bar(stat="identity", fill="gray") + geom_line() + labs(x="Principal Components",y="Proportion of Variance")
varPlot
ggsave("elbow_wave1_est.pdf", height=5, width=5, units='in')

### Cumulative Plot ###
variances[,1] <- cumsum(variances[,1])/sum(variances[,1])
varPlot <- ggplot(variances, aes(principal_component, variances)) + geom_bar(stat="identity", fill="gray") + geom_line()
varPlot

### Clustering ###
Stest <- Mclust(x[,1:10] , G=2)
classify = Stest$classification

Y = data.frame(cbind(EDEPI, classify))
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
f2<-c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n<-round(n2*f2[j])
  preval.t = NULL
  for(i in 1:500){
    s = srswor(n, n2)
    s = c(s1,s)
    k = which(s==1)
    #length(k)
    #rm(list="s")	
    
    ### Sample for GS ###
    new<-Y[k,]
    #comp<- cbind(temp[,1][s==0], temp[,2][s==0], temp[,3][s==0])
    colnames(new) <- c("T","P")
    new<- as.data.frame(new)	
    
    ### Calculation of Lambda ###
    (tb<-table(P=new$P, T=new$T))
    lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
    lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])
    
    ### Prevalence and Variance ###
    preval.t <- c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  }
  preval = rbind(preval, c(length(k), preval.t))
}

rmean = rowMeans(preval[,-1])
rsd = apply(preval[,-1],1,sd)	
tab.val = cbind(rep(dim(furn)[1],3), preval[,1], rmean, rsd, rep(NA,3), rep(NA,3))

preval = t(preval[,-1])

### Prevalence Box Plots ###
colnames(preval) = c("5%","10%","20%")

plt = ggplot(data=melt(preval),aes(as.character(X2),value,fill=as.factor(X2))) + geom_boxplot() + scale_fill_manual(values = c("#99FFCC", "#99CCFF", "#CCCCCC")) + theme(legend.position="none") + scale_x_discrete(limits=c("5%","10%","20%"))
plt = plt + stat_summary(fun.y=mean, geom="point", shape=17, size=2, col="red") + geom_hline(yintercept=(length(which(EDEPI == 1))/length(EDEPI)),col="red",linetype="dashed")
plt + labs(x="Percent of Negative Screened Individuals",y="Predicted Prevalence") + ylim(0.0,0.4)
ggsave("plot_mclust_wave1_est.pdf", height=5, width=5, units='in')



##################
## Screening Test: Phase-II
################## 

### Read Data ###
furn = read.csv("demob_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

### Principle Components ###
pca.furn <- prcomp(furn1, center = T, scale = T)
x = pca.furn$x

### Scree Plot ###
library(ggplot2)
variances <- data.frame(variances=pca.furn$sdev**2, principal_component=1:length(pca.furn$sdev))
variances[,1] <- variances[,1]/sum(variances[,1])
varPlot <- ggplot(variances, aes(principal_component, variances)) + geom_bar(stat="identity", fill="gray") + geom_line() + labs(x="Principal Components",y="Proportion of Variance")
varPlot
ggsave("elbow_wave2_est.pdf", height=5, width=5, units='in')


### Cumulative Plot ###
variances[,1] <- cumsum(variances[,1])/sum(variances[,1])
varPlot <- ggplot(variances, aes(principal_component, variances)) + geom_bar(stat="identity", fill="gray") + geom_line()
varPlot

### Clustering ###
Stest <- Mclust(x[,1:10] , G=2)
classify = Stest$classification

Y = data.frame(cbind(EDEPI, classify))
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
f2<-c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n<-round(n2*f2[j])
  preval.t = NULL
  for(i in 1:500){
    s = srswor(n, n2)
    s = c(s1,s)
    k = which(s==1)

    ### Sample for GS ###
    new<-Y[k,]
    colnames(new) <- c("T","P")
    new<- as.data.frame(new)	
    
    ### Calculation of Lambda ###
    (tb<-table(P=new$P, T=new$T))
    lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
    lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])
    
    ### Prevalence and Variance ###
    preval.t <- c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  }
  preval = rbind(preval, c(length(k), preval.t))
}

### Incidence Calculation ###

pt2 = rowMeans(preval[,-1]);vpt2 = apply(preval[,-1],1,var);pt1 = tab.val[,3]; vpt1 = tab.val[,4]^2;
Nt2 = dim(furn)[1]
Nt1 = dim(furn)[1] + round(tab.val[1,1]*pt1)
incid = 1- (Nt2*(1-pt2))/(Nt1*(1-pt1))
sd.incid = sqrt(((1-pt2)^2*vpt1 + vpt2*(1 - vpt1))/(1 - pt1)^2)
tab.val = rbind(tab.val, cbind(rep(dim(furn)[1],3), preval[,1], pt2, sqrt(vpt2), incid, sd.incid))

preval = t(preval[,-1])

### Prevalence Box Plots ###
colnames(preval) = c("5%","10%","20%")

plt = ggplot(data=melt(preval),aes(as.character(X2),value,fill=as.factor(X2))) + geom_boxplot() + scale_fill_manual(values = c("#99FFCC", "#99CCFF", "#CCCCCC")) + theme(legend.position="none") + scale_x_discrete(limits=c("5%","10%","20%"))
plt = plt + stat_summary(fun.y=mean, geom="point", shape=17, size=2, col="red") + geom_hline(yintercept=(length(which(EDEPI == 1))/length(EDEPI)),col="red",linetype="dashed")
plt + labs(x="Percent of Negative Screened Individuals",y="Predicted Prevalence") + ylim(0.0,0.4)
ggsave("plot_mclust_wave2_est.pdf", height=5, width=5, units='in')


##################
## Screening Test: Phase-III
################## 

### Read Data ###
furn = read.csv("democ_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

### Principle Components ###
pca.furn <- prcomp(furn1, center = T, scale = T)
x = pca.furn$x

### Scree Plot ###
library(ggplot2)
variances <- data.frame(variances=pca.furn$sdev**2, principal_component=1:length(pca.furn$sdev))
variances[,1] <- variances[,1]/sum(variances[,1])
varPlot <- ggplot(variances, aes(principal_component, variances)) + geom_bar(stat="identity", fill="gray") + geom_line() + labs(x="Principal Components",y="Proportion of Variance")
varPlot
ggsave("elbow_wave3_est.pdf", height=5, width=5, units='in')

### Cumulative Plot ###
variances[,1] <- cumsum(variances[,1])/sum(variances[,1])
varPlot <- ggplot(variances, aes(principal_component, variances)) + geom_bar(stat="identity", fill="gray") + geom_line()
varPlot

### Clustering ###
Stest <- Mclust(x[,1:10] , G=2)
classify = Stest$classification

Y = data.frame(cbind(EDEPI, classify))
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
f2<-c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n<-round(n2*f2[j])
  preval.t = NULL
  for(i in 1:500){
    s = srswor(n, n2)
    s = c(s1,s)
    k = which(s==1)
    #length(k)
    #rm(list="s")	
    
    ### Sample for GS ###
    new<-Y[k,]
    #comp<- cbind(temp[,1][s==0], temp[,2][s==0], temp[,3][s==0])
    colnames(new) <- c("T","P")
    new<- as.data.frame(new)	
    
    ### Calculation of Lambda ###
    (tb<-table(P=new$P, T=new$T))
    lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
    lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])
    
    ### Prevalence and Variance ###
    preval.t <- c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  }
  preval = rbind(preval, c(length(k), preval.t))
}

### Incidence Calculation ###

pt2 = rowMeans(preval[,-1]);vpt2 = apply(preval[,-1],1,var);pt1 = tab.val[4:6,3]; vpt1 = tab.val[4:6,4]^2;
Nt2 = dim(furn)[1]
Nt1 = dim(furn)[1] + round(tab.val[4,1]*pt1)
incid = 1- (Nt2*(1-pt2))/(Nt1*(1-pt1))
sd.incid = sqrt(((1-pt2)^2*vpt1 + vpt2*(1 - vpt1))/(1 - pt1)^2)
tab.val = rbind(tab.val, cbind(rep(dim(furn)[1],3), preval[,1], pt2, sqrt(vpt2), incid, sd.incid))

colnames(tab.val) = c("#obs","Sample","Predicted Prevalence","SD Prevalence","Incidence Rate","SD Incidence")
rownames(tab.val) = c(rep("Wave-I",3),rep("Wave-II",3),rep("Wave-III",3))
write.csv(tab.val,"Stagewise_Prevalence_Mclust_est.csv")
preval = t(preval[,-1])

### Prevalence Box Plots ###
colnames(preval) = c("5%","10%","20%")

plt = ggplot(data=melt(preval),aes(as.character(X2),value,fill=as.factor(X2))) + geom_boxplot() + scale_fill_manual(values = c("#99FFCC", "#99CCFF", "#CCCCCC")) + theme(legend.position="none") + scale_x_discrete(limits=c("5%","10%","20%"))
plt = plt + stat_summary(fun.y=mean, geom="point", shape=17, size=2, col="red") + geom_hline(yintercept=(length(which(EDEPI == 1))/length(EDEPI)),col="red",linetype="dashed")
plt + labs(x="Percent of Negative Screened Individuals",y="Predicted Prevalence") + ylim(0.0,0.4)
ggsave("plot_mclust_wave3_est.pdf", height=5, width=5, units='in')