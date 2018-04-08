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

tab.val = NULL

### Read Data ###
furn = read.csv("demoa_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

##############
### Without PCA
##############

### Clustering ###
clusters = furn1 %>%  dist %>% hclust %>% as.dendrogram
clustercut = cut(clusters,h=27)

### Dendrogram ###
pdf('dendogram.hclust_wave1_est.pdf')
clustercut$upper %>% set("branches_k_color", value = c("hotpink2","lightseagreen"), k = 2) %>% set("branches_lwd", 2.0) %>% set("labels", "") %>% hang.dendrogram(hang = -1) %>% plot()
clusters %>% rect.dendrogram(k = 2, border = 8, lty = 5, lwd = 2, prop_k_height = -1.0)
legend("topright",c("MDD","NON-MDD"), fill=c("hotpink2","lightseagreen"))
dev.off()

### Ordering Sample ###
clustercut = cutree(clusters, k=2)

Y = data.frame(cbind(EDEPI, clustercut))
Y = Y[order(Y[,2]),]

### Calculating the screening test sensitivity and specificity
sensitivity.screen = specifisity.screen = NULL
tab.screen = table(Y$clustercut, Y$EDEPI)
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
f2=c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n=round(n2*f2[j])
  preval.t = NULL
  for(i in 1:500){
    s = srswor(n, n2)
    s = c(s1,s)
    k = which(s==1)
    length(k)
    rm(list="s")	
    
    ### Sample for GS ###
    new=Y[k,]
    colnames(new) = c("T","P")
    new= as.data.frame(new)	
    
    ### Calculation of Lambda ###
    (tb=table(P=new$P, T=new$T))
    lt1= tb[1,1]/(tb[1,1]+tb[1,2])
    lt2= tb[2,1]/ (tb[2,1]+tb[2,2])
    
    ### Prevalence and Variance ###
    preval.t = c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  }
  preval = rbind(preval, c(length(k),preval.t))
}

rmean = rowMeans(preval[,-1])
rsd = apply(preval[,-1],1,sd)	
tab.val = cbind(rep(dim(furn)[1],3), preval[,1], rmean, rsd, rep(NA,3), rep(NA,3))

preval = t(preval[,-1])

### Prevalence Box Plots ###
colnames(preval) = c("05%","10%","20%")


plt = ggplot(data=melt(preval),aes(as.character(X2),value,fill=as.factor(X2))) + geom_boxplot() + scale_fill_brewer(palette="Blues") + theme(legend.position="none") 
plt = plt + stat_summary(fun.y=mean, geom="point", shape=17, size=2, col="red") + geom_hline(yintercept=0.1595,col="red",linetype="dashed")
plt + labs(title="Wave-I Boxplots for Prevalence(HCLUST)",x="Percent of Negative Screened Individuals",y="Predicted Prevalence") + ylim(0.0,0.4)
ggsave("plot_hclust_wave1_est.pdf", height=5, width=5, units='in')



##################
## Screening Test: Wave-II
##################


### Read Data ###
furn = read.csv("demob_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

##############
### Without PCA
##############

### Clustering ###
clusters = furn1 %>%  dist %>% hclust %>% as.dendrogram
clustercut = cut(clusters,h=25)

### Dendrogram ###
pdf('dendogram.hclust_wave2_est.pdf')
clustercut$upper %>% set("branches_k_color", value = c("hotpink2","lightseagreen"), k = 2) %>% set("branches_lwd", 2.0) %>% set("labels", "") %>% hang.dendrogram(hang = -1) %>% plot(main="Wave-II Screening Test Dendrogram")
clusters %>% rect.dendrogram(k = 2, border = 8, lty = 5, lwd = 2, prop_k_height = -1.6)
legend("topright",c("MDD","NON-MDD"), fill=c("hotpink2","lightseagreen"))
dev.off()

### Ordering Sample ###
clustercut = cutree(clusters, k=2)

Y = data.frame(cbind(EDEPI, clustercut))
Y = Y[order(Y[,2]),]

### Calculating the screening test sensitivity and specificity
sensitivity.screen = specifisity.screen = NULL
tab.screen = table(Y$clustercut, Y$EDEPI)
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
f2=c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n=round(n2*f2[j])
  preval.t = NULL
  for(i in 1:500){
    s = srswor(n, n2)
    s = c(s1,s)
    k = which(s==1)
    length(k)
    rm(list="s")	
    
    ### Sample for GS ###
    new=Y[k,]
    colnames(new) = c("T","P")
    new= as.data.frame(new)	
    
    ### Calculation of Lambda ###
    (tb=table(P=new$P, T=new$T))
    lt1= tb[1,1]/(tb[1,1]+tb[1,2])
    lt2= tb[2,1]/ (tb[2,1]+tb[2,2])
    
    ### Prevalence and Variance ###
    preval.t = c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  }
  preval = rbind(preval, c(length(k),preval.t))
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
colnames(preval) = c("05%","10%","20%")

plt = ggplot(data=melt(preval),aes(as.character(X2),value,fill=as.factor(X2))) + geom_boxplot() + scale_fill_brewer(palette="Blues") + theme(legend.position="none") 
plt = plt + stat_summary(fun.y=mean, geom="point", shape=17, size=2, col="red") + geom_hline(yintercept=0.1022,col="red",linetype="dashed")
plt + labs(title="Wave-II Boxplots for Prevalence(HCLUST)",x="Percent of Negative Screened Individuals",y="Predicted Prevalence") + ylim(0.0,0.4)
ggsave("plot_hclust_wave2_est.pdf", height=5, width=5, units='in')



##################
## Screening Test: Wave-III
##################

### Read Data ###
furn = read.csv("democ_sim.csv",header = T);

### Variables for Principal Component ###
furn1 = furn[,-dim(furn)[2]]
EDEPI = furn$EDEPI

##############
### Without PCA
##############

### Clustering ###
clusters = furn1 %>%  dist %>% hclust %>% as.dendrogram
clustercut = cut(clusters,h=20)

### Dendrogram ###
pdf('dendogram.hclust_wave3_est.pdf')
clustercut$upper %>% set("branches_k_color", value = c("hotpink2","lightseagreen"), k = 2) %>% set("branches_lwd", 2.0) %>% set("labels", "") %>% hang.dendrogram(hang = -1) %>% plot(main="Wave-III Screening Test Dendrogram")
clusters %>% rect.dendrogram(k = 2, border = 8, lty = 5, lwd = 2, prop_k_height = -1.58)
legend("topright",c("MDD","NON-MDD"), fill=c("hotpink2","lightseagreen"))
dev.off()

### Ordering Sample ###
clustercut = cutree(clusters, k=2)

Y = data.frame(cbind(EDEPI, clustercut))
Y = Y[order(Y[,2]),]

### Calculating the screening test sensitivity and specificity
sensitivity.screen = specifisity.screen = NULL
tab.screen = table(Y$clustercut, Y$EDEPI)
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
f2=c(0.05,0.1,0.2)
for(j in 1:length(f2)){
  n=round(n2*f2[j])
  preval.t = NULL
  for(i in 1:500){
    s = srswor(n, n2)
    s = c(s1,s)
    k = which(s==1)
    length(k)
    rm(list="s")	
    
    ### Sample for GS ###
    new=Y[k,]
    colnames(new) = c("T","P")
    new= as.data.frame(new)	
    
    ### Calculation of Lambda ###
    (tb=table(P=new$P, T=new$T))
    lt1= tb[1,1]/(tb[1,1]+tb[1,2])
    lt2= tb[2,1]/ (tb[2,1]+tb[2,2])
    
    ### Prevalence and Variance ###
    preval.t = c(preval.t, pi1*lt1 +(1- pi1)*lt2)
  }
  preval = rbind(preval, c(length(k),preval.t))
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
write.csv(tab.val,"Stagewise_Prevalence_Hclust_est.csv")

preval = t(preval[,-1])

### Prevalence Box Plots ###
colnames(preval) = c("05%","10%","20%")

plt = ggplot(data=melt(preval),aes(as.character(X2),value,fill=as.factor(X2))) + geom_boxplot() + scale_fill_brewer(palette="Blues") + theme(legend.position="none") 
plt = plt + stat_summary(fun.y=mean, geom="point", shape=17, size=2, col="red") + geom_hline(yintercept=0.1569,col="red",linetype="dashed")
plt + labs(title="Wave-III Boxplots for Prevalence(HCLUST)",x="Percent of Negative Screened Individuals",y="Predicted Prevalence") + ylim(0.0,0.4)
ggsave("plot_hclust_wave3_est.pdf", height=5, width=5, units='in')