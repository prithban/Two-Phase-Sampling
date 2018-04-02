########################################################
#### Code for New Simulation Study
########################################################

#########################
###Informants score degrades with time with fixed incidence rate (new model 3).
#########################
rm(list=ls())

library("sampling")
set.seed(100)
### Assign  baseline
m<-1 #baseline

### Normal group
x<-as.vector(rnorm(900,2,1.5))  ##feature sample
D<-rep(2,length(x))             ##response 
x1<-cbind(x,D)                  ##normal group data

### Demented group
x<-as.vector(rnorm(100,-2,4))   ##feature sample
D <-rep(1, length(x))           ##response
x2<-cbind(x, D)                 ##demented group data

### Final data
z<-rbind(x1,x2)

### Screening test for data partition 
cl <- kmeans(z[,1], 2, iter.max = 200, nstart = 10, algorithm = "MacQueen")

### Plots for visualization
plot(z[,1], col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)

### Predicted partition
Y<-as.integer(cl$cluster)

center<- sort(cl$centers);
if (center[1]!= cl$centers[1])  {for (i in 1:length(Y)) {if (Y[i]==2) Y[i]=1 else Y[i]=2} }

### Checking the correct classification
sum(z[,2]== Y)

### COmbining original labels and predicted labels
temp<-cbind(z, Y)
temp<- temp[order(temp[,3]),]

### Calculating different quantities of two-phase sampling

## Assigning proportion of observations chosen from each group
f2<-0.1;f1<-1;                                

## Sample sizes and sample object
n1<- sum(Y==2);n2<- sum(Y==1);s1<-rep(1,n2);

## MLE of demented group
pi1<-n2/length(Y);

## Sample size from non-demented by screening test
n<-round(n1*f2)

### simple random sampling from non-demented group
s=srswor(n, n1)
s<- c(s1,s)

### Final sample for GS test
new <- temp[which(s==1),]
comp<- temp[which(s==0),]

### Data for GS test
colnames(new) <- c("X","D","Y")
new<- as.data.frame(new)
tb<-table(new$Y, new$D)

### Storing results for baseline
result <- matrix(NA, nrow = 11, ncol=9, byrow=TRUE, dimnames = list(c("baseline", "wave1", "wave2", "wave3","wave4", "wave5", "wave6","wave7", "wave8", "wave9", "wave10"), c("prev", "stdpre","incid", "stdinc","N", "Sensit" ,"Specif" ,"tpev", "tinc")))

## Calculating lambda_t1 and lambda_t2 
lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])

## Storing estimated prevalence
result[m,1]<- pi1*lt1 +(1- pi1)*lt2

## Storing the variance
result[m,2]<- sqrt((pi1*lt1*(1-lt1)/f1+ (1-pi1)*lt2*(1-lt2)/f2+pi1*(1-pi1)*((lt1-lt2)^2))/length(Y))

## Storing values of incidence rate and sample size 
result[m,3]<- 0; result[m,4]<- 0; result[m,5]<- length(Y); 

## Storing sensitivity, specificity, true prevalence, true incidence
result[m,6]<- tb[1,1]/(tb[1,1]+tb[2,1]);result[m,7]<- tb[2,2]/(tb[1,2]+tb[2,2]); result[m,8]<- sum(temp[,2]==1)/nrow(temp); result[m,9]<-NA;

#rm(.Random.seed); 


### Analysis of rest of the 10 waves
while(m <= 10)
{
m<-m+1 # First wave

### Data Construction for the 
s<-ifelse(as.vector(new[,2])==2,0,1)

new<-cbind(new[,1][s==0], new[,2][s==0], new[,3][s==0])

### Dataset for the wave
temp<-rbind(new,comp)

### incidence rate
theta<-0.05; 

n1<-round(sum(temp[,2]==2)*theta); 

### new demented people
j<-1; 
for (i in 1:n1) {
  ind<- -1;
  while (ind<0)
  {if (temp[j,2]==2) {temp[j,1]<- rnorm(1,-2,4);  temp[j,2]<-1; ind<-1}; 
    j<- j+1;}
}    

### degradation of informant score
temp[,1]<-temp[,1]+((-1)^(temp[,2]+1))*runif(1,0,1.5)*rbinom(1,1,0.5)

### Screening test for data partition 
cl <- kmeans(temp[,1], 2, iter.max = 200, nstart = 10, algorithm = "MacQueen")

### Plot for visualization
plot(temp[,1], col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)

### Prediction for patients in screening test
Y<-as.integer(cl$cluster)

### Correction for predictions
center<- sort(cl$centers);
if (center[1]!= cl$centers[1])  {for (i in 1:length(Y)) {if (Y[i]==2) Y[i]=1 else Y[i]=2} }

### data preparation for GS test
sum(temp[,2]== Y) ##correct classification
temp[,3]<-Y
temp<- temp[order(temp[,3]), ]

### Measures for prediction of prevalence and incedence
n1<- sum(Y==2);n2<- sum(Y==1);s1<-rep(1,n2);
pi1<-n2/length(Y);
n<-round(n1*f2)

### simple random sampling
s=srswor(n, n1)
s<- c(s1,s)

### the sample for GS test
new <- temp[which(s==1),]
comp<- temp[which(s==0),]

### data frame formation
colnames(new) <- c("X","D","Y")
new<- as.data.frame(new)

### calculation of different quantities
tb<-table(new$Y, new$D)

lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])

## Storing the results
result[m,1]<- pi1*lt1 +(1- pi1)*lt2
result[m,2]<- sqrt((pi1*lt1*(1-lt1)/f1+ (1-pi1)*lt2*(1-lt2)/f2+pi1*(1-pi1)*((lt1-lt2)^2))/length(Y))
result[m,5]<- length(Y);
result[m,3]<- 1-result[m,5]*(1-result[m,1])/(result[m-1,5]*(1-result[m-1,1]))
result[m,4]<- sqrt(((1- result[m,1])^2)*(result[m-1,2]+ result[m,2]- result[m,2]*result[m-1,2])/((1- result[m-1,1])^2))
result[m,6]<- tb[1,1]/(tb[1,1]+tb[2,1]);result[m,7]<- tb[2,2]/(tb[1,2]+tb[2,2]);
result[m,8]<- sum(temp[,2]==1)/nrow(temp); result[m,9]<-theta;
}

## saving the table
write.table(result, "new_sim_info_degrade_fixed_incidence.csv")

#########################
####End of New Simulation 3
#########################




##########################
### Informants score degrades with time with variable incidence rate (model 6).
#########################

rm(list=ls())

set.seed(109)
### Assign  baseline
m<-1 #baseline

### Normal group
x<-as.vector(rnorm(900,2,1.5))  ##feature sample
D<-rep(2,length(x))             ##response 
x1<-cbind(x,D)                  ##normal group data

### Demented group
x<-as.vector(rnorm(100,-2,4))   ##feature sample
D <-rep(1, length(x))           ##response
x2<-cbind(x, D)                 ##demented group data

### Final data
z<-rbind(x1,x2)

### Screening test for data partition 
cl <- kmeans(z[,1], 2, iter.max = 200, nstart = 10, algorithm = "MacQueen")

### Plots for visualization
plot(z[,1], col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)

### Predicted partition
Y<-as.integer(cl$cluster)

center<- sort(cl$centers);
if (center[1]!= cl$centers[1])  {for (i in 1:length(Y)) {if (Y[i]==2) Y[i]=1 else Y[i]=2} }

### Checking the correct classification
sum(z[,2]== Y)

### COmbining original labels and predicted labels
temp<-cbind(z, Y)
temp<- temp[order(temp[,3]),]

### Calculating different quantities of two-phase sampling

## Assigning proportion of observations chosen from each group
f2<-0.1;f1<-1;                                

## Sample sizes and sample object
n1<- sum(Y==2);n2<- sum(Y==1);s1<-rep(1,n2);

## MLE of demented group
pi1<-n2/length(Y);

## Sample size from non-demented by screening test
n<-round(n1*f2)

### simple random sampling from non-demented group
s=srswor(n, n1)
s<- c(s1,s)

### Final sample for GS test
new <- temp[which(s==1),]
comp<- temp[which(s==0),]

### Data for GS test
colnames(new) <- c("X","D","Y")
new<- as.data.frame(new)
tb<-table(new$Y, new$D)

### Storing results for baseline
result <- matrix(NA, nrow = 11, ncol=9, byrow=TRUE, dimnames = list(c("baseline", "wave1", "wave2", "wave3","wave4", "wave5", "wave6","wave7", "wave8", "wave9", "wave10"), c("prev", "stdpre","incid", "stdinc","N", "Sensit" ,"Specif" ,"tpev", "tinc")))

## Calculating lambda_t1 and lambda_t2 
lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])

## Storing estimated prevalence
result[m,1]<- pi1*lt1 +(1- pi1)*lt2

## Storing the variance
result[m,2]<- sqrt((pi1*lt1*(1-lt1)/f1+ (1-pi1)*lt2*(1-lt2)/f2+pi1*(1-pi1)*((lt1-lt2)^2))/length(Y))

## Storing values of incidence rate and sample size 
result[m,3]<- 0; result[m,4]<- 0; result[m,5]<- length(Y); 

## Storing sensitivity, specificity, true prevalence, true incidence
result[m,6]<- tb[1,1]/(tb[1,1]+tb[2,1]);result[m,7]<- tb[2,2]/(tb[1,2]+tb[2,2]); result[m,8]<- sum(temp[,2]==1)/nrow(temp); result[m,9]<-NA;

#rm(.Random.seed);

while(m<=10){

m<-m+1 # First wave

### Dataset for the wave
s<-ifelse(as.vector(new[,2])==2,0,1)

new<-cbind(new[,1][s==0], new[,2][s==0], new[,3][s==0])

### Dataset for the wave
temp<-rbind(new,comp)

### New patients in Demented group
theta<-runif(1,0.04,0.06); # variable incidence rate
  
n1<-round(sum(temp[,2]==2)*theta);

j<-1;
for (i in 1:n1) {
  ind<- -1;
  while (ind<0)
  {if (temp[j,2]==2) {temp[j,1]<- rnorm(1,-2,4);  temp[j,2]<-1; ind<-1};
    j<- j+1;}
}    #newly demented people

### degradation of information score
temp[,1]<-temp[,1]+((-1)^(temp[,2]+1))*runif(1,0,1.5)*rbinom(1,1,0.5)

### Screening test for data partition 
cl <- kmeans(temp[,1], 2, iter.max = 200, nstart = 10)

### Plot for visualization
plot(temp[,1], col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)

### Prediction for patients in screening test
Y<-as.integer(cl$cluster)

### Correction for predictions
center<- sort(cl$centers);
if (center[1]!= cl$centers[1])  {for (i in 1:length(Y)) {if (Y[i]==2) Y[i]=1 else Y[i]=2} }

### data preparation for GS test
sum(temp[,2]== Y) ##correct classification
temp[,3]<-Y
temp<- temp[order(temp[,3]), ]

### Measures for prediction of prevalence and incedence
n1<- sum(Y==2);n2<- sum(Y==1);s1<-rep(1,n2);
pi1<-n2/length(Y);
n<-round(n1*f2)

## simple random sampling
s=srswor(n, n1)
s<- c(s1,s)

### the sample for GS test
new<-cbind(temp[,1][s==1], temp[,2][s==1], temp[,3][s==1])
comp<- cbind(temp[,1][s==0], temp[,2][s==0], temp[,3][s==0])

### data of GS test
colnames(new) <- c("X","D","Y")
new<- as.data.frame(new)
tb<-table(new$Y, new$D)

lt1<- tb[1,1]/(tb[1,1]+tb[1,2])
lt2<- tb[2,1]/ (tb[2,1]+tb[2,2])

### Storing the result
result[m,1]<- pi1*lt1 +(1- pi1)*lt2
result[m,2]<- ((pi1*lt1*(1-lt1)/f1+ (1-pi1)*lt2*(1-lt2)/f2+pi1*(1-pi1)*((lt1-lt2)^2))/length(Y))
result[m,5]<- length(Y);
result[m,3]<- 1-result[m,5]*(1-result[m,1])/(result[m-1,5]*(1-result[m-1,1]))
result[m,4]<- (((1- result[m,1])^2)*(result[m-1,2]+ result[m,2]- result[m,2]*result[m-1,2])/((1- result[m-1,1])^2))
result[m,6]<- tb[1,1]/(tb[1,1]+tb[2,1]);result[m,7]<- tb[2,2]/(tb[1,2]+tb[2,2]);
result[m,8]<- sum(temp[,2]==1)/nrow(temp); result[m,9]<-theta;

}

## saving the simulation results
write.table(result, "new_sim_info_degrade_variable_incidence.csv")
