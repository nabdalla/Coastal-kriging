library(spBayes)
library(classInt)
library(RColorBrewer)
library(MBA)
library(sp)
library(geoR)
library(fields)
library(maptools)
library(rgdal)
library(lattice)
library(geosphere)
library(nlme)
library(LearnBayes)
library(maps)
library(ggplot2)
library(plyr)
library(RgoogleMaps)
library(MASS)
library(grid)
library(raster)
library(rgeos)
library(plotrix)
library(mvtnorm)
library(mapproj)

#Data of GuLF STUDY
data<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project1/data 2.csv")
dim(data)
THC<-subset(data, data$Reported_Analyte=="Total Hydrocarbons")
plot(density(log(THC$Result2)))
hist(log(THC$Result2))
summary(log(THC$Result2))
#map plot
dev.off()

sortthc<-THC[order(THC$device_longitude),]
#Sinusoidal projection
mapping<-mapproject(sortthc$device_longitude, sortthc$device_latitude, proj= "sinusoidal", parameters=NULL, orientation=NULL)
sortthc$x<-mapping$x^3
sortthc$y<-mapping$y
plot(sortthc$x,sortthc$y)
head(sortthc)
dim(sortthc)
#log transform for normality
sortthc$log<-log(sortthc$Result2)
#Radius of earth
radius<-6371

#Total distance of the coast in km using the approximation
sum(sortthc$tthc)*radius
#Total Euclidean distance of the coast in km 
radius*sqrt((max(sortthc$x)-min(sortthc$x))^2+(max(sortthc$y)-min(sortthc$y))^2)
longthc<-sortthc$device_longitude
latthc<-sortthc$device_latitude
#Calculation of the parameter corresponding to line segment approximations to the coast
sortthc$tthc<-NULL
for ( j in 2:length(sortthc$log)){
  sortthc$tthc[j-1]<-sqrt((sortthc$x[j]-sortthc$x[j-1])^2+(sortthc$y[j]-sortthc$y[j-1])^2)
}
#Remove duplicated locations
sortthc$thc<-NULL
for ( j in 2:length(sortthc$log)){
  sortthc$thc[j-1]<-sqrt((longthc[j]-longthc[j-1])^2+(longthc[j]-longthc[j-1])^2)
}
thc<-c(0,sortthc$thc)
sortthc2<-sortthc
sortthc2<-sortthc2[!duplicated(sortthc2$thc),]
dim(sortthc2)
sortthc2$tthc2<-cumsum(sortthc2$tthc)
#Distance matrix in km
distthc<-radius*abs(outer(sortthc2$tthc2,sortthc2$tthc2,"-"))
#sortthc2<-sortthc2[which(radius*sortthc2$x > -1),]
summary(sortthc2$x)
dim(sortthc2)
#MLE of the parameters to use as fixed values in models 2a and 2b
euc<-likfit(coords=cbind(radius*sortthc2$x, radius*sortthc2$y), data=sortthc2$log, ini.cov.pars=c(1,1))
sigmathc<-euc$cov.pars[1]
phithc<-euc$cov.pars[2]
tauthc<-euc$nugget
alphathc<-tauthc/sigmathc
#Model 2b
x0<-rep(1,length(sortthc2$log))
n.samples=5000
p=1
beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)
sigma.sq.prior.shape <- 2
sigma.sq.prior.rate <-2
#trainthc=c(seq(1,29,1),seq(43,length(sortthc2$log),1)) 
trainthc=c(seq(21,length(sortthc2$log),1))
testthc<-seq(1,20,1)
n1=length(trainthc)
n2<-length(testthc)
#trainthc<-seq(1,38,1)
exact1thc <- bayesGeostatExact(sortthc2$log[trainthc]~1, n.samples=n.samples,
                               beta.prior.mean=beta.prior.mean,
                               beta.prior.precision=beta.prior.precision,
                               coords=cbind(radius*sortthc2$tthc2[trainthc], rep(0,length(sortthc2$tthc2[trainthc]))), phi=phithc, alpha=alphathc,
                               sigma.sq.prior.shape=sigma.sq.prior.shape,
                               sigma.sq.prior.rate=sigma.sq.prior.rate)
print(summary(exact1thc$p.samples))
beta0thc.mean<-mean(exact1thc$p.samples[,1])
sigmathc.post.mean<-mean(exact1thc$p.samples[,2])
tauthc.mean<-mean(exact1thc$p.samples[,3])
beta0thc<-(exact1thc$p.samples[,1])
sigmathc.post<-(exact1thc$p.samples[,2])
tauthc<-(exact1thc$p.samples[,3])
#Prediction using test sample
tpthc<-c(sortthc2$tthc2[trainthc],sortthc2$tthc2[testthc])
length(tpthc)
dist2thc<-abs(outer(radius*tpthc,radius*tpthc,"-"))
v2thc<-vv.p2thc<-muhat2thc<-yhatthc<-vector("list")
for(i in 1:n.samples){
  v2thc[[i]]<-sigmathc.post[i]*exp(-phithc*dist2thc)+tauthc[i]*diag(n2+n1)
  vv.p2thc[[i]]<- v2thc[[i]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v2thc[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2thc[[i]][1:n1,1:n1]))%*%(t(v2thc[[i]][(n1+1):(n2+n1), 1:n1]))
  muhat2thc[[i]]<-t(beta0thc[i]%*%t(rep(1, (n2))))+(v2thc[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2thc[[i]][1:n1,1:n1]))%*%t((t(sortthc2$log[trainthc])-beta0thc[i]%*%t(x0[1:n1])))
  set.seed(i)
  yhatthc[[i]]<-rmvnorm(1,muhat2thc[[i]],round(vv.p2thc[[i]],4))
  
}
summary(unlist(muhat2thc), quantile)
mean(unlist(muhat2thc))
yhatthc.mat<-matrix(unlist(yhatthc), nrow=n.samples, ncol=n2, byrow=T)
predthc<-NULL
for (i in 1:n2){
  predthc[i]<-mean(yhatthc.mat[,i])
}
#DIC calculation
pdfthc<-NULL
for(i in 1:n.samples){
  pdfthc[i]<-dmvnorm(sortthc2$log[testthc],muhat2thc[[i]],round(vv.p2thc[[i]],4))
}
v2thcmean<-sigmathc.post.mean*exp(-phithc*dist2thc)+tauthc.mean*diag(n1+n2)
vv.p2thcmean<- v2thcmean[(n1+1):(n2+n1), (n1+1):(n2+n1)]-((v2thcmean[(n1+1):(n2+n1), 1:n1])%*%
               chol2inv(chol(v2thcmean[1:n1,1:n1])))%*%(t(v2thcmean[(n1+1):(n2+n1), 1:n1]))
muhat2thcmean<-t(beta0thc.mean%*%t(rep(1, (n2))))+((v2thcmean[(n1+1):(n2+n1), 1:n1])%*%
               chol2inv(chol(v2thcmean[1:n1,1:n1])))%*%t((t(sortthc2$log[trainthc])-beta0thc.mean%*%t(x0[1:n1])))

D.barthc <- -2*log(mean(pdfthc))
D.hatthc <- -2*log(dmvnorm(sortthc2$log[testthc],muhat2thcmean,vv.p2thcmean))
pDthc <- D.barthc - D.hatthc
DICthc <- pDthc+D.barthc
#Mean square error calculation
msethc<-sum((exp(sortthc2$log[testthc])-exp(predthc))^2)/length(sortthc2$log[testthc])

#using Euclidean distance method

exactthc_euc <- bayesGeostatExact(sortthc2$log[trainthc]~1, n.samples=n.samples,
                                 beta.prior.mean=beta.prior.mean,
                                 beta.prior.precision=beta.prior.precision,
                                 coords=cbind(radius*sortthc2$x[trainthc], radius*sortthc2$y[trainthc]), phi=phithc, alpha=alphathc,
                                 sigma.sq.prior.shape=sigma.sq.prior.shape,
                                 sigma.sq.prior.rate=sigma.sq.prior.rate)
print(summary(exactthc_euc$p.samples))
beta0thc_euc<-(exactthc_euc$p.samples[,1])
sigmathc.post_euc<-(exactthc_euc$p.samples[,2])
tauthc_euc<-(exactthc_euc$p.samples[,3])
beta0thc_euc.mean<-mean(exactthc_euc$p.samples[,1])
sigmathc.post_euc.mean<-mean(exactthc_euc$p.samples[,2])
tauthc_euc.mean<-mean(exactthc_euc$p.samples[,3])
listthc<-data.frame(latitude=radius*c(sortthc2$x[trainthc],sortthc2$x[testthc]),
                    longitude=radius*c(sortthc2$y[trainthc], sortthc2$y[testthc]))  
matthc <- as.matrix(dist(listthc[,c('longitude','latitude')],upper=TRUE, diag=TRUE))
v2thc_euc<-vv.p2thc_euc<-muhat2thc_euc<-yhatthc_euc<-vector("list")
for(i in 1:n.samples){
  v2thc_euc[[i]]<-sigmathc.post_euc[i]*exp(-phithc*matthc)+tauthc_euc[i]*diag(n2+n1)
  vv.p2thc_euc[[i]]<- v2thc_euc[[i]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v2thc_euc[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2thc_euc[[i]][1:n1,1:n1]))%*%(t(v2thc_euc[[i]][(n1+1):(n2+n1), 1:n1]))
  muhat2thc_euc[[i]]<-t(beta0thc_euc[i]%*%t(rep(1, (n2))))+(v2thc_euc[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2thc_euc[[i]][1:n1,1:n1]))%*%t((t(sortthc2$log[trainthc])-beta0thc_euc[i]%*%t(x0[1:n1])))
  set.seed(i)
  yhatthc_euc[[i]]<-rmvnorm(1,muhat2thc_euc[[i]],vv.p2thc_euc[[i]])
}
summary(unlist(muhat2thc_euc), quantile)
yhatthc.mat_euc<-matrix(unlist(yhatthc_euc), nrow=n.samples, ncol=n2, byrow=T)
predthc_euc<-NULL
for (i in 1:n2){
  predthc_euc[i]<-mean(yhatthc.mat_euc[,i])
}

pdfthc_euc<-NULL
for(i in 1:n.samples){
  pdfthc_euc[i]<-dmvnorm(sortthc2$log[testthc],muhat2thc_euc[[i]],(vv.p2thc_euc[[i]]))
}
v2thcmean_euc<-sigmathc.post_euc.mean*exp(-phithc*matthc)+tauthc_euc.mean*diag(n1+n2)
vv.p2thcmean_euc<- v2thcmean_euc[(n1+1):(n2+n1), (n1+1):(n2+n1)]-((v2thcmean_euc[(n1+1):(n2+n1), 1:n1])%*%
                  chol2inv(chol(v2thcmean_euc[1:n1,1:n1])))%*%(t(v2thcmean_euc[(n1+1):(n2+n1), 1:n1]))
muhat2thcmean_euc<-t(beta0thc_euc.mean%*%t(rep(1, (n2))))+((v2thcmean_euc[(n1+1):(n2+n1), 1:n1])%*%
                  chol2inv(chol(v2thcmean_euc[1:n1,1:n1])))%*%t((t(sortthc2$log[trainthc])-beta0thc_euc.mean%*%t(x0[1:n1])))

D.barthc_euc <- -2*log(mean(pdfthc_euc))
D.hatthc_euc <- -2*log(dmvnorm(sortthc2$log[testthc],muhat2thcmean_euc,vv.p2thcmean_euc))
pDthc_euc <- D.barthc_euc - D.hatthc_euc
DICthc_euc <- pDthc_euc+D.barthc_euc
msethc_euc<-sum((exp(sortthc2$log[testthc])-exp(predthc_euc))^2)/length(sortthc2$log[testthc])

#Model 1b
priors.2 <- list("beta.flat",
                 "phi.Unif"=c(0.1, 30), "sigma.sq.IG"=c(2,2),
                 "tau.sq.IG"=c(2, 2))
starting2 <- list("phi"=phithc, "sigma.sq"=sigmathc, "tau.sq"=tauthc)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
sp.exactthc <-spLM(sortthc2$log[trainthc]~1,
                   coords=cbind(radius*sortthc2$tthc2[trainthc], rep(0,n1)), tuning=tuning, starting=starting2, n.samples=5000,
                   priors=priors.2,cov.model="exponential")
burn.in <- 0.5*n.samples
m.11thc <- spRecover(sp.exactthc , start=burn.in, verbose=FALSE)
spDiag(m.11thc)
#Prediction using test data
m.predthc <- spPredict(m.11thc, pred.covars=matrix(rep(1,n2),n2,1), 
                       pred.coords=cbind(radius*sortthc2$tthc2[testthc], rep(0,(n2))),
                       start=0.5*n.samples)
y.hatthc <- as.numeric(apply(m.predthc$p.y.predictive.samples, 1, mean))
beta0spthc<-mean(m.11thc$p.beta.recover.samples[,1])
sigmaspthc<-mean(m.11thc$p.theta.recover.samples[,1])
tauspthc<-mean(m.11thc$p.theta.recover.samples[,2])
phispthc<-mean(m.11thc$p.theta.recover.samples[,3])
beta0spthc.post<-(m.11thc$p.beta.recover.samples[,1])
sigmaspthc.post<-(m.11thc$p.theta.recover.samples[,1])
tauspthc.post<-(m.11thc$p.theta.recover.samples[,2])
phispthc.post<-(m.11thc$p.theta.recover.samples[,3])
round(summary(m.11thc$p.theta.recover.samples)$quantiles,3)
round(summary(m.11thc$p.beta.recover.samples)$quantiles,3)
#For DIC calculatiosn
v2thcsp<-vv.p2thcsp<-vector("list")
pdfthcsp<-NULL
for(i in 1:(n.samples*0.5)){
  v2thcsp[[i]]<-m.11thc$p.theta.recover.samples[i,1]*exp(-m.11thc$p.theta.recover.samples[i,3]*dist2thc)+m.11thc$p.theta.recover.samples[i,2]*diag(n2+n1)
  vv.p2thcsp[[i]]<- v2thcsp[[i]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v2thcsp[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2thcsp[[i]][1:n1,1:n1]))%*%(t(v2thcsp[[i]][(n1+1):(n2+n1), 1:n1]))
  pdfthcsp[i]<-dmvnorm(sortthc2$log[testthc],m.predthc$p.y.predictive.samples[,i],vv.p2thcsp[[i]])
}
v2thcmeansp<-sigmaspthc*exp(-phispthc*dist2thc)+tauspthc*diag(n1+n2)
vv.p2thcmeansp<- v2thcmeansp[(n1+1):(n2+n1), (n1+1):(n2+n1)]-((v2thcmeansp[(n1+1):(n2+n1), 1:n1])%*%
                 chol2inv(chol(v2thcmeansp[1:n1,1:n1])))%*%(t(v2thcmeansp[(n1+1):(n2+n1), 1:n1]))
muhat2thcmeansp<-t(beta0spthc%*%t(rep(1, (n2))))+((v2thcmeansp[(n1+1):(n2+n1), 1:n1])%*%
                chol2inv(chol(v2thcmeansp[1:n1,1:n1])))%*%t((t(sortthc2$log[trainthc])-beta0spthc%*%t(x0[1:n1])))

D.barthcsp <- -2*log(mean(pdfthcsp))
D.hatthcsp <- -2*log(dmvnorm(sortthc2$log[testthc],muhat2thcmeansp,vv.p2thcmeansp))
pDthcsp <- D.barthcsp - D.hatthcsp
DICthcsp <- pDthcsp+D.barthcsp
#MSE calculations
msethc2<-sum((exp(sortthc2$log[testthc])-exp(y.hatthc))^2)/length(sortthc2$log[testthc])
#Euclidean distance method
sp.exactthc_euc <-spLM(sortthc2$log[trainthc]~1,
                      coords=cbind(radius*sortthc2$x[trainthc], radius*sortthc2$y[trainthc]), tuning=tuning, starting=starting2, n.samples=5000,
                      priors=priors.2,cov.model="exponential")
mthc_euc <- spRecover(sp.exactthc_euc , start=burn.in, verbose=FALSE)
m.predthc_euc <- spPredict(mthc_euc, pred.covars=matrix(rep(1,n2),n2,1), 
                          pred.coords=cbind(sortthc2$x[testthc], sortthc2$y[testthc]),
                          start=0.5*n.samples)
y.hatthc_euc <- as.numeric(apply(m.predthc_euc$p.y.predictive.samples, 1, mean))
beta0spthc_euc<-mean(mthc_euc$p.beta.recover.samples[,1])
sigmaspthc_euc<-mean(mthc_euc$p.theta.recover.samples[,1])
tauspthc_euc<-mean(mthc_euc$p.theta.recover.samples[,2])
phispthc_euc<-mean(mthc_euc$p.theta.recover.samples[,3])
beta0spthc_euc.post<-(mthc_euc$p.beta.recover.samples[,1])
sigmaspthc_euc.post<-(mthc_euc$p.theta.recover.samples[,1])
tauspthc_euc.post<-(mthc_euc$p.theta.recover.samples[,2])
phispthc_euc.post<-(mthc_euc$p.theta.recover.samples[,3])
round(summary(mthc_euc$p.theta.recover.samples)$quantiles,3)
round(summary(mthc_euc$p.beta.recover.samples)$quantiles,3)

v2thcsp_euc<-vv.p2thcsp_euc<-vector("list")
pdfthcsp_euc<-NULL
for(i in 1:(n.samples*0.5)){
  v2thcsp_euc[[i]]<-mthc_euc$p.theta.recover.samples[i,1]*exp(-mthc_euc$p.theta.recover.samples[i,3]*matthc)+mthc_euc$p.theta.recover.samples[i,2]*diag(n2+n1)
  vv.p2thcsp_euc[[i]]<- v2thcsp_euc[[i]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v2thcsp_euc[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2thcsp_euc[[i]][1:n1,1:n1]))%*%(t(v2thcsp_euc[[i]][(n1+1):(n2+n1), 1:n1]))
  pdfthcsp_euc[i]<-dmvnorm(sortthc2$log[testthc],m.predthc_euc$p.y.predictive.samples[,i],vv.p2thcsp_euc[[i]])
}
v2thcmeansp_euc<-sigmaspthc_euc*exp(-phispthc_euc*matthc)+tauspthc_euc*diag(n1+n2)
vv.p2thcmeansp_euc<- v2thcmeansp_euc[(n1+1):(n2+n1), (n1+1):(n2+n1)]-((v2thcmeansp_euc[(n1+1):(n2+n1), 1:n1])%*%
                    chol2inv(chol(v2thcmeansp_euc[1:n1,1:n1])))%*%(t(v2thcmeansp_euc[(n1+1):(n2+n1), 1:n1]))
muhat2thcmeansp_euc<-t(beta0spthc_euc%*%t(rep(1, (n2))))+((v2thcmeansp_euc[(n1+1):(n2+n1), 1:n1])%*%
                     chol2inv(chol(v2thcmeansp_euc[1:n1,1:n1])))%*%t((t(sortthc2$log[trainthc])-beta0spthc_euc%*%t(x0[1:n1])))

D.barthcsp_euc <- -2*log(mean(pdfthcsp_euc))
D.hatthcsp_euc <- -2*log(dmvnorm(sortthc2$log[testthc],muhat2thcmeansp_euc,vv.p2thcmeansp_euc))
pDthcsp_euc <- D.barthcsp_euc - D.hatthcsp_euc
DICthcsp_euc <- pDthcsp_euc+D.barthcsp_euc

sum((exp(sortthc2$log[testthc])-exp(y.hatthc_euc))^2)/length(sortthc2$log[testthc])

#G+P=D
vd1s<-vd1b<-vd2s<-vd2b<-vector("list")
yhatd1s<-yhatd1b<-yhatd2s<-yhatd2b<-matrix(0,nrow=(0.5*n.samples), ncol=n1)
for(j in 1:(0.5*n.samples)){
  vd1s[[j]]<-sigmaspthc_euc.post[j]*exp(-phispthc_euc.post[j]*matthc[1:n1,1:n1])+tauspthc_euc.post[j]*diag(n1)
  yhatd1s[j,1:n1]<-rmvnorm(1,rep(beta0spthc_euc.post[j],n1),vd1s[[j]])
  vd1b[[j]]<-sigmaspthc.post[j]*exp(-phispthc.post[j]*dist2thc[1:n1,1:n1])+tauspthc.post[j]*diag(n1)
  yhatd1b[j,1:n1]<-rmvnorm(1,rep(beta0spthc.post[j],n1),vd1b[[j]])
  vd2s[[j]]<-sigmathc.post_euc[j]*exp(-phithc*matthc[1:n1,1:n1])+tauthc_euc[j]*diag(n1)
  yhatd2s[j,1:n1]<-rmvnorm(1,rep(beta0thc_euc[j],n1),vd2s[[j]])
  vd2b[[j]]<-sigmathc.post[j]*exp(-phithc*dist2thc[1:n1,1:n1])+tauthc[j]*diag(n1)
  yhatd2b[j,1:n1]<-rmvnorm(1,rep(beta0thc[j],n1),vd2b[[j]])
}

murepd1s<-apply(yhatd1s,2 , mean)
murepd1b<-apply(yhatd1b,2 , mean)
murepd2s<-apply(yhatd2s,2 , mean)
murepd2b<-apply(yhatd2b,2 , mean)
Pd1s<-apply(yhatd1s,2 , var)
Pd1b<-apply(yhatd1b,2 , var)
Pd2s<-apply(yhatd2s,2 , var)
Pd2b<-apply(yhatd2b,2 , var)
Gd1s<-sum((sortthc2$log[trainthc]-murepd1s)^2)
Gd1b<-sum((sortthc2$log[trainthc]-murepd1b)^2)
Gd2s<-sum((sortthc2$log[trainthc]-murepd2s)^2)
Gd2b<-sum((sortthc2$log[trainthc]-murepd2b)^2)
Dd1s<-sum(Pd1s)+Gd1s
Dd1b<-sum(Pd1b)+Gd1b
Dd2s<-sum(Pd2s)+Gd2s
Dd2b<-sum(Pd2b)+Gd2b
Dd1s
Dd1b
Dd2s
Dd2b
#plot THC

quartz()
par(mfrow=c(1,2))
plot(exp(y.hatthc), exp(sortthc2$log[testthc]),xlab="Predicted", ylab="True",
     main="Model 1b",
     cex.main=0.9,cex=0.9)
abline(coef=c(0,1))
plot(exp(predthc), exp(sortthc2$log[testthc]),xlab="Predicted", ylab="True",
     main="Model 2b",
     cex.main=0.9,cex=0.9)
abline(coef=c(0,1))

col.pal <- col.br(5)
quantthc2b<-classIntervals(predthc, style="quantile")
quant.colthc2b<-findColours(quantthc2b, col.pal)
quantthc1b<-classIntervals(y.hatthc, style="quantile")
quant.colthc1b<-findColours(quantthc1b, col.pal)
quantthc<-classIntervals(sortthc2$log[testthc], style="quantile")
quant.colthc<-findColours(quantthc, col.pal)
quartz()
par(mfrow=c(1,3))
plot(radius*cbind(sortthc2$x[testthc], sortthc4$y[testthc]), col=quant.colthc, xlab="Easting", ylab="Northing",main="True THC values")
plot(radius*cbind(sortthc2$x[testthc], sortthc4$y[testthc]), col=quant.colthc1b,xlab="Easting", ylab="Northing", main="Predicted values model 1b")
plot(radius*cbind(sortthc2$x[testthc], sortthc4$y[testthc]), col=quant.colthc2b,xlab="Easting", ylab="Northing", main="Predicted values model 2b")

#Google map with predictions
newx<-approx((sortthc$device_longitude),(sortthc$device_latitude),n=100)
plot(newx$x,newx$y)
mappingnew<-mapproject(newx$x, newx$y, proj= "sinusoidal", parameters=NULL, orientation=NULL)

#Calculation of the parameter corresponding to line segment approximations to the coast
newtthc<-NULL
for ( j in 2:100){
  newtthc[j-1]<-sqrt((mappingnew$x[j]-mappingnew$x[j-1])^2+(mappingnew$y[j]-mappingnew$y[j-1])^2)
}

newtthc<-c(0,newtthc)
newtthc2<-cumsum(newtthc)
m.predthc <- spPredict(m.11thc, pred.covars, 
                       pred.coords=cbind(radius*newtthc2, rep(0,(100))),
                       start=0.5*n.samples)
y.hatthc.new <- as.numeric(apply(m.predthc$p.y.predictive.samples, 1, mean))

col.br=colorRampPalette(c("blue", "cyan", "yellow", "red"))
col.pal <- col.br(5)
quantthc<-classIntervals(exp(y.hatthc.new) , style="quantile")
quant.colthc<-findColours(quantthc, col.pal)
quartz()
par(mfrow=c(2,1), oma=c(3,3,2,2))
#layout(rbind(1,2,3),heights=c(8,7,1))
MyMapthc <- GetMap.bbox(lonR = c(-89.38,-89.35), latR = c(30,30.35),  
                        size=c(640,640), maptype = "hybrid")
PlotOnStaticMap(MyMapthc)
convert_pointsthcnew <- LatLon2XY.centered(MyMapthc, newx$y, newx$x)

points(convert_pointsthcnew$newX, convert_pointsthcnew$newY, col = quant.colthc, pch=19,cex=0.3)
plot(radius*mappingnew$x,radius*mappingnew$y, col=quant.colthc, cex=0.5,ylab="Northing", xlab="Easting")
     mtext(text="Northing",side=2,line=2)
     mtext(text="Easting",side=1,line=2)
legend("topleft",fill=attr(quant.colthc,"palette"),legend=
         names(attr(quant.colthc,"table")),bg="white", cex=1, bty="n")

quartz()
par(mfrow=c(2,1), oma=c(2,2,2,2))

data(wrld_simpl)
x <- list(x=-90:-89, y = 30:31, z = outer(1:2, 1:2, "+"))
image(x=c(-89.5,-89), y = c(30,30.5), z = outer(1:1, 1:1, "+"), 
      xlab = "lon", ylab = "lat")
map("state", add=TRUE)
text(x=-89.4, y=30.4, "MS", cex=0.9)
outline <- map("usa", plot=FALSE) # returns a list of x/y coords
xrange <- range(outline$x, na.rm=TRUE) # get bounding box
yrange <- range(outline$y, na.rm=TRUE)
xbox <- xrange + c(-2, 2)
ybox <- yrange + c(-2, 2)
# create the grid path in the current device
polypath(c(outline$x, NA, c(xbox, rev(xbox))),
         c(outline$y, NA, rep(ybox, each=2)),
         col="light blue", rule="evenodd")
points(THC$device_longitude,THC$device_latitude+0.03, col="red", cex=0.5, xlab="Longitude", ylab="Latitude")
#points(convert_pointsthcnew$newX, convert_pointsthcnew$newY, col = quant.colthc, pch=19,cex=0.3)
plot(radius*mappingnew$x,radius*mappingnew$y, col=quant.colthc, cex=0.5, xlab="lon", ylab="lat")
legend("topleft",fill=attr(quant.colthc,"palette"),legend=
         names(attr(quant.colthc,"table")),bg="white", cex=1, bty="n")

