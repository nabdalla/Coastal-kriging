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
THC<-subset(data, data$Reported_Analyte=="Total Hydrocarbons")
sortthc<-THC[order(THC$device_longitude),]
#Sinusoidal projection
longthc<-sortthc$device_longitude
latthc<-sortthc$device_latitude
mapping<-mapproject(sortthc$device_longitude, sortthc$device_latitude, proj= "sinusoidal", parameters=NULL, orientation=NULL)
sortthc$x<-mapping$x
sortthc$y<-mapping$y
#log transform for normality
sortthc$log<-log(sortthc$Result2)
#Radius of earth
radius<-6371
line.seg<-function(x,y,data){
  xy<-cbind(x,y)
  xy<-xy[order(x,-y),]
  data<-data[order(x,-y)]
  t<-rep(0,length(x))
  for ( j in 2:length(x)){
    t[j]<-sqrt((xy[j,1]-xy[j-1,1])^2+(xy[j,2]-xy[j-1,2])^2)
  }
  t<-cumsum(t)
  return(list(t=t,data=data))
}
tcalc<-line.seg(sortthc$x,sortthc$y,sortthc$log)
sortthc$tthc<-tcalc$t
sortthc$log<-tcalc$data
#Remove duplicated locations
sortthc$thc<-NULL
for ( j in 2:length(sortthc$log)){
  sortthc$thc[j-1]<-sqrt((longthc[j]-longthc[j-1])^2+(longthc[j]-longthc[j-1])^2)
}
thc<-c(0,sortthc$thc)
sortthc2<-sortthc
sortthc2<-sortthc2[!duplicated(sortthc2$thc),]
#MLE of the parameters to use as fixed values in models 2a and 2b
euc<-likfit(coords=cbind(radius*sortthc2$x, radius*sortthc2$y), data=sortthc2$log, ini.cov.pars=c(1,1))
sigmathc<-euc$cov.pars[1]
phithc<-euc$cov.pars[2]
tauthc<-euc$nugget
alphathc<-tauthc/sigmathc
x0<-rep(1,length(sortthc2$log))

trainthc=c(seq(21,length(sortthc2$log),1))
testthc<-seq(1,20,1)
n1=length(trainthc)
n2<-length(testthc)
#SpLm
priors1 <- list("beta.flat",
                 "phi.Unif"=c(0.1, 30), "sigma.sq.IG"=c(2,2),
                 "tau.sq.IG"=c(2, 2))
starting1 <- list("phi"=phithc, "sigma.sq"=sigmathc, "tau.sq"=tauthc)
tuning1 <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
Geostat.full<-function(x,y,data,n.samples){
  sp <-spLM(data~1, coords=cbind(radius*x, radius*y), tuning=tuning1, starting=starting1, n.samples=5000,
                     priors=priors1,cov.model="exponential")
  burn.in <- 0.5*n.samples
  m <- spRecover(sp , start=burn.in, verbose=FALSE)
  beta0<-(m$p.beta.recover.samples[,1])
  sigma<-(m$p.theta.recover.samples[,1])
  tau<-(m$p.theta.recover.samples[,2])
  phi<-(m$p.theta.recover.samples[,3])
  return(list(m=m,beta0=beta0, sigma=sigma, tau=tau,phi=phi))
}
model1b<-Geostat.full(x=radius*sortthc2$tthc[trainthc],y=rep(0,n1),data=sortthc2$log[trainthc],
                      n.samples=5000)
beta01b<-model1b$beta0
sigma1b<-model1b$sigma
tau1b<-model1b$tau
phi1b<-model1b$phi
m1b<-model1b$m
modelfsk<-Geostat.full(x=radius*sortthc2$x[trainthc],y=radius*sortthc2$y[trainthc],data=sortthc2$log[trainthc],5000)
beta0fsk<-modelfsk$beta0
sigmafsk<-modelfsk$sigma
taufsk<-modelfsk$tau
phifsk<-modelfsk$phi
mfsk<-modelfsk$m
#geostat.exact
Geostat.exact<-function(x,y,data,n.samples,beta.prior.mean,beta.prior.precision,
                        sigma.sq.prior.shape, sigma.sq.prior.rate,phi,alpha ){
  exact <- bayesGeostatExact(data~1, n.samples=n.samples,
                                 beta.prior.mean=beta.prior.mean,
                                 beta.prior.precision=beta.prior.precision,
                                 coords=cbind(x, y), phi=phi, alpha=alpha,
                                 sigma.sq.prior.shape=sigma.sq.prior.shape,
                                 sigma.sq.prior.rate=sigma.sq.prior.rate)
  beta0<-(exact$p.samples[,1])
  sigma<-(exact$p.samples[,2])
  tau<-(exact$p.samples[,3])
  return(list(beta0=beta0, sigma=sigma, tau=tau))
}
model2b<-Geostat.exact(x=radius*sortthc2$tthc[trainthc],y=rep(0,n1),data=sortthc2$log[trainthc],5000,0,0,
              2,2,phi=phithc,alpha=alphathc)
beta02b<-model2b$beta0
sigma2b<-model2b$sigma
tau2b<-model2b$tau
modelssk<-Geostat.exact(x=radius*sortthc2$x[trainthc],y=radius*sortthc2$y[trainthc],data=sortthc2$log[trainthc],5000,0,0,
                        2,2,phi=phithc,alpha=alphathc)
beta0ssk<-modelssk$beta0
sigmassk<-modelssk$sigma
taussk<-modelssk$tau
#Prediction and goodness of fit using test sample simplified model

goodnessoffit.simple<-function(datatrain, datatest,sigma,phi,tau,beta0, n.samples,
                               xtrain, xtest,ytrain,ytest,coastal,simple,m){
  if(coastal=="TRUE"){
    ttest=ttrain=NULL
    #x y of projected coordinates
      xytrain<-cbind(xtrain,ytrain)
      xytest<-cbind(xtest,ytest)
      xytrain<-xytrain[order(xtrain,-ytrain),]
      datatrain<-datatrain[order(xtrain,-ytrain)]
      xytest<-xytest[order(xtest,-ytest),]
      datatest[order(xtest,-ytest)]
      ttrain<-rep(0,length(xtrain))
      ttest<-rep(0,length(xtest))
      for ( j in 2:length(xtrain)){
        ttrain[j]<-sqrt((xytrain[j,1]-xytrain[j-1,1])^2+(xytrain[j,2]-xytrain[j-1,2])^2)
      }
      ttrain<-cumsum(ttrain)
      for ( j in 2:length(xtest)){
        ttest[j]<-sqrt((xytest[j,1]-xytest[j-1,1])^2+(xytest[j,2]-xytest[j-1,2])^2)
      }
      ttest<-cumsum(ttest)
      ttot<-c(ttrain,ttest)
      dist<-abs(outer(radius*ttot,radius*ttot,"-"))
  } else{
  listt<-data.frame(latitude=radius*c(xtrain,xtest),
                      longitude=radius*c(ytrain, ytest))  
  dist <- as.matrix(dist(listt[,c('longitude','latitude')],upper=TRUE, diag=TRUE))}
  v<-vv<-muhat<-yhat<-vd<-vector("list")
  pdf<-NULL
  n1<-length(ttrain)
  n2<-length(ttest)
  yhatd<-matrix(0,nrow=(0.5*n.samples), ncol=n1)
  beta0.mean<-mean(beta0)
  sigma.mean<-mean(sigma)
  tau.mean<-mean(tau)
  for(i in 1:(0.5*n.samples)){
    v[[i]]<-sigma[i]*exp(-phi[i]*dist)+tau[i]*diag(n2+n1)
    vv[[i]]<- v[[i]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v[[i]][(n1+1):(n2+n1), 1:n1])%*%
      chol2inv(chol(v[[i]][1:n1,1:n1]))%*%(t(v[[i]][(n1+1):(n2+n1), 1:n1]))
    muhat[[i]]<-t(beta0[i]%*%t(rep(1, (n2))))+(v[[i]][(n1+1):(n2+n1), 1:n1])%*%
      chol2inv(chol(v[[i]][1:n1,1:n1]))%*%t((t(datatrain)-beta0[i]*t(x0[1:n1])))
    set.seed(i)
    yhat[[i]]<-rmvnorm(1,muhat[[i]],round(vv[[i]],4))
    pdf[i]<-dmvnorm(datatest,muhat[[i]],round(vv[[i]],4))
    vd[[i]]<-sigma[i]*exp(-phi[i]*dist[1:n1,1:n1])+tau[i]*diag(n1)
    yhatd[i,1:n1]<-rmvnorm(1,rep(beta0[i],n1),vd[[i]])
  }
  if(simple=="TRUE"){
    yhat.mat<-matrix(unlist(yhat), nrow=(n.samples*0.5), ncol=n2, byrow=T)
    pred<-NULL
    for (i in 1:n2){
      pred[i]<-mean(yhat.mat[,i])
    }
  }else{
  m.pred <- spPredict(m, pred.covars=matrix(rep(1,n2),n2,1), 
                      pred.coords=cbind(radius*xtest, radius*ytest),
                      start=0.5*n.samples)
  pred <- as.numeric(apply(m.pred$p.y.predictive.samples, 1, mean))}
 
  v.mean<-sigma.mean*exp(-phi[i]*dist)+tau.mean*diag(n1+n2)
  vv.mean<- v.mean[(n1+1):(n2+n1), (n1+1):(n2+n1)]-((v.mean[(n1+1):(n2+n1), 1:n1])%*%
               chol2inv(chol(v.mean[1:n1,1:n1])))%*%(t(v.mean[(n1+1):(n2+n1), 1:n1]))
  muhat.mean<-t(beta0.mean%*%t(rep(1, (n2))))+((v.mean[(n1+1):(n2+n1), 1:n1])%*%
               chol2inv(chol(v.mean[1:n1,1:n1])))%*%t((t(datatrain)-beta0.mean%*%t(x0[1:n1])))
  D.bar <- -2*log(mean(pdf))
  D.hat <- -2*log(dmvnorm(datatest,muhat.mean,vv.mean))
  pD <- D.bar - D.hat
  DIC <- pD+D.bar
  mse<-sum((exp(datatest)-exp(pred))^2)/n2
  mud<-apply(yhatd,2 , mean)
  Pd<-apply(yhatd,2 , var)
  Gd<-sum((datatrain-mud)^2)
  Dd<-sum(Pd)+Gd
  return(list(pred=pred,DIC=DIC,MSE=mse,Dd=Dd))
}
gof2b<-goodnessoffit.simple(datatrain=sortthc2$log[trainthc], datatest=sortthc2$log[testthc],sigma=sigma2b,phi=rep(phithc,n.samples),tau=tau2b,
              beta0=beta02b,n.samples=5000,xtrain=sortthc2$x[trainthc],xtest=sortthc2$x[testthc],ytrain=sortthc2$y[trainthc],
              ytest=sortthc2$y[testthc],coastal="TRUE",simple="TRUE",m=NULL)
DIC2b<-gof2b$DIC
D2b<-gof2b$Dd
pred2b<-gof2b$pred
mse2b<-gof2b$MSE

gofssk<-goodnessoffit.simple(datatrain=sortthc2$log[trainthc], datatest=sortthc2$log[testthc],sigma=sigmassk,phi=rep(phithc,n.samples),tau=taussk,
                     beta0=beta0ssk,n.samples=5000,xtrain=sortthc2$x[trainthc],xtest=sortthc2$x[testthc],ytrain=sortthc2$y[trainthc],
                     ytest=sortthc2$y[testthc],coastal="FALSE",simple="TRUE",m=NULL)
DICssk<-gofssk$DIC
Dssk<-gofssk$Dd
predssk<-gofssk$pred
msessk<-gofssk$MSE

gof1b<-goodnessoffit.simple(datatrain=sortthc2$log[trainthc], datatest=sortthc2$log[testthc],sigma=sigma1b,phi=phi1b,tau=tau1b,
                            beta0=beta01b,n.samples=5000,xtrain=sortthc2$x[trainthc],xtest=sortthc2$x[testthc],ytrain=sortthc2$y[trainthc],
                            ytest=sortthc2$y[testthc],coastal="TRUE",simple="FALSE",m=m1b)
DIC1b<-gof1b$DIC
D1b<-gof1b$Dd
pred1b<-gof1b$pred
mse1b<-gof1b$MSE
goffsk<-goodnessoffit.simple(datatrain=sortthc2$log[trainthc], datatest=sortthc2$log[testthc],sigma=sigma1b,phi=phi1b,tau=tau1b,
                            beta0=beta01b,n.samples=5000,xtrain=sortthc2$x[trainthc],xtest=sortthc2$x[testthc],ytrain=sortthc2$y[trainthc],
                            ytest=sortthc2$y[testthc],coastal="FALSE",simple="FALSE",m=mfsk)
DICfsk<-goffsk$DIC
Dfsk<-goffsk$Dd
predfsk<-goffsk$pred
msefsk<-goffsk$MSE
