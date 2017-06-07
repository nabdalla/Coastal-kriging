num.sim=100
set.seed(111)
#curve
t<-runif(num.sim, 0, pi*2)
t<-sort(t)
x<-cos(t)
y<-sin(t)
xy<-cbind(x, y)
plot(x, y)
#line segment approximation
t2<-NULL
for ( j in 2:length(t)){
  t2[j-1]<-sqrt((abs(x[j]-x[j-1]))^2+(abs(y[j]-y[j-1]))^2)
}
t2<-cumsum(c(t[1],t2))
#t2<-(t2-mean(t2))/sqrt(var(t2))
#t<-(t-mean(t))/sqrt(var(t))
plot(t,t2)
abline(c(0,1))
#variance of the response variable
distance<-abs(outer(t ,t,"-"))
tau.p<-0.1
sigma.p<-1
vv2<-sigma.p*exp(-1*distance)
beta.p<-0
#simulation fo the response variable
z<-NULL
x0<-rep(1,num.sim)
set.seed(111)
z<-rmvnorm(1, beta.p%*%t(cbind(x0)),vv2+tau.p*diag(num.sim) )
#Variogram to get initial values for the full model and fixed values for the simplified model
curve<-data.frame(t(z))
vario <- variog(coords = cbind(t, rep(0,num.sim)), data = curve)
initial.values <-expand.grid(seq(0,1,l=5), seq(0,1,l=5))
fit <- variofit(vario, ini.cov.pars = initial.values,  cov.model = "exponential")
sigma<-fit$cov.pars[1]
phi<-fit$cov.pars[2]
tau<-fit$nugget
alpha <- tau/fit$cov.pars[1]
dev.off()
plot(vario)
lines(fit)
#Variogram to get initial values for the simple kriging
curve<-data.frame(t(z))
varios <- variog(coords = cbind(x, y), data = curve)
fits <- variofit(varios, ini.cov.pars = initial.values,  cov.model = "exponential")
sigmas<-fits$cov.pars[1]
phis<-fits$cov.pars[2]
taus<-fits$nugget+0.01
alphas <- taus/fits$cov.pars[1]
dev.off()
plot(varios)
lines(fits)
#prior distributions parameters for beta and sigma
p<-1
beta.prior.mean <- as.matrix(rep(0, times=p))
beta.prior.precision <- matrix(0, nrow=p, ncol=p)
sigma.sq.prior.shape <- 2
sigma.sq.prior.rate <-2
#choose the training and hold out samples
test<-seq(1,num.sim,4)
train<-(seq(1,num.sim,1)[-test])
n1<-length(train)
n2<-length(test)
#arc length model 2a
exact1 <- bayesGeostatExact(z[train]~1, n.samples=n.samples,
                            beta.prior.mean=beta.prior.mean,
                            beta.prior.precision=beta.prior.precision,
                            coords=cbind(t[train], rep(0,n1)), phi=phi, alpha=alpha,
                            sigma.sq.prior.shape=sigma.sq.prior.shape,
                            sigma.sq.prior.rate=sigma.sq.prior.rate)
#line segment approximation model 2b
exact2 <- bayesGeostatExact(z[train]~1, n.samples=n.samples,
                            beta.prior.mean=beta.prior.mean,
                            beta.prior.precision=beta.prior.precision,
                            coords=cbind(t2[train], rep(0,n1)), phi=phi, alpha=alpha,
                            sigma.sq.prior.shape=sigma.sq.prior.shape,
                            sigma.sq.prior.rate=sigma.sq.prior.rate)
exactsimple <- bayesGeostatExact(z[train]~1, n.samples=n.samples,
                                 beta.prior.mean=beta.prior.mean,
                                 beta.prior.precision=beta.prior.precision,
                                 coords=cbind(x[train], y[train]), phi=phis, alpha=alphas,
                                 sigma.sq.prior.shape=sigma.sq.prior.shape,
                                 sigma.sq.prior.rate=sigma.sq.prior.rate)
print(summary(exact1$p.samples))
print(summary(exact2$p.samples))
print(summary(exactsimple$p.samples))

beta02a<-(exact1$p.samples[,1])
sigma2a<-(exact1$p.samples[,2])
tau2a<-(exact1$p.samples[,3])
beta02b<-(exact2$p.samples[,1])
sigma2b<-(exact2$p.samples[,2])
tau2b<-(exact2$p.samples[,3])
beta0s<-(exactsimple$p.samples[,1])
sigmas<-(exactsimple$p.samples[,2])
taus<-(exactsimple$p.samples[,3])
#Predictoin using test sample model 2a
tpt<-c(t[train],t[test])
dista<-abs(outer(tpt,tpt,"-"))
v2<-vv.p2<-muhat2<-yhat<-uarc<-larc<-vector("list")
for(i in 1:n.samples){
  v2[[i]]<-sigma2a[i]*exp(-phi*dista)+tau2a[i]*diag(n1+n2)
  vv.p2[[i]]<- v2[[i]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v2[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2[[i]][1:n1,1:n1]))%*%(t(v2[[i]][(n1+1):(n2+n1), 1:n1]))
  muhat2[[i]]<-t(beta02a[i]%*%t(rep(1, (n2))))+(v2[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2[[i]][1:n1,1:n1]))%*%t((t(z[train])-beta02a[i]%*%t(x0[1:n1])))
  set.seed(i)
  yhat[[i]]<-rmvnorm(1,muhat2[[i]],round(vv.p2[[i]],5))
  uarc[[i]]<-yhat[[i]]+1.96*sqrt(diag(vv.p2[[i]]))
  larc[[i]]<-yhat[[i]]-1.96*sqrt(diag(vv.p2[[i]]))
}
yhat.mat<-matrix(unlist(yhat), nrow=n.samples, ncol=n2, byrow=T)
uarc.mat<-matrix(unlist(uarc), nrow=n.samples, ncol=n2, byrow=T)
larc.mat<-matrix(unlist(larc), nrow=n.samples, ncol=n2, byrow=T)

predarc<-upperarc<-lowerarc<-NULL
for (i in 1:n2){
  predarc[i]<-mean(yhat.mat[,i])
  upperarc[i]<-mean(uarc.mat[,i])
  lowerarc[i]<-mean(larc.mat[,i])
}

#prediction test data model 2b
tp<-c(t2[train],t2[test])
distb<-abs(outer(tp,tp,"-"))
v<-vv.p<-muhat<-yhat2<-uline<-lline<-vector("list")
for(i in 1:n.samples){
  v[[i]]<-sigma2b[i]*exp(-phi*distb)+tau2b[i]*diag(n1+n2)
  vv.p[[i]]<- v[[i]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v[[i]][1:n1,1:n1]))%*%(t(v[[i]][(n1+1):(n2+n1), 1:n1]))
  muhat[[i]]<-t(beta02b[i]%*%t(rep(1, (n2))))+(v[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v[[i]][1:n1,1:n1]))%*%t((t(z[train])-beta02b[i]%*%t(x0[1:n1])))
  set.seed(i)
  yhat2[[i]]<-rmvnorm(1,muhat[[i]],vv.p[[i]])
  uline[[i]]<-yhat2[[i]]+1.96*sqrt(diag(vv.p[[i]]))
  lline[[i]]<-yhat2[[i]]-1.96*sqrt(diag(vv.p[[i]]))
}
yhatline.mat<-matrix(unlist(yhat2), nrow=n.samples, ncol=n2, byrow=T)
uline.mat<-matrix(unlist(uline), nrow=n.samples, ncol=n2, byrow=T)
lline.mat<-matrix(unlist(lline), nrow=n.samples, ncol=n2, byrow=T)

predline<-upperline<-lowerline<-NULL
for (i in 1:n2){
  predline[i]<-mean(yhatline.mat[,i])
  upperline[i]<-mean(uline.mat[,i])
  lowerline[i]<-mean(lline.mat[,i])
  
}

#prediction test data model simple
dist<-as.matrix(dist(data.frame(cbind(x,y)),upper=TRUE, diag=TRUE,method="euclidean"))
v<-vv.p<-muhat<-yhat2<-uline<-lline<-vector("list")
for(i in 1:n.samples){
  v[[i]]<-sigmas[i]*exp(-phis*dist)+taus[i]*diag(n1+n2)
  vv.p[[i]]<- v[[i]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v[[i]][1:n1,1:n1]))%*%(t(v[[i]][(n1+1):(n2+n1), 1:n1]))
  muhat[[i]]<-t(beta0s[i]%*%t(rep(1, (n2))))+(v[[i]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v[[i]][1:n1,1:n1]))%*%t((t(z[train])-beta0s[i]%*%t(x0[1:n1])))
  set.seed(i)
  yhat2[[i]]<-rmvnorm(1,muhat[[i]],vv.p[[i]])
  uline[[i]]<-yhat2[[i]]+1.96*sqrt(diag(vv.p[[i]]))
  lline[[i]]<-yhat2[[i]]-1.96*sqrt(diag(vv.p[[i]]))
}
yhatline.mats<-matrix(unlist(yhat2), nrow=n.samples, ncol=n2, byrow=T)
uline.mats<-matrix(unlist(uline), nrow=n.samples, ncol=n2, byrow=T)
lline.mats<-matrix(unlist(lline), nrow=n.samples, ncol=n2, byrow=T)

predlines<-upperlines<-lowerlines<-NULL
for (i in 1:n2){
  predlines[i]<-mean(yhatline.mats[,i])
  upperlines[i]<-mean(uline.mats[,i])
  lowerlines[i]<-mean(lline.mats[,i])
  
}

plotCI(z[test], predline,ui=upperline, li=lowerline)
abline(0,1)
plotCI(z[test], predarc,ui=upperarc, li=lowerarc)
abline(0,1)
plotCI(z[test], predlines,ui=upperlines, li=lowerlines)
abline(0,1)
#Full hierarchical models 1a and 1b 
starting <- list("phi"=phi, "sigma.sq"=sigma, "tau.sq"=tau)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                 "phi.Unif"=c(1, 3/0.1), "sigma.sq.IG"=c(2,2),
                 "tau.sq.IG"=c(2, 2))
#model 1a
sp.exact1 <-spLM(z[train]~1,
                 coords=cbind(t[train], rep(0,n1)), tuning=tuning, starting=starting, n.samples=n.samples,
                 priors=priors.1,cov.model="exponential")
#model 1b
sp.exact2 <-spLM(z[train]~1,
                 coords=cbind(t2[train], rep(0,n1)), tuning=tuning, starting=starting, n.samples=n.samples,
                 priors=priors.1,cov.model="exponential")
sp.exactsimple <-spLM(z[train]~1,
                      coords=cbind(x[train], y[train]), tuning=tuning, starting=starting, n.samples=n.samples,
                      priors=priors.1,cov.model="exponential")

burn.in <- 0.5*n.samples
m.11 <- spRecover(sp.exact1 , start=burn.in, verbose=FALSE)
m.12<- spRecover(sp.exact2 , start=burn.in, verbose=FALSE)
m.1simple<- spRecover(sp.exactsimple , start=burn.in, verbose=FALSE)

round(summary(m.11$p.theta.recover.samples)$quantiles,3)
round(summary(m.12$p.theta.recover.samples)$quantiles,3)
round(summary(m.11$p.beta.recover.samples)$quantiles,3)
round(summary(m.12$p.beta.recover.samples)$quantiles,3)

##Produce the posterior summaries
beta0sp1<-mean(m.11$p.beta.recover.samples[,1])
beta0sp2<-mean(m.12$p.beta.recover.samples[,1])
#Prediction using test data
pred.covars <- cbind(rep(1,(n2) ))
m.pred1 <- spPredict(m.11, pred.covars, 
                     pred.coords=cbind(t[test], rep(0,(n2))),
                     start=0.5*n.samples)
m.pred2 <- spPredict(m.12, pred.covars, 
                     pred.coords=cbind(t2[test], rep(0,(n2))),
                     start=0.5*n.samples)
m.predsimple <- spPredict(m.1simple, pred.covars, 
                          pred.coords=cbind(x[test], y[test]),
                          start=0.5*n.samples)
y.hatsimple <- as.numeric(apply(m.predsimple$p.y.predictive.samples, 1, mean))
y.hat1 <- as.numeric(apply(m.pred1$p.y.predictive.samples, 1, mean))
sigma.full1<-m.11$p.theta.recover.samples[,1]
tau.full1<-m.11$p.theta.recover.samples[,2]
phi.full1<-m.11$p.theta.recover.samples[,3]
beta.full1<-m.11$p.beta.recover.samples[,1]
y.hat2 <- as.numeric(apply(m.pred2$p.y.predictive.samples, 1, mean))
sigma.full2<-m.12$p.theta.recover.samples[,1]
tau.full2<-m.12$p.theta.recover.samples[,2]
phi.full2<-m.12$p.theta.recover.samples[,3]
beta.full2<-m.12$p.beta.recover.samples[,1]
#Credible intervals 
uline1<-uline2<-lline1<-lline2<-uarcfull<-larcfull<-ulinefull<-llinefull<-variance<-variance2<-NULL
v1<-vv.p1<-v2<-vv.p2<-vector("list")
for(j in 1:(n.samples*0.5)){
  v1[[j]]<-sigma.full1[j]*exp(-phi.full1[j]*dista)+tau.full1[j]*diag(n1+n2)
  vv.p1[[j]]<- v1[[j]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v1[[j]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v1[[j]][1:n1,1:n1]))%*%(t(v1[[j]][(n1+1):(n2+n1), 1:n1]))
  v2[[j]]<-sigma.full2[j]*exp(-phi.full2[j]*distb)+tau.full2[j]*diag(n1+n2)
  vv.p2[[j]]<- v2[[j]][(n1+1):(n2+n1), (n1+1):(n2+n1)]-(v2[[j]][(n1+1):(n2+n1), 1:n1])%*%
    chol2inv(chol(v2[[j]][1:n1,1:n1]))%*%(t(v2[[j]][(n1+1):(n2+n1), 1:n1]))
  uline1[[i]]<-1.96*sqrt(diag(vv.p1[[j]]))
  lline1[[i]]<- -1.96*sqrt(diag(vv.p1[[j]]))
  uline2[[i]]<-1.96*sqrt(diag(vv.p2[[j]]))
  lline2[[i]]<- -1.96*sqrt(diag(vv.p2[[j]]))
}

uline.mat1<-matrix(unlist(uline1), nrow=n.samples*0.5, ncol=n2, byrow=T)
lline.mat1<-matrix(unlist(lline1), nrow=n.samples*0.5, ncol=n2, byrow=T)
uline.mat2<-matrix(unlist(uline2), nrow=n.samples*0.5, ncol=n2, byrow=T)
lline.mat2<-matrix(unlist(lline2), nrow=n.samples*0.5, ncol=n2, byrow=T)
upperline1<-lowerline1<-upperline2<-lowerline2<-NULL
for (i in 1:n2){
  upperline1[i]<- mean(uline.mat1[,i])+ y.hat1[i]
  lowerline1[i]<- y.hat1[i]+mean(lline.mat1[,i])
  upperline2[i]<-mean(uline.mat2[,i])+ y.hat2[i]
  lowerline2[i]<- y.hat2[i]+mean(lline.mat2[,i])
}
#G+P=D
v1a<-v1b<-v2a<-v2b<-vs<-vector("list")
yhat1a<-yhat1b<-yhat2a<-yhat2b<-yhats<-matrix(0,nrow=(0.5*n.samples), ncol=num.sim)
for(j in 1:(0.5*n.samples)){
    v1a[[j]]<-sigma.full1[j]*exp(-phi.full1[j]*dista)+tau.full1[j]*diag(num.sim)
    yhat1a[j,1:num.sim]<-rmvnorm(1,rep(beta.full1[j],num.sim),v1a[[j]])
    v1b[[j]]<-sigma.full2[j]*exp(-phi.full2[j]*distb)+tau.full2[j]*diag(num.sim)
    yhat1b[j,1:num.sim]<-rmvnorm(1,rep(beta.full2[j],num.sim),v1b[[j]])
    v2a[[j]]<-sigma2a[j+2400]*exp(-phi*dista)+tau2a[j+2400]*diag(num.sim)
    yhat2a[j,1:num.sim]<-rmvnorm(1,rep(beta02a[j+2400],num.sim),v2a[[j]])
    v2b[[j]]<-sigma2b[j+2400]*exp(-phi*distb)+tau2b[j+2400]*diag(num.sim)
    yhat2b[j,1:num.sim]<-rmvnorm(1,rep(beta02b[j+2400],num.sim),v2b[[j]])
    vs[[j]]<-sigmas[j+2400]*exp(-phi*dist)+taus[j+2400]*diag(num.sim)
    yhats[j,1:num.sim]<-rmvnorm(1,rep(beta0s[j+2400],num.sim),vs[[j]])
}
plot(density(yhat2b))
murep1a<-apply(yhat1a,2 , mean)
murep1b<-apply(yhat1b,2 , mean)
murep2a<-apply(yhat2a,2 , mean)
murep2b<-apply(yhat2b,2 , mean)
mureps<-apply(yhats,2 , mean)
P1a<-apply(yhat1a,2 , var)
P1b<-apply(yhat1b,2 , var)
P2a<-apply(yhat2a,2 , var)
P2b<-apply(yhat2b,2 , var)
Ps<-apply(yhats,2 , var)
G1a<-sum((z[c(train,test)]-murep1a)^2)
G1b<-sum((z[c(train,test)]-murep1b)^2)
G2a<-sum((z[c(train,test)]-murep2a)^2)
G2b<-sum((z[c(train,test)]-murep2b)^2)
Gs<-sum((z[c(train,test)]-mureps)^2)
D1a<-sum(P1a)+G1a
D1b<-sum(P1b)+G1b
D2a<-sum(P2a)+G2a
D2b<-sum(P2b)+G2b
Ds<-sum(Ps)+Gs
D1a
D1b
D2a
D2b
Ds
#Kullback-Liberur
kl1a<-kl1b<-kl2a<-kl2b<-kls<-NULL
v0<-vv2[c(train,test),c(train,test)]+tau.p*diag(num.sim)
for(j in 1:(0.5*n.samples)){
  v1a[[j]]<-sigma.full1[j]*exp(-phi.full1[j]*dista)+tau.full1[j]*diag(num.sim)
  kl1a[j]<-tr(solve(v1a[[j]])%*%(v0))+t(rep(beta.full1[j],num.sim))%*%
    solve(v1a[[j]])%*%(rep(beta.full1[j],num.sim))-num.sim+log(det(v1a[[j]]))-log(det(v0))
  v1b[[j]]<-sigma.full2[j]*exp(-phi.full2[j]*distb)+tau.full2[j]*diag(num.sim)
  kl1b[j]<-tr(solve(v1b[[j]])%*%(v0))+t(rep(beta.full2[j],num.sim))%*%
    solve(v1b[[j]])%*%(rep(beta.full2[j],num.sim))-num.sim+log(det(v1b[[j]]))-log(det(v0))
  v2a[[j]]<-sigma2a[j+2400]*exp(-phi*dista)+tau2a[j+2400]*diag(num.sim)
  kl2a[j]<-tr(solve(v2a[[j]])%*%(v0))+t(rep(beta02a[j+2400],num.sim))%*%
    solve(v2a[[j]])%*%rep(beta02a[j+2400],num.sim)-num.sim+log(det(v2a[[j]]))-log(det(v0))
  v2b[[j]]<-sigma2b[j+2400]*exp(-phi*distb)+tau2b[j+2400]*diag(num.sim)
  kl2b[j]<-tr(solve(v2b[[j]])%*%(v0))+t(rep(beta02b[j+2400],num.sim))%*%
    solve(v2b[[j]])%*%rep(beta02b[j+2400],num.sim)-num.sim+log(det(v2b[[j]]))-log(det(v0))
  vs[[j]]<-sigmas[j+2400]*exp(-phi*dist)+taus[j+2400]*diag(num.sim)
  kls[j]<-tr(solve(vs[[j]])%*%(v0))+t(rep(beta0s[j+2400],num.sim))%*%
    solve(vs[[j]])%*%rep(beta0s[j+2400],num.sim)-num.sim+log(det(vs[[j]]))-log(det(v0))}
print(meankl1a<-0.5*mean(kl1a))
print(meankl1b<-0.5*mean(kl1b))
print(meankl2a<-0.5*mean(kl2a))
print(meankl2b<-0.5*mean(kl2b))
print(meankls<-0.5*mean(kls))
plot(density(kl2b))
#Plots with 95%CI
quartz()
layout(matrix(c(1,1,1,1,1,0,2,2,2,2,2,3,3,3,3,3,0,4,4,4,4,4,0,0,0,5,5,5,5,5,0,0,0), 3, 11, byrow = TRUE))
plotCI(z[test], y.hat1,ui=upperline1, li=lowerline1,xlab="True", ylab="Predicted",
       main="Model 1a",
       cex.main=0.9,cex=0.9)
abline(0,1)
plotCI(z[test], y.hat2,ui=upperline2, li=lowerline2,xlab="True", ylab="Predicted",
       main="Model 1b",
       cex.main=0.9,cex=0.9)
abline(0,1)
plotCI(z[test], predarc,ui=upperarc, li=lowerarc,xlab="True", ylab="Predicted",
       main="Model 2a",
       cex.main=0.9,cex=0.9)
abline(0,1)
plotCI(z[test], predline,ui=upperline, li=lowerline,xlab="True", ylab="Predicted",
       main="Model 2b",
       cex.main=0.9,cex=0.9)
abline(0,1)
plotCI(z[test], predlines,ui=upperlines, li=lowerlines,xlab="True", ylab="Predicted",
       main="Simple Kriging",
       cex.main=0.9,cex=0.9)
abline(0,1)
plot(z[test],y.hatsimple)
abline(0,1)
#mean square error
mse2a<-sum((z[test]-predline)^2)/n2
mse2b<-sum((z[test]-predarc)^2)/n2
mse2s<-sum((z[test]-predlines)^2)/n2
mse1a<-sum((z[test]-y.hat1)^2)/n2
mse1b<-sum((z[test]-y.hat2)^2)/n2
