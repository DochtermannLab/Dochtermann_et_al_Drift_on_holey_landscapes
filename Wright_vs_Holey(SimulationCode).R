#### General Run Conditions ####
N=7500 #Number of individuals
n.traits=10 #Number of traits
h2=.8 #heritability
n.gen=100 #Number of generations
nsim=250 #Number of simulations

#### IVC Model (Wright LS) ####
prop=0.5 #proportion selected
eigen.stor.Wright<-matrix(NA,nsim,n.traits)
#Fitness landscape/ISS
for(j in 1:nsim){
Wright.LS<-rlkjcorr(1, n.traits , eta=1 )
ISS.var<-rep(1,n.traits) #this currently makes no difference because I'm using truncation selection
Wright.ISS<-cor2cov(ISS.var,Wright.LS)
optimum=matrix(rep(0,n.traits),ncol=1)

#starting values
mu.G=rep(0,n.traits)
mu.E=rep(0,n.traits)
Sigma.G=diag(n.traits)
Sigma.E=diag(n.traits)

eigen.stor<-matrix(NA,n.gen,n.traits)
#fit.stor <- matrix(NA,n.gen,2)

#main model
for (i in 1:n.gen){
  G<-mvrnorm(N,mu=mu.G,Sigma=Sigma.G)
  eigen.stor[i,]<-c(eigen(cov(G))$values)
  E<-mvrnorm(N,mu=mu.E,Sigma=Sigma.E)
  P<-h2*G+(1-h2)*E
  fit<-apply(P,1,fitness.AL.Lande1980,omega=Wright.ISS,optimum=optimum)
  #fit.stor[i,1] <- mean(fit)
  #fit.stor[i,2] <- sd(fit)
  pop<-cbind(fit,P,G,E)
  surv.pop<-trunc.select(pop,prop=prop)
  
  Sigma.G<-cov(surv.pop[,c((n.traits+2):(n.traits*2+1))])
  mu.G<-colMeans(surv.pop[,c((n.traits+2):(n.traits*2+1))])
  
}
eigen.stor.Wright[j,]<-eigen.stor[n.gen,]
}
Wright.12<-eigen.stor.Wright[,1]/eigen.stor.Wright[,2]

#### IVC Model (Holey LS) P = 0.5 ####
p=0.5 #holey landscape proportion parameter

eigen.stor.Holey<-matrix(NA,nsim,n.traits)

#Fitness landscape/ISS
HLS.phen<-permutations(2,n.traits,c(0,1),repeats.allowed = T)
phen.combs<-dim(HLS.phen)[1]
HLS.t<-c(t(HLS.phen))

for(j in 1:nsim){
HLS.fit<-sample(c(0,1),phen.combs,prob=c(1-p,p),replace=T)

#starting values
mu.G=rep(0,n.traits)
mu.E=rep(0,n.traits)
Sigma.G=diag(n.traits)
Sigma.E=diag(n.traits)

eigen.stor<-matrix(NA,n.gen,n.traits)

#main model
for (i in 1:n.gen){
  G<-mvrnorm(N,mu=mu.G,Sigma=Sigma.G)
  eigen.stor[i,]<-c(eigen(cov(G))$values)
  E<-mvrnorm(N,mu=mu.E,Sigma=Sigma.E)
  P<-h2*G+(1-h2)*E
  P2<-make.binary(P)#convert P to 0 or 1
  fit<-apply(P2,1,fitness.HL,mat=HLS.t,HLS.fit=HLS.fit)
  pop<-cbind(fit,P,G,E)
  surv.pop<-pop[which(pop[,1]==1),]
  
  Sigma.G<-cov(surv.pop[,c((n.traits+2):(n.traits*2+1))])
  mu.G<-colMeans(surv.pop[,c((n.traits+2):(n.traits*2+1))])
  
}
eigen.stor.Holey[j,]<-eigen.stor[n.gen,]
}
Holey.12<-eigen.stor.Holey[,1]/eigen.stor.Holey[,2]

#save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/Genetics.Exeter/IVCsimulation2.RData")

#### IVC Model (Holey LS) P = 0.2 ####
p=0.2
eigen.stor.Holey2<-matrix(NA,nsim,n.traits)

#Fitness landscape/ISS
for(j in 1:nsim){
  HLS.fit<-sample(c(0,1),phen.combs,prob=c(1-p,p),replace=T)
  
  #starting values
  mu.G=rep(0,n.traits)
  mu.E=rep(0,n.traits)
  Sigma.G=diag(n.traits)
  Sigma.E=diag(n.traits)
  
  eigen.stor<-matrix(NA,n.gen,n.traits)
  
  #main model
  for (i in 1:n.gen){
    G<-mvrnorm(N,mu=mu.G,Sigma=Sigma.G)
    eigen.stor[i,]<-c(eigen(cov(G))$values)
    E<-mvrnorm(N,mu=mu.E,Sigma=Sigma.E)
    P<-h2*G+(1-h2)*E
    P2<-make.binary(P)#convert P to threshold 0 or 1
    fit<-apply(P2,1,fitness.HL,mat=HLS.t,HLS.fit=HLS.fit)
    pop<-cbind(fit,P,G,E)
    surv.pop<-pop[which(pop[,1]==1),]
    
    Sigma.G<-cov(surv.pop[,c((n.traits+2):(n.traits*2+1))])
    mu.G<-colMeans(surv.pop[,c((n.traits+2):(n.traits*2+1))])
    
  }
  eigen.stor.Holey2[j,]<-eigen.stor[n.gen,]
}
Holey.12.2<-eigen.stor.Holey2[,1]/eigen.stor.Holey2[,2]

#save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/Genetics.Exeter/IVCsimulation2.RData")

#### IVC Model (Holey LS) P = 0.8 ####
p=0.8
eigen.stor.Holey8<-matrix(NA,nsim,n.traits)

#Fitness landscape/ISS
for(j in 1:nsim){
  HLS.fit<-sample(c(0,1),phen.combs,prob=c(1-p,p),replace=T)
  
  #starting values
  mu.G=rep(0,n.traits)
  mu.E=rep(0,n.traits)
  Sigma.G=diag(n.traits)
  Sigma.E=diag(n.traits)
  
  eigen.stor<-matrix(NA,n.gen,n.traits)
  
  #main model
  for (i in 1:n.gen){
    G<-mvrnorm(N,mu=mu.G,Sigma=Sigma.G)
    eigen.stor[i,]<-c(eigen(cov(G))$values)
    E<-mvrnorm(N,mu=mu.E,Sigma=Sigma.E)
    P<-h2*G+(1-h2)*E
    P2<-make.binary(P)#convert P to threshold 0 or 1
    fit<-apply(P2,1,fitness.HL,mat=HLS.t,HLS.fit=HLS.fit)
    pop<-cbind(fit,P,G,E)
    surv.pop<-pop[which(pop[,1]==1),]
    
    Sigma.G<-cov(surv.pop[,c((n.traits+2):(n.traits*2+1))])
    mu.G<-colMeans(surv.pop[,c((n.traits+2):(n.traits*2+1))])
    
  }
  eigen.stor.Holey8[j,]<-eigen.stor[n.gen,]
}
Holey.12.8<-eigen.stor.Holey8[,1]/eigen.stor.Holey8[,2]
#save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/Genetics.Exeter/IVCsimulation2.RData")


#### IVC Model (Drift only) ####
prop=0.5 #proportion surviving (same as Wright LS sims)
N=7500
n.traits=10
n.gen=100
nsim=250

eigen.stor.drift<-matrix(NA,nsim,n.traits)
#Fitness landscape/ISS
for(j in 1:nsim){
  
  #starting values
  mu.G=rep(0,n.traits)
  Sigma.G=diag(n.traits)
  
  eigen.stor<-matrix(NA,n.gen,n.traits)
  
  #main model
  for (i in 1:n.gen){
    G<-mvrnorm(N,mu=mu.G,Sigma=Sigma.G)
    eigen.stor[i,]<-c(eigen(cov(G))$values)
    pop<-G
    surv.pop <- pop[sample(N,N*prop,replace=F),] #randomly sample population for how survives
    Sigma.G<-cov(surv.pop[,c(1:n.traits)])
    mu.G<-colMeans(surv.pop[,c(1:n.traits)])
  }
  eigen.stor.drift[j,]<-eigen.stor[n.gen,]
}
hist(eigen.stor.drift[,2]/eigen.stor.drift[,1])

#### Model Summary ####

eigen.stor.Wright.std1<-eigen.stor.Wright/apply(eigen.stor.Wright,1,max)
eigen.stor.Holey.std1<-eigen.stor.Holey/apply(eigen.stor.Holey,1,max)
eigen.stor.Holey2.std1<-eigen.stor.Holey2/apply(eigen.stor.Holey2,1,max)
eigen.stor.Holey8.std1<-eigen.stor.Holey8/apply(eigen.stor.Holey8,1,max)
eigen.stor.Drift.std1<-eigen.stor.drift/apply(eigen.stor.drift,1,max)

comb.ls2<-as.data.frame(rbind(matrix(eigen.stor.Wright.std1[,2],ncol=1),
                              matrix(eigen.stor.Drift.std1[,2],ncol=1),
                              matrix(eigen.stor.Holey2.std1[,2],ncol=1),
                             matrix(eigen.stor.Holey.std1[,2],ncol=1),
                             matrix(eigen.stor.Holey8.std1[,2],ncol=1)))
names(comb.ls2)<-c("EigenRatio")
comb.ls2$LS<-c(rep("Wright",nsim),rep("Drift",nsim),rep("Holey p = 0.2",nsim),rep("Holey p = 0.5",nsim),
               rep("Holey p = 0.8",nsim))
comb.ls2$LS <- as.factor(comb.ls2$LS)
comb.ls2$LS <- factor(comb.ls2$LS,levels(comb.ls2$LS)[c(2,3,4,5,1)])


sum.ls2<-cbind(tapply(comb.ls2$EigenRatio, comb.ls2$LS, mean),
              tapply(comb.ls2$EigenRatio, comb.ls2$LS, quantile, probs=0.025),
              tapply(comb.ls2$EigenRatio, comb.ls2$LS, quantile, probs=0.975))
summary(aov(EigenRatio~LS,data=comb.ls2))
TukeyHSD(aov(EigenRatio~LS,data=comb.ls2))


#### 2 x 2 plot ####
#tiff('2by2.tiff',width=8,height=8,units='in',res=600)
par(mfrow=c(2,2),pty='s',mar=c(6,7.5,2.5,.5))

sim.mean<-colMeans(eigen.stor.Holey2/apply(eigen.stor.Holey2,1,max))
sim.sd<-apply(eigen.stor.Holey2/apply(eigen.stor.Holey2,1,max),2,sd)
plot(NA,xlim=c(1,10),ylim=c(0,1), xlab=" ", 
     ylab=" ",cex.axis=1.25)
arrows(x0 = c(1:10), y0=sim.mean-sim.sd, x1 = c(1:10),
       y1=sim.mean+sim.sd,length = .075,angle=90,code=3)
mtext("Holey (p = 0.2)",3,line=.5,cex=1.35)
mtext(expression(paste(frac(lambda[i],"max("*lambda*")"))),2,
      line=2.5,las=2,cex=1.5)
points(y=colMeans(eigen.stor.Holey2/apply(eigen.stor.Holey2,1,max)),
       x=c(1:10),cex=2.5,lwd=1.5,pch=21,bg="white")

sim.mean<-colMeans(eigen.stor.Holey/apply(eigen.stor.Holey,1,max))
sim.sd<-apply(eigen.stor.Holey/apply(eigen.stor.Holey,1,max),2,sd)
plot(NA,xlim=c(1,10),ylim=c(0,1), xlab=" ", 
     ylab=" ",cex.axis=1.25)
arrows(x0 = c(1:10), y0=sim.mean-sim.sd, x1 = c(1:10),
       y1=sim.mean+sim.sd,length = .075,angle=90,code=3)
mtext("Holey (p = 0.5)",3,line=.5,cex=1.35)
points(y=colMeans(eigen.stor.Holey/apply(eigen.stor.Holey,1,max)),
       x=c(1:10),cex=2.5,lwd=1.5,pch=21,bg="white")


sim.mean<-colMeans(eigen.stor.Holey8/apply(eigen.stor.Holey8,1,max))
sim.sd<-apply(eigen.stor.Holey8/apply(eigen.stor.Holey8,1,max),2,sd)
plot(NA,xlim=c(1,10),ylim=c(0,1), xlab=" ", 
     ylab=" ",cex.axis=1.25)
arrows(x0 = c(1:10), y0=sim.mean-sim.sd, x1 = c(1:10),
       y1=sim.mean+sim.sd,length = .075,angle=90,code=3)
mtext(expression(atop("Dimension", paste("(ordered by ",lambda[i],")"))),1,
      line=5,cex=1.5)
mtext("Holey (p = 0.8)",3,line=.5,cex=1.35)
mtext(expression(paste(frac(lambda[i],"max("*lambda*")"))),2,
      line=2.5,las=2,cex=1.5)
points(y=colMeans(eigen.stor.Holey8/apply(eigen.stor.Holey8,1,max)),
       x=c(1:10),cex=2.5,lwd=1.5,pch=21,bg="white")

sim.mean<-colMeans(eigen.stor.Wright/apply(eigen.stor.Wright,1,max))
sim.sd<-apply(eigen.stor.Wright/apply(eigen.stor.Wright,1,max),2,sd)
plot(NA,xlim=c(1,10),ylim=c(0,1), xlab=" ", 
     ylab=" ",cex.axis=1.25)
arrows(x0 = c(1:10), y0=sim.mean-sim.sd, x1 = c(1:10),
       y1=sim.mean+sim.sd,length = .075,angle=90,code=3)
mtext("Wrightian",3,line=.5,cex=1.35)
mtext(expression(atop("Dimension", paste("(ordered by ",lambda[i],")"))),1,
      line=5,cex=1.5)
points(y=colMeans(eigen.stor.Wright/apply(eigen.stor.Wright,1,max)),
       x=c(1:10),cex=2.5,lwd=1.5,pch=21,bg="white")
#dev.off()
#### End ####
#save.image("C:/Users/Ned Dochtermann/Dropbox/Working/Projects/HoleyLandscapes/IVCsimulation4.RData")

