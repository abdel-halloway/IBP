rm(list=ls(all=TRUE))
graphics.off()

this.dir<-dirname(parent.frame(2)$ofile)
work.dir<-paste0(this.dir,'/Figures/')
dir.create(work.dir,recursive=TRUE)
setwd(work.dir)

IBP.criterion<-function(params,b.,k.,h.,u.,vigilance){
  if(vigilance=="optimized"){
    r<-params['r']; K<-params['K']
    d<-params['d']; c<-params['c']
    m<-params['m']; k<-params['k']; b<-params['b']
    
    if(k.>b.){
      IBP.crit<-r*(c/d-k./(m*K))
    }
    else{
      IBP.crit<-r/(2*b.)*(b.+k.)*(c/d-(k.+b.)/(2*m*K))
    }
  }else{
    r<-params['r']; K<-params['K']
    d<-params['d']; c<-params['c']
    m<-params['m']; k<-params['k']; b<-params['b']
    
    IBP.crit<-r*c*(1-u.)*(1/d-(k.+b.*u.)/(m*(c-d*h.)*K))
  }
  if(IBP.crit<0){
    IBP.crit<-0
  }
  return(IBP.crit)
}

IBP.criterion.vec<-Vectorize(IBP.criterion,c("b.","k."))

parameters<-c(r=0.80,K=100,d=0.12,c=1,m=0.04,k=15,b=30)

pdf("Figure 3.pdf",width=6,height=6)
par(mfrow=c(2,2))
b.range<-k.range<-seq(0,200,1)
plot.new()
plot.window(xlim=c(min(b.range),max(b.range)),ylim=c(0,4))
axis(1,pos=0)
axis(2,pos=0)
lines(b.range,IBP.criterion.vec(parameters,b.range,parameters['k'],0,0,'optimized'),type='l',lwd=1.5)
abline(h=1,lty=3)
title(ylab=expression(IBP~criterion~~'('*italic(over(P^'*',N^'*'))*')'),line=1.5)
title(xlab=expression(Effectiveness~of~Vigilance~'('*italic(b)*')'))

plot.new()
plot.window(xlim=c(min(k.range),max(k.range)),ylim=c(0,4))
axis(1,pos=0)
axis(2,pos=0)
lines(k.range,IBP.criterion.vec(parameters,parameters['b'],k.range,0,0,'optimized'),type='l',lwd=1.5)
abline(h=1,lty=3)
title(ylab=expression(IBP~criterion~~'('*italic(over(P^'*',N^'*'))*')'),line=1.5)
title(xlab=expression(Inverse~Predator~Lethality~'('*italic(k)*')'))

plot.new()
plot.window(xlim=c(min(b.range),max(b.range)),ylim=c(0,0.5))
axis(1,pos=0)
axis(2,pos=0)
u.star.b<-(b.range-parameters['k'])/(2*b.range)
lines(b.range,u.star.b,type='l',lwd=1.5)
abline(h=1,lty=3)
title(ylab=expression(Optimal~Vigilance~"("*italic(u^"*")*")"),line=1.5)
title(xlab=expression(Effectiveness~of~Vigilance~'('*italic(b)*')'))

plot.new()
plot.window(xlim=c(min(k.range),max(k.range)),ylim=c(0,0.5))
axis(1,pos=0)
axis(2,pos=0)
u.star.k<-(parameters['b']-k.range)/(2*parameters['b'])
u.star.k[u.star.k<0]<-0
lines(k.range,u.star.k,type='l',lwd=1.5)
abline(h=1,lty=3)
title(ylab=expression(Optimal~Vigilance~"("*italic(u^"*")*")"),line=1.5)
title(xlab=expression(Inverse~Predator~Lethality~'('*italic(k)*')'))
dev.off()
