rm(list=ls(all=TRUE))
graphics.off()

this.dir<-dirname(parent.frame(2)$ofile)
work.dir<-paste0(this.dir,'/New Figures')
dir.create(work.dir,recursive=TRUE)
setwd(work.dir)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Fig 1. Type I and Type II w/o vig. (h=0, h=1, h=2... etc)
# --General Pred, Prey Isocline and P* as fn of N,h,u
# 
# Fig 2. Conceptual diagram of qualitative effects of vigilance 
# --For Meg
# 
# Fig 3. Fixed Vigilance 
# --General Pred, Prey Isocline and P* as a fn of N,h,u
# 
# Fig 4. Type I w/ and w/o optimized vig
# --How does opt vig change the eq point...
# 
# Fig 5. 3D vig w/ N v P
# --u* as a fn of N,P

prey.isocline<-function(params,N.,h.,u.,vigilance){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  u<-u.; h<-h.
  
  a<-m/(k+b*u)
  
  if(vigilance=="fixed"){
    x<-r/a*(1-u)*(1+a*h*N.)*(1-N./K)
  }else if(vigilance=="optimized"){
    x<-r/(4*m*b)*(b+k+m*h*N.)^2*(1-N./K)
  }

  return(x)
}

pred.isocline<-function(params,N.,h.,u.,vigilance){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  u<-u.; h<-h.
  
  a<-m/(k+b*u)
  
  if(vigilance=="fixed"){
    x<-d/(a*(c-d*h))
  }else if(vigilance=="optimized"){
    x<-r*m/b*(c/d*N.)^2*(1-N./K)
  }

  return(x)
}

u_star<-function(N.,P.,params){
  r<-params['r']; K<-params['K']
  m<-params['m']; k<-params['k']; b<-params['b']
  
  u_star<-sqrt(m*P./(r*b*(1-N./K)))-k/b
  if(u_star>1 | u_star==Inf){
    u_star<-1
  }else if(u_star<0){
    u_star<-0
  }
  
  return(u_star)
}

N_star<-function(params,h.,u.,vigilance){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  u<-u.; h<-h.
  
  a<-m/(k+b*u)
  
  if(vigilance=="fixed"){
    N_star<-d/(a*(c-d*h))
  }else if(vigilance=="optimized"){
    N_star<-d*(b+k)/(m)/(2*c-d*h)
  }
  
  return(N_star)
}

P_star<-function(params,h.,u.,vigilance){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  u<-u.; h<-h.
  
  a<-m/(k+b*u)
  
  if(vigilance=="fixed"){
    P_star<-r/a*(1-u)*(1+d*h/(c-d*h))*(1-d/(a*K*(c-d*h)))
  }else if(vigilance=="optimized"){
    P_star<-r/(4*m*b)*(1-d*(b+k)/(m*K)/(2*c-d*h))*(b+k+d*h*(b+k)/(2*c-d*h))^2
  }
  
  return(P_star)
}

u_star<-function(N.,P.,params){
  r<-params['r']; K<-params['K']
  m<-params['m']; k<-params['k']; b<-params['b']
  
  u_star<-sqrt(m*P./(r*b*(1-N./K)))-k/b
  if(u_star<0 | is.nan(u_star) | is.na(u_star)){
    u_star<-0
  }else if(u_star>1 | is.infinite(u_star)){
    u_star<-1
  }
  
  return(u_star)
}

plot.IBP<-function(params,N.,h.,u.,type,vigilance,...){
  if(vigilance=="optimized"){
    pred.no.vig.params<-replace(params,"u",0); pred.total.vig.params<-replace(params,"u",1)
    pred.N<-seq(pred.isocline(pred.no.vig.params,N.,h.,pred.no.vig.params["u"],"fixed"),pred.isocline(pred.total.vig.params,N.,h.,pred.total.vig.params["u"],"fixed"),0.01)
    prey.iso<-prey.isocline(params,N.,h.,u.,vigilance)
    pred.iso<-pred.isocline(params,pred.N,h.,u.,vigilance)
    lines(N.,prey.iso,lty=2,col="blue",...)
    if(type=="one"){
      N_one_star<-expression(frac(dk,cm))
      N_two_star<-expression(frac(d(k+b),cm))
    }else if(type=="two"){
      N_one_star<-expression(frac(dk,m(c-dh)))
      N_two_star<-expression(frac(d(k+bu),m(c-dh)))
    }
    segments(pred.N[1],0,pred.N[1],pred.iso[1],col="red",lty=2)
    segments(tail(pred.N,1),tail(pred.iso,1),tail(pred.N,1),1e9,col="red",lty=2)
    lines(pred.N,pred.iso,col="red",lty=2)
    segments(0,0,500,500,lty=3)
  }else{
    prey.iso<-prey.isocline(params,N.,h.,u.,vigilance)
    pred.iso<-pred.isocline(params,N.,h.,u.,vigilance)
    
    plot(N.,prey.iso,type='l',xlim=c(0,params['K']),xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",axes=FALSE,col="darkblue",...)
    abline(v=pred.iso,lty=1,col="darkred")
    abline(a=0,b=1,lty=3)
  }
}

pred.lines.IBP<-function(params,N.,h.,u.,type,vigilance,...){
  if(vigilance=="optimized"){
    pred.no.vig.params<-replace(params,"u",0); pred.total.vig.params<-replace(params,"u",1)
    pred.N<-seq(pred.isocline(pred.no.vig.params,N.,h.,pred.no.vig.params["u"],"fixed"),pred.isocline(pred.total.vig.params,N.,h.,pred.total.vig.params["u"],"fixed"),0.01)
    prey.iso<-prey.isocline(params,N.,h.,u.,vigilance)
    pred.iso<-pred.isocline(params,pred.N,h.,u.,vigilance)
    segments(pred.N[1],0,pred.N[1],pred.iso[1],col="red",lty=2)
    segments(tail(pred.N,1),tail(pred.iso,1),tail(pred.N,1),1e9,col="red",lty=2)
    lines(pred.N,pred.iso,col="red",lty=2)
  }else{
    prey.iso<-prey.isocline(params,N.,h.,u.,vigilance)
    pred.iso<-pred.isocline(params,N.,h.,u.,vigilance)
    
    plot(N.,prey.iso,type='l',xlim=c(0,params['K']),xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",axes=FALSE,col="darkblue",...)
    abline(v=pred.iso,lty=1,col="darkred")
    abline(a=0,b=1,lty=3)
  }
}

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

u.isoleg<-function(params,N.,h.,u.){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  
  P<-r/(b*m)*(b*u.+k+m*h.*N.)^2*(1-N./K)
  return(P)
}

prey.isocline.vec<-Vectorize(prey.isocline,c("N.","h.","u."))
pred.isocline.vec<-Vectorize(pred.isocline,c("N.","h.","u."))
P_star.vec<-Vectorize(P_star,c("h.","u."))
u_star.vec<-Vectorize(u_star,c("N.","P."))
u.isoleg.vec<-Vectorize(u.isoleg,c("N.","h.","u."))
IBP.criterion.vec<-Vectorize(IBP.criterion,c("b.","k."))

parameters<-c(r=0.13,K=100,d=0.04,c=1,m=0.04,k=15,b=60)
color.gradient<-colorRampPalette(c("gold3","dimgray","blue"))
prey.gradient<-colorRampPalette(c('blue1','blue4'))
pred.gradient<-colorRampPalette(c('firebrick1','firebrick4'))
gold.gradient<-colorRampPalette(c('darkgoldenrod1','darkgoldenrod4'))

##-------Fig. 1------##
png("Figure 1.png",width=10,height=5,bg="white",res=500,units="in")
par(mfrow=c(1,2))
h<-u<-0
N<-seq(0,parameters['K'],0.1)
prey.cols<-prey.gradient(2)
pred.cols<-pred.gradient(2)
plot(N,prey.isocline.vec(replace(parameters,"K",Inf),N,h,u,"fixed"),type="l",xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",col=prey.cols[1],ylim=c(0,parameters['K']),xlim=c(0,parameters['K']),axes=FALSE,lwd=1.5)
axis(1,pos=0)
axis(2,pos=0)
segments(pred.isocline(replace(parameters,"d",0.1),N,h,u,"fixed"),h,pred.isocline(replace(parameters,"d",0.1),N,h,u,"fixed"),150,col=pred.cols[1],lwd=1.5)
segments(0,0,500,500,lty=3)

h<-u<-0
N<-seq(0,parameters['K'],0.1)
plot(N,prey.isocline.vec(parameters,N,h,u,"fixed"),type="l",xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",col=prey.cols[2],ylim=c(0,parameters['K']),xlim=c(0,parameters['K']),axes=FALSE,lwd=1.5)
axis(1,pos=0)
axis(2,pos=0)
segments(pred.isocline(replace(parameters,"d",0.1),N,h,u,"fixed"),h,pred.isocline(replace(parameters,"d",0.1),N,h,u,"fixed"),150,col=pred.cols[2],lwd=1.5)
segments(0,0,500,500,lty=3)
dev.off()

##-------Fig. 2------##
png("Figure 2.png",width=5,height=5,bg="white",res=500,units="in")
h.range<-seq(0,21.25,4.25)
u<-0
N<-seq(0,parameters['K'],0.1)
prey.cols<-prey.gradient(length(h.range))
pred.cols<-pred.gradient(length(h.range))
plot(N,prey.isocline.vec(parameters,N,h.range[1],u,"fixed"),type="l",xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",col=prey.cols[1],ylim=c(0,parameters['K']),xlim=c(0,parameters['K']),axes=FALSE)
axis(1,pos=0)
axis(2,pos=0)
segments(pred.isocline(parameters,N,h.range[1],u,"fixed"),0,pred.isocline(parameters,N,h.range[1],u,"fixed"),150,col=pred.cols[1])
for(i in seq_along(h.range)[-1]){
  lines(N,prey.isocline(parameters,N,h.range[i],u,"fixed"),col=prey.cols[i])
  segments(pred.isocline(parameters,N,h.range[i],u,"fixed"),0,pred.isocline(parameters,N,h.range[i],u,"fixed"),150,col=pred.cols[i])
}
h.fine<-seq(h.range[1],tail(h.range,1),0.01)
segments(0,0,500,500,lty=3)
lines(N_star(parameters,h.fine,u,"fixed"),P_star(parameters,h.fine,u,"fixed"),lty=4,lwd=1.5)
dev.off()

##-------Fig. 3------##
png("Figure 3.png",width=7.5,height=7.5,bg="white",res=500,units="in")
par(mfrow=c(2,2))
h<-0
u.range<-c(0,0.375,0.75,0.99)
prey.cols<-prey.gradient(length(u.range))
pred.cols<-pred.gradient(length(u.range))
N<-seq(0,parameters['K'],0.1)
for(i in seq_along(u.range)){
  plot(N,prey.isocline.vec(parameters,N,h,u.range[i],"fixed"),type="l",xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",col=prey.cols[i],ylim=c(0,parameters['K']),xlim=c(0,parameters['K']),axes=FALSE,lwd=1.5)
  axis(1,pos=0)
  axis(2,pos=0)
  segments(pred.isocline(parameters,N,h,u.range[i],"fixed"),0,pred.isocline(parameters,N,h,u.range[i],"fixed"),150,col=pred.cols[i],lty=1,lwd=1.5)
  u.fine<-seq(0,u.range[i],0.01)
  segments(0,0,500,500,lty=3)
  lines(N_star(parameters,h,u.fine,"fixed"),P_star(parameters,h,u.fine,"fixed"),lty=4,lwd=1.5)
  points(N_star(parameters,h,tail(u.fine,1),"fixed"),P_star(parameters,h,tail(u.fine,1),"fixed"),pch='*')
}
dev.off()

##-------Fig. 4------##
png("Figure 4.png",width=7.5,height=7.5,bg="white",res=500,units="in")
par(mfrow=c(2,2))
plot.new()
plot.window(xlim=c(0,parameters['K']),ylim=c(0,parameters['K']))
axis(1,pos=0)
axis(2,pos=0)
cont.N<-seq(0,parameters['K'],1)
cont.P<-seq(0,parameters['K'],1)
u.star<-outer(cont.N,cont.P,"u_star.vec",params=parameters)
lev.num=20
cont.cols<-terrain.colors(lev.num+1)
transpar.cont.cols<-add.alpha(cont.cols,alpha=0.5)
.filled.contour(cont.N,cont.P,u.star,levels=seq(0,1,length.out=lev.num),col=cont.cols)
title(xlab='Prey Biomass (N)',ylab='Predator Biomass (P)')

plot.new()
h<-0; u.range<-seq(0,1,0.25)
cols<-color.gradient(length(u.range))
plot.window(xlim=c(0,parameters['K']),ylim=c(0,parameters['K']))
.filled.contour(cont.N,cont.P,u.star,levels=seq(0,1,length.out=lev.num),col=transpar.cont.cols)
lines(N,prey.isocline.vec(parameters,N,h,u.range[1],"fixed"),lty=1,col=rgb(0,0,1,0.5))
# polygon(c(0,0,parameters['K']),c(0,u.isoleg(parameters,0,0,u.range[1]),0),density=20,col=cols[1],angle=135)
# polygon(c(0,parameters['K'],parameters['K'],0),c(u.isoleg(parameters,0,0,tail(u.range,1)),0,1000,1000),density=20,col=tail(cols,1),angle=45)
lines(N,prey.isocline.vec(parameters,N,h,u.range[2],"fixed"),lty=2,col=rgb(0,0,1,1))
axis(1,pos=0)
axis(2,pos=0)
title(xlab="Prey Biomass (N)",ylab="Predator Biomass (P)")

plot.new()
plot.window(xlim=c(0,parameters['K']),ylim=c(0,parameters['K']))
# polygon(c(0,0,parameters['K']),c(0,u.isoleg(parameters,0,0,u.range[1]),0),density=20,col=cols[1],angle=135)
# polygon(c(0,parameters['K'],parameters['K'],0),c(u.isoleg(parameters,0,0,tail(u.range,1)),0,1000,1000),density=20,col=tail(cols,1),angle=45)
axis(1,pos=0)
axis(2,pos=0)
.filled.contour(cont.N,cont.P,u.star,levels=seq(0,1,length.out=lev.num),col=transpar.cont.cols)
segments(pred.isocline.vec(parameters,N[1],h,u[1],"fixed"),0,pred.isocline.vec(parameters,N[1],h,u[1],"fixed"),150,col=rgb(1,0,0,0.5),lty=1)
pred.lines.IBP(parameters,N,h,u[0],"one","optimized")
title(xlab="Prey Biomass (N)",ylab="Predator Biomass (P)")

N<-seq(0,parameters['K'],1)
plot.new()
plot.window(xlim=c(0,parameters['K']),ylim=c(0,parameters['K']))
# polygon(c(0,0,parameters['K']),c(0,u.isoleg(parameters,0,0,u.range[1]),0),density=20,col=cols[1],angle=135)
# polygon(c(0,parameters['K'],parameters['K'],0),c(u.isoleg(parameters,0,0,tail(u.range,1)),0,1000,1000),density=20,col=tail(cols,1),angle=45)
.filled.contour(cont.N,cont.P,u.star,levels=seq(0,1,length.out=lev.num),col=transpar.cont.cols)
lines(N,prey.isocline.vec(parameters,N,h,u[1],"fixed"),col=rgb(0,0,1,0.5))
segments(pred.isocline.vec(parameters,N[1],h,u[1],"fixed"),0,pred.isocline.vec(parameters,N[1],h,u[1],"fixed"),150,col=rgb(1,0,0,0.5),lty=1)
axis(1,pos=0)
axis(2,pos=0)
plot.IBP(parameters,N,h,u[0],"one","optimized",ylim=c(0,parameters['K']))
title(xlab="Prey Biomass (N)",ylab="Predator Biomass (P)")
dev.off()

##-------Fig. 5------##
png("Figure 5.png",width=10,height=5,bg="white",res=500,units="in")
par(mfrow=c(1,2))
b.range<-k.range<-seq(0,200,1)
d.range<-seq(0.03,0.07,0.005)
cols<-gold.gradient(length(d.range))
plot.new()
plot.window(xlim=c(min(b.range),max(b.range)),ylim=c(0,4))
axis(1,pos=0)
axis(2,pos=0)
for(i in seq_along(d.range)){
  lines(b.range,IBP.criterion.vec(replace(parameters,'d',d.range[i]),b.range,parameters['k'],0,0,'optimized'),type='l',col=cols[i],lwd=1.5)
}
abline(h=1,lty=3)
title(ylab=expression(IBP~criterion~~'('*italic(over(P^'*',N^'*'))*')'),line=1.5)
title(xlab=expression(Effectiveness~of~Vigilance~'('*italic(b)*')'))

plot.new()
plot.window(xlim=c(min(k.range),max(k.range)),ylim=c(0,4))
axis(1,pos=0)
axis(2,pos=0)
for(i in seq_along(d.range)){
  lines(k.range,IBP.criterion.vec(replace(parameters,'d',d.range[i]),parameters['b'],k.range,0,0,'optimized'),type='l',col=cols[i],lwd=1.5)
}
abline(h=1,lty=3)
title(ylab=expression(IBP~criterion~~'('*italic(over(P^'*',N^'*'))*')'),line=1.5)
title(xlab=expression(Inverse~Predator~Lethality~'('*italic(k)*')'))
dev.off()
