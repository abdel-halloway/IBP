rm(list=ls(all=TRUE))
graphics.off()

this.dir<-dirname(parent.frame(2)$ofile)
work.dir<-paste0(this.dir,'/Figures/')
dir.create(work.dir,recursive=TRUE)
setwd(work.dir)

prey.isocline<-function(params,N.,h.,u.,vigilance){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  u<-u.; h<-h.
  
  a<-m/(k+b*u)
  
  if(vigilance=="r"){
    x<-r/a*(1-u)*(1+a*h*N.)*(1-N./K)
  }else if(vigilance=="optimized"){
    x<-r/(4*m*b)*(b+k+m*h*N.)^2*(1-N./K)
  }else if(vigilance=="K"){
  	x<-r/a*(1+a*h*N.)*(1-u-N./K)
  }

  return(x)
}

pred.isocline<-function(params,N.,h.,u.,vigilance){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  u<-u.; h<-h.
  
  a<-m/(k+b*u)
  
  if(vigilance=="r" | vigilance=="K"){
    x<-d/(a*(c-d*h))
  }else if(vigilance=="optimized"){
    x<-r*m/b*(c/d*N.)^2*(1-N./K)
  }

  return(x)
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
  
  if(vigilance=="r"){
    P_star<-r/a*(1-u)*(1+d*h/(c-d*h))*(1-d/(a*K*(c-d*h)))
  }else if(vigilance=="optimized"){
    P_star<-r/(4*m*b)*(1-d*(b+k)/(m*K)/(2*c-d*h))*(b+k+d*h*(b+k)/(2*c-d*h))^2
  }else if(vigilance=="K"){
  	P_star<-r/a*(1+d*h/(c-d*h))*(1-u-d/(a*K*(c-d*h)))
	}
  
  return(P_star)
}

prey.isocline.vec<-Vectorize(prey.isocline,c("N.","h.","u."))
pred.isocline.vec<-Vectorize(pred.isocline,c("N.","h.","u."))
P_star.vec<-Vectorize(P_star,c("h.","u."))

parameters<-c(r=0.08,K=100,d=0.04,c=1,m=0.04,k=15,b=60)
prey.gradient<-colorRampPalette(c('blue1','blue4'))
pred.gradient<-colorRampPalette(c('firebrick1','firebrick4'))

pdf("Figure 1.pdf",width=10,height=10,bg="white")
par(mfrow=c(2,2),cex=1.25)
h<-0
u.range<-c(0,0.375,0.80,0.99)
N<-seq(0,parameters['K'],0.1)
prey.cols<-prey.gradient(length(u.range))
pred.cols<-pred.gradient(length(u.range))
prey.iso<-prey.isocline.vec(parameters,N,h,u.range[1],"r")
prey.iso[which(prey.iso<0)]=NA
plot(N,prey.iso,type="l",xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",
     col=prey.cols[1],ylim=c(0,parameters['K']),xlim=c(0,parameters['K']),axes=FALSE,lwd=2.5)
axis(1,pos=0)
axis(2,pos=0)
segments(pred.isocline(parameters,N,h,u.range[1],"r"),0,pred.isocline(parameters,N,h,u.range[1],"r"),150,col=pred.cols[1],lwd=2.5)
for(i in seq_along(u.range)[-1]){
  prey.iso<-prey.isocline.vec(parameters,N,h,u.range[i],"r")
  prey.iso[which(prey.iso<0)]=NA
  lines(N,prey.iso,col=prey.cols[i],lwd=2.5)
  segments(pred.isocline(parameters,N,h,u.range[i],"r"),0,pred.isocline(parameters,N,h,u.range[i],"r"),150,col=pred.cols[i],lwd=2.5)
}
u.fine<-seq(u.range[1],tail(u.range,1),0.01)

P.star<-P_star(parameters,h,u.fine,"r")
P.star[which(P.star<0)]=0

segments(0,0,500,500,lty=3,lwd=2.5)
lines(N_star(parameters,h,u.fine,"fixed"),P.star,lty=4,lwd=3)
mtext('A)',side=3, adj=0, cex=2, line=1)

plot(u.fine,P.star/N_star(parameters,h,u.fine,"fixed"),xlab="",ylab="",type='l',lwd=2.5)
title(ylab='Predator-Prey Biomass Ratio')
title(xlab=expression(Vigilance~'('*italic(u)*')'))
abline(h=1,lty=3)
mtext('B)',side=3, adj=0, cex=2, line=1)

parameters<-c(r=0.30,K=200,d=0.1,c=1,m=0.04,k=15,b=60,h=0,u=0)
h<-0
u.range<-c(0,0.10,6.5/14,0.93)
N<-seq(0,parameters['K'],0.1)
prey.cols<-prey.gradient(length(u.range))
pred.cols<-pred.gradient(length(u.range))
prey.iso<-prey.isocline.vec(parameters,N,h,u.range[1],"K")
prey.iso[which(prey.iso<0)]=NA
plot(N,prey.iso,type="l",xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",
	 col=prey.cols[1],ylim=c(0,parameters['K']),xlim=c(0,parameters['K']),axes=FALSE,lwd=2.5)
axis(1,pos=0)
axis(2,pos=0)
segments(pred.isocline(parameters,N,h,u.range[1],"K"),0,pred.isocline(parameters,N,h,u.range[1],"K"),250,col=pred.cols[1],lwd=2.5)
for(i in seq_along(u.range)[-1]){
	prey.iso<-prey.isocline.vec(parameters,N,h,u.range[i],"K")
	prey.iso[which(prey.iso<0)]=NA
	lines(N,prey.iso,col=prey.cols[i],lwd=2.5)
	segments(pred.isocline(parameters,N,h,u.range[i],"K"),0,pred.isocline(parameters,N,h,u.range[i],"K"),250,col=pred.cols[i],lwd=2.5)
}
u.fine<-seq(u.range[1],1,0.01)

P.star<-P_star(parameters,h,u.fine,"K")
P.star[which(P.star<0)]=0

segments(0,0,500,500,lty=3,lwd=2.5)
lines(N_star(parameters,h,u.fine,"fixed"),P.star,lty=4,lwd=3)
mtext('C)',side=3, adj=0, cex=2, line=1)

plot(u.fine,P.star/N_star(parameters,h,u.fine,"fixed"),xlab="",ylab="",type='l',lwd=2.5)
title(ylab='Predator-Prey Biomass Ratio')
title(xlab=expression(Vigilance~'('*italic(u)*')'))
abline(h=1,lty=3)
mtext('D)',side=3, adj=0, cex=2, line=1)
dev.off()