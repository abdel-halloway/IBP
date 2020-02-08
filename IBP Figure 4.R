rm(list=ls(all=TRUE))
graphics.off()

this.dir<-dirname(parent.frame(2)$ofile)
work.dir<-paste0(this.dir,'/Figures/')
dir.create(work.dir,recursive=TRUE)
setwd(work.dir)

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

prey.isocline<-function(params,N.,vigilance){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  u<-params['u']; h<-params['h']
  
  a<-m/(k+b*u)
  
  if(vigilance=="fixed"){
    x<-r/a*(1+a*h*N.)*(1-u-N./K)
  }else if(vigilance=="optimized"){
    x<-ifelse(0>=(1-k/b-N./K),r/a*(1-N./K),x<-r*b/(4*m)*(1-N./K+k/b)^2)
  }
  
  return(x)
}

pred.isocline<-function(params,N.,vigilance){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']
  u<-params['u']; h<-params['h']
  
  a<-m/(k+b*u)
  
  if(vigilance=="fixed"){
    x<-d/(a*(c-d*h))
  }else if(vigilance=="optimized"){
    x<-r*m/b*(c/d*N.)^2
  }
  
  return(x)
}

u_star<-function(N.,P.,params){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']; h<-params['h']
  
  u_star<-sqrt(m*P./(r*b))-(k+m*h*N.)/b
  
  if(u_star<0 | is.nan(u_star) | is.na(u_star)){
    u_star<-0
  }
  if(u_star>1 | is.infinite(u_star)){
    u_star<-1
  }
  
  return(u_star)
}

plot.IBP<-function(params,N.,type,vigilance,...){
  if(vigilance=="optimized"){
    pred.no.vig.params<-replace(params,"u",0); pred.total.vig.params<-replace(params,"u",1)
    pred.N<-seq(pred.isocline(pred.no.vig.params,N.,"fixed"),pred.isocline(pred.total.vig.params,N.,"fixed"),0.01)
    prey.iso<-prey.isocline(params,N.,vigilance)
    pred.iso<-pred.isocline(params,pred.N,vigilance)
    lines(N.,prey.iso,lty=1,col="blue",...)
    if(type=="one"){
      N_one_star<-expression(frac(dk,cm))
      N_two_star<-expression(frac(d(k+b),cm))
    }else if(type=="two"){
      N_one_star<-expression(frac(dk,m(c-dh)))
      N_two_star<-expression(frac(d(k+bu),m(c-dh)))
    }
    segments(pred.N[1],0,pred.N[1],pred.iso[1],col="red",lty=1,...)
    segments(tail(pred.N,1),tail(pred.iso,1),tail(pred.N,1),1e9,col="red",lty=1,...)
    lines(pred.N,pred.iso,col="red",lty=1,...)
    segments(0,0,500,500,lty=3,...)
  }else{
    prey.iso<-prey.isocline(params,N.,vigilance)
    pred.iso<-pred.isocline(params,N.,vigilance)
    
    plot(N.,prey.iso,type='l',xlim=c(0,params['K']),xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",axes=FALSE,col="darkblue",...)
    abline(v=pred.iso,lty=1,col="darkred",...)
    abline(a=0,b=1,lty=3,...)
  }
}

G_N<-function(params,v.,N.,P.){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']; h<-params['h']
  
  r*(1-v.-N./K)-m*P./(k+b*v.+m*h*N.)
}

G_P<-function(params,v.,N.,P.){
  r<-params['r']; K<-params['K']
  d<-params['d']; c<-params['c']
  m<-params['m']; k<-params['k']; b<-params['b']; h<-params['h']
  
  c*m*N./(k+b*v.+m*h*N.)-d
}

pdf("Figure 4.pdf",bg="white",width=10,height=5*1.05)
par(mar=c(5.1,4.1,2.1,1.1))
m<-matrix(c(1,1,2,2,3,
            1,1,2,2,3
            ),ncol=5,byrow=TRUE)
layout(mat=m, heights=c(0.225,0.225),widths=c(0.225,0.225,0.25,0.25,0.1))

parameters<-c(r=0.08,K=200,d=0.04,c=1,m=0.04,k=15,b=60,h=0,u=0)
Timesteps<-1e4
InitN<-1
InitP<-1
Initu<-u_star(InitN,InitP,parameters)

dt<-0.1

u.dynam<-N.dynam<-P.dynam<-numeric(length=Timesteps+1)

u.dynam[1]<-Initu; N.dynam[1]<-InitN; P.dynam[1]<-InitP

for(t in 1:Timesteps){
  N.dynam[t+1]<-N.dynam[t]*exp(G_N(parameters,u.dynam[t],N.dynam[t],P.dynam[t])*dt)
  P.dynam[t+1]<-P.dynam[t]*exp(G_P(parameters,u.dynam[t],N.dynam[t],P.dynam[t])*dt)
  u.dynam[t+1]<-u_star(N.dynam[t+1],P.dynam[t+1],parameters)
}

par(cex=1)
prey.isocline.vec<-Vectorize(prey.isocline,c("N."))
pred.isocline.vec<-Vectorize(pred.isocline,c("N."))
u_star.vec<-Vectorize(u_star,c("N.","P."))

cont.N<-seq(0,parameters['K'],0.2)
cont.P<-seq(0,parameters['K'],0.2)
u.star<-outer(cont.N,cont.P,"u_star.vec",params=parameters)
lev.num=12
cont.cols<-terrain.colors(lev.num)
transpar.cont.cols<-add.alpha(cont.cols,alpha=0.5)
col.breaks<--1:11/10
col.breaks[2]<-1e-15; col.breaks[12]<-1-1e-15

N<-seq(0,parameters['K'],1)

plot.window(xlim=c(0,parameters['K']),ylim=c(0,parameters['K']))
image(cont.N,cont.P,u.star,col=transpar.cont.cols,xlab="Prey Biomass (N)",ylab="Predator Biomass (P)",breaks=col.breaks)
lines(N,prey.isocline.vec(parameters,N,"fixed"),col=rgb(0,0,1,0.5),lty=2,lwd=2)
segments(pred.isocline.vec(parameters,N[1],"fixed"),0,pred.isocline.vec(parameters,N[1],"fixed"),250,col=rgb(1,0,0,0.5),lty=2,lwd=2)
plot.IBP(parameters,N,"one","optimized",ylim=c(0,parameters['K']),lwd=2)
x.incs<-y.incs<-11
x.seq<-seq(0,parameters['K'],length.out=x.incs); y.seq<-seq(0,parameters['K'],length.out=y.incs)
x.dir<-y.dir<-matrix(nrow=x.incs,ncol=y.incs)
box.length<-min(diff(x.seq),diff(y.seq))
for(i in seq_along(x.seq)){
	for(j in seq_along(y.seq)){
		x<-x.seq[i]; y<-y.seq[j]
		u<-u_star(x,y,parameters)
		x.dir[i,j]<-G_N(parameters,u,x,y); y.dir[i,j]<-arrow.length<-G_P(parameters,u,x,y)
		arrows(x,y,x+x.dir[i,j]*30,y+y.dir[i,j]*30,length = 0.025)
	}
}
mtext('A)',side=3, adj=0, cex=1.5, line=0.5)

par(mar=c(5.1,4.1,2.1,4.1))
plot(u.dynam,type='l',axes=FALSE,lwd=2.5,xlab="Time",ylab="",ylim=c(0,1))
axis(1,pos=0, at=0:5*2e3,labels = c("0","2e3","4e3","6e3","8e3","1e4"))
axis(4,pos=Timesteps)
mtext('B)',side=3, adj=0, cex=1.5, line=0.5)
mtext(expression('Optimal vigilance ('*italic(u)^'*'*')'),side=4,line=2.1)
par(new=TRUE)
plot(P.dynam/N.dynam,type='l',axes=FALSE,lwd=2.5,xlab="",ylab="Predator-Prey Biomass Ratio",col="purple")
axis(2,pos=0)

par(mar=c(1.1,0.1,1.1,3.1),cex.axis=0.75)
x<-seq(0,10)
y<-seq(-0.1,1.1,0.01)
z<-matrix(y,nrow=length(x),ncol=length(y),byrow=TRUE)
image(x,y,z,col=transpar.cont.cols,axes=FALSE,xlab='',ylab='',breaks=col.breaks)
box()
axis(4,pos=10.5,at=0:11/10-0.05,labels=c("0",paste0(0:9/10,'-',1:10/10),'1'),las=2)
dev.off()