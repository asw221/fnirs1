
betafe = as.matrix(read.table("sub_betadraws0.dat"))
BASIS = as.matrix(read.table("HRF.dat"))

Xbeta = hrf%*%t(beta[,15:21])
mXbeta = apply(Xbeta,1,mean)
sdXbeta = apply(Xbeta,1,sd)
Zm = apply(abs((Xbeta-mXbeta) /sdXbeta),2,max,na.rm=T)
qa = quantile(Zm,0.95,na.rm=T)
Ia = matrix(0,nrow(Xbeta),2)
Ia[,1] = (mXbeta - qa*sdXbeta)
Ia[,2] = (mXbeta + qa*sdXbeta)
matplot(Ia,type="l",lty=1,col=2)
abline(h=0)
lines(mXbeta,col=1)


pdf("/Users/tdj/seminars/Mexico2016/figures/jointband1.pdf")
xx = (1:nrow(Ia))/20
xx = c(xx,rev(xx))
yy = c(Ia[,1],rev(Ia[,2]))
matplot((1:nrow(Ia))/20,Ia,type="n",lty=1,col=2,main="E(HRF | Y) + Joint Bands, Stimulus 1",ylab="Hemodynamic Response Function",xlab="Seconds")
polygon(xx,yy,col="gray75",border="gray75")
abline(h=0)
lines((1:nrow(Ia))/20,mXbeta,col=1,lwd=2)
lines((1:nrow(Ia))/20,mXbeta,col=1,lwd=2)
segments(17.95,-0.001,22.55,-0.001,col=1,lwd=3)
text(20.25,-.002,"4.6 Sec")
#segments(17.75,-0.0005,22.90,-0.0005,col=1,lwd=3)
#text(20.325,-.001,"5.15 Sec")
dev.off()





Xbeta = BASIS%*%t(betafe[,8:14])
mXbeta = apply(Xbeta,1,mean)
sdXbeta = apply(Xbeta,1,sd)
Zm = apply(abs((Xbeta-mXbeta) /sdXbeta),2,max,na.rm=T)
qa = quantile(Zm,0.95,na.rm=T)
Ia = matrix(0,nrow(Xbeta),2)
Ia[,1] = (mXbeta - qa*sdXbeta)
Ia[,2] = (mXbeta + qa*sdXbeta)
matlines(Ia,type="l",lty=1,col=3)
abline(h=0)
lines(mXbeta,col=4,lty=1)
