Y = scan("mean_Y00.dat")
X = as.matrix(read.table("X00.dat"))
W = as.matrix(read.table("W00.dat"))
V = as.matrix(read.table("V00.dat"))
des = as.matrix(read.table("design3934.dat"))

eta = as.matrix(read.table("etadraws00.dat"))
beta = as.matrix(read.table("betadraws00.dat"))
delta = as.matrix(read.table("deltadraws00.dat"))
fit = scan("mean_fit00.dat")
res = scan("mean_res00.dat")
Wdelta = scan("Wdelta00.dat")
Z = scan("Z00.dat")
x11()
P  = 150
lower = P+1
Y = scan("mean_Y00.dat")
X = as.matrix(read.table("X00.dat"))
W = as.matrix(read.table("W00.dat"))
V = as.matrix(read.table("V00.dat"))
des = as.matrix(read.table("design3934.dat"))

eta = as.matrix(read.table("etadraws00.dat"))
beta = as.matrix(read.table("betadraws00.dat"))
delta = as.matrix(read.table("deltadraws00.dat"))
fit = scan("mean_fit00.dat")
res = scan("mean_res00.dat")
Wdelta = scan("Wdelta00.dat")
Z = scan("Z00.dat")
X11()
P  = 125
lower = P+1
upper = length(Y)
plot(Y,type="l",col="gray75")
lines(1:length(Y),Wdelta,col=7,lwd=1)
lines(1:length(Y),fit,col=4,lwd=1)
lines(lower:length(Y),Z,col=5,lwd=1)
lines(1:length(Y),X%*%apply(beta,2,mean),col=2,lwd=2)
lines(V%*%apply(eta,2,mean),col=5,lwd=2)
matlines(des-1.5,lty=1,col=c(rep(1,4),rep(2,3),rep(3,3)))

plot(apply(delta,2,mean),pch=19,cex=0.25);abline(h=0)

plot(fit,res,pch=19,cex=0.25)

res = res[126:length(res)]
std.res = (res-mean(res))/sd(res)
qqnorm(std.res);abline(0,1,col=2,lwd=2)

upper = length(Y)
plot(Y,type="l",col="gray75")
lines(1:length(Y),Wdelta,col=7,lwd=1)
lines(1:length(Y),fit,col=4,lwd=1)
lines(1:length(Y),X%*%apply(beta,2,mean),col=2,lwd=2)
lines(V%*%apply(eta,2,mean),col=5,lwd=2)
matlines(des-1.5,lty=1,col=c(rep(1,4),rep(2,3),rep(3,3)))

plot(apply(delta,2,mean),pch=19,cex=0.25);abline(h=0)

plot(fit,res,pch=19,cex=0.25)


std.res = (res-mean(res))/sd(res)
qqnorm(std.res);abline(0,1,col=2,lwd=2)


plot((2*(Mod(fft(res)/sqrt(upper-P))[1:((upper-P)/2)])^2),type="h")



ch = c(12,14,15,17,23,25,27,29)

#1857
nn = 1857
a = 1:23676
b = 23701:46400
c = 46402:70670
tmp = readMat("20160328_1857.mat")
for (i in 1:length(ch)) {
write(tmp$HbO[a,ch[i]],paste("Y",nn,"_",ch[i],"a.dat",sep=""))
write(tmp$HbO[b,ch[i]],paste("Y",nn,"_",ch[i],"b.dat",sep=""))
write(tmp$HbO[c,ch[i]],paste("Y",nn,"_",ch[i],"c.dat",sep=""))
}

#1641
nn = 1641
a = 100:23900
b = 23950:45741
c = 45756:68900
tmp = readMat("20160211_1641.mat")
for (i in 1:length(ch)) {
write(tmp$HbO[a,ch[i]],paste("Y",nn,"_",ch[i],"a.dat",sep=""))
write(tmp$HbO[b,ch[i]],paste("Y",nn,"_",ch[i],"b.dat",sep=""))
write(tmp$HbO[c,ch[i]],paste("Y",nn,"_",ch[i],"c.dat",sep=""))
}

#1706
nn = 1706
a = 1:22592
b =22602:45728
c = 46001:68550
tmp = readMat("20151123_1706.mat")
for (i in 1:length(ch)) {
write(tmp$HbO[a,ch[i]],paste("Y",nn,"_",ch[i],"a.dat",sep=""))
write(tmp$HbO[b,ch[i]],paste("Y",nn,"_",ch[i],"b.dat",sep=""))
write(tmp$HbO[c,ch[i]],paste("Y",nn,"_",ch[i],"c.dat",sep=""))
}

#0856
nn = 856
a = 1:21598
b = 21607:43700
c = 43715:66100
tmp = readMat("20160419_0856.mat")
for (i in 1:length(ch)) {
write(tmp$HbO[a,ch[i]],paste("Y0",nn,"_",ch[i],"a.dat",sep=""))
write(tmp$HbO[b,ch[i]],paste("Y0",nn,"_",ch[i],"b.dat",sep=""))
write(tmp$HbO[c,ch[i]],paste("Y0",nn,"_",ch[i],"c.dat",sep=""))
}

#1659
nn = 1659
a = 1:23050
b = 23051:45300
c = 45500:67350
tmp = readMat("20160330_1659.mat")
for (i in 1:length(ch)) {
write(tmp$HbO[a,ch[i]],paste("Y",nn,"_",ch[i],"a.dat",sep=""))
write(tmp$HbO[b,ch[i]],paste("Y",nn,"_",ch[i],"b.dat",sep=""))
write(tmp$HbO[c,ch[i]],paste("Y",nn,"_",ch[i],"c.dat",sep=""))
}


