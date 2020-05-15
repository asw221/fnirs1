 la = c(2,3,6,10,1,5,8,4,7,9)
 lb = c(3,4,7,2,5,9,1,6,8,10)
 lc = c(1,6,7,2,3,5,8,4,9,10)
 
 a = read.table("design0856a.dat")
 write.table(a[,la],"design0856a.dat",row.names=F,col.names=F)

 a = read.table("design0856b.dat")
 write.table(a[,lb],"design0856b.dat",row.names=F,col.names=F)

 a = read.table("design0856c.dat")
 write.table(a[,lc],"design0856c.dat",row.names=F,col.names=F)

 a = read.table("design1641a.dat")
 write.table(a[,la],"design1641a.dat",row.names=F,col.names=F)

 a = read.table("design1641b.dat")
 write.table(a[,lb],"design1641b.dat",row.names=F,col.names=F)

 a = read.table("design1641c.dat")
 write.table(a[,lc],"design1641c.dat",row.names=F,col.names=F)

 a = read.table("design1659a.dat")
 write.table(a[,la],"design1659a.dat",row.names=F,col.names=F)

 a = read.table("design1659b.dat")
 write.table(a[,lb],"design1659b.dat",row.names=F,col.names=F)

 a = read.table("design1659c.dat")
 write.table(a[,lc],"design1659c.dat",row.names=F,col.names=F)

 a = read.table("design1706a.dat")
 write.table(a[,la],"design1706a.dat",row.names=F,col.names=F)

 a = read.table("design1706b.dat")
 write.table(a[,lb],"design1706b.dat",row.names=F,col.names=F)

 a = read.table("design1706c.dat")
 write.table(a[,lc],"design1706c.dat",row.names=F,col.names=F)

 a = read.table("design1857a.dat")
 write.table(a[,la],"design1857a.dat",row.names=F,col.names=F)

 a = read.table("design1857b.dat")
 write.table(a[,lb],"design1857b.dat",row.names=F,col.names=F)

 a = read.table("design1857c.dat")
 write.table(a[,lc],"design1857c.dat",row.names=F,col.names=F)

#######################3333333


 a = read.table("design0856a.dat")
 a = cbind(apply(a[,1:4],1,sum),apply(a[,5:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design0856a.dat",row.names=F,col.names=F)

 a = read.table("design0856b.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:6],1,sum),apply(a[,7:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design0856b.dat",row.names=F,col.names=F)

 a = read.table("design0856c.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design0856c.dat",row.names=F,col.names=F)

 a = read.table("design1641a.dat")
 a = cbind(apply(a[,1:4],1,sum),apply(a[,5:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1641a.dat",row.names=F,col.names=F)

 a = read.table("design1641b.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:6],1,sum),apply(a[,7:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1641b.dat",row.names=F,col.names=F)

 a = read.table("design1641c.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1641c.dat",row.names=F,col.names=F)

 a = read.table("design1659a.dat")
 a = cbind(apply(a[,1:4],1,sum),apply(a[,5:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1659a.dat",row.names=F,col.names=F)

 a = read.table("design1659b.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:6],1,sum),apply(a[,7:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1659b.dat",row.names=F,col.names=F)

 a = read.table("design1659c.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1659c.dat",row.names=F,col.names=F)

 a = read.table("design1706a.dat")
 a = cbind(apply(a[,1:4],1,sum),apply(a[,5:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1706a.dat",row.names=F,col.names=F)

 a = read.table("design1706b.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:6],1,sum),apply(a[,7:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1706b.dat",row.names=F,col.names=F)

 a = read.table("design1706c.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1706c.dat",row.names=F,col.names=F)

 a = read.table("design1857a.dat")
 a = cbind(apply(a[,1:4],1,sum),apply(a[,5:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1857a.dat",row.names=F,col.names=F)

 a = read.table("design1857b.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:6],1,sum),apply(a[,7:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1857b.dat",row.names=F,col.names=F)

 a = read.table("design1857c.dat")
 a = cbind(apply(a[,1:3],1,sum),apply(a[,4:7],1,sum),apply(a[,8:10],1,sum))
# a = a[seq(from=1,by=5,to=nrow(a)),]
 write.table(a,"design1857c.dat",row.names=F,col.names=F)


a = scan("Y0856_12a.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y0856_12a.dat")

a = scan("Y0856_12b.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y0856_12b.dat")

a = scan("Y0856_12c.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y0856_12c.dat")

a = scan("Y1641_12a.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1641_12a.dat")

a = scan("Y1641_12b.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1641_12b.dat")

a = scan("Y1641_12c.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1641_12c.dat")

a = scan("Y1659_12a.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1659_12a.dat")

a = scan("Y1659_12b.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1659_12b.dat")

a = scan("Y1659_12c.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1659_12c.dat")

a = scan("Y1706_12a.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1706_12a.dat")

a = scan("Y1706_12b.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1706_12b.dat")

a = scan("Y1706_12c.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1706_12c.dat")

a = scan("Y1857_12a.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1857_12a.dat")

a = scan("Y1857_12b.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1857_12b.dat")

a = scan("Y1857_12c.dat")
a = a[seq(1,by=5,to=length(a))]
write(a,"Y1857_12c.dat")


