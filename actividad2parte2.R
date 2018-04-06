library(SpiecEasi)
library(Matrix)
library(reshape2)
library(igraph)
library(mDINGO)
library(metagenomeSeq)

p=read.table(file="p1.txt",header = TRUE)
n=read.table(file="n1.txt",sep="\t",header = TRUE)
Y_dis = read.table(file="otudi1.txt",header = TRUE)
Y_con = read.table(file="otucon1.txt",header = TRUE)
YTotal = read.table(file="otutotal1.txt",header = TRUE)
X = YTotal[1:578,2:128]



se.mb<-spiec.easi(Y_dis,method='mb',lambda.min.ratio=1e-2,nlambda=20,
                  icov.select.params=list(rep.num=20,ncores=2))
se.gl<-spiec.easi(Y_dis,method='glasso',lambda.min.ratio=1e-2,nlambda=20,
                  icov.select.params=list(rep.num=20,ncores=2))
sp.est<-sparcc(Y_dis)
## define arbitrary correlation threshold, see also sparccboot
sp.net<-abs(sp.est$Cor)>=0.15
diag(sp.net)<-0

ig.mb<-adj2igraph(symBeta(getOptBeta(se.mb),mode='maxabs'))
ig.gl<-adj2igraph(se.gl$opt.cov*se.gl$refit)
ig.sp<-adj2igraph(sp.est$Cor*sp.net)


## uncommenct for side-by-side plotting
# par(mfrow=c(2,2))
#plotnet(ig.mb,main="MB")
vsize <- rowMeans(clr(amgut1.filt, 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB-Dissease")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso-Dissease")
plot(ig.sp, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc-Dissease")



#dd.gl <- degree.distribution(ig.gl)
#dd.mb <- degree.distribution(ig.mb)
#dd.sparcc <- degree.distribution(ig.sp)

#plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,.5), type='b', 
#     ylab="Frequency", xlab="Degree", main="Degree Distributions")
#points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
#points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
#legend("topright", c("MB-Control", "glasso-Control", "sparcc-Control"), 
#       col=c("forestgreen", "red", "black"), pch=1, lty=1)






#se.mb<-spiec.easi(Y_con,method='mb',lambda.min.ratio=1e-2,nlambda=20,
#                  icov.select.params=list(rep.num=20,ncores=2))
#se.gl<-spiec.easi(Y_con,method='glasso',lambda.min.ratio=1e-2,nlambda=20,
 #                 icov.select.params=list(rep.num=20,ncores=2))
#sp.est<-sparcc(Y_con)
## define arbitrary correlation threshold, see also sparccboot
#sp.net<-abs(sp.est$Cor)>=0.15
#diag(sp.net)<-0

#ig.mb<-adj2igraph(symBeta(getOptBeta(se.mb),mode='maxabs'))
#ig.gl<-adj2igraph(se.gl$opt.cov*se.gl$refit)
#ig.sp<-adj2igraph(sp.est$Cor*sp.net)


## uncommenct for side-by-side plotting
# par(mfrow=c(2,2))
#plotnet(ig.mb,main="MB")
#vsize <- rowMeans(clr(amgut1.filt, 1))+6
#am.coord <- layout.fruchterman.reingold(ig.mb)

#par(mfrow=c(1,3))
#plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB-Control")
#plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso-Control")
#plot(ig.sp, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc-Control")