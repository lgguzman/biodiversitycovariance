
library(SpiecEasi)
library(Matrix)
library(reshape2)
library(igraph)
library(mDINGO)
library(metagenomeSeq)
#availcores = 1#parallel::detectCores() - 1
data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1 <- round(amgut1.filt.n * min(depths))
p = dim(amgut1)[2]
p
n = dim(amgut1)[1]
graphList = gen_graphs_sf(p, 1, .6, .05)
write.table(p,file="p3.txt",sep="\t",col.names=NA)
#write.table(n,file="n3.txt",sep="\t",col.names=NA)
#l_fr = layout_with_fr(graphList$disease)
par(mfrow = c(2,2), mar=c(.5,.5,1,.5))
plot(graphList$disease, vertex.size= 0, main = "Disease Graph", layout = l_fr, vertex.label = NA)
plot(graphList$control, vertex.size= 0, main = "Control Graph", layout = l_fr, vertex.label = NA)


prec_dis = igraph2prec(graphList$disease)
cor_dis = cov2cor(prec2cov(prec_dis))
prec_con = igraph2prec(graphList$control)
cor_con = cov2cor(prec2cov(prec_con))
Y_dis = synth_comm_from_counts(amgut1, mar=2, distr='zinegbin', Sigma=cor_dis, n=n)
Y_con = synth_comm_from_counts(amgut1, mar=2, distr='zinegbin', Sigma=cor_con, n=n)
write.table(Y_dis,file="otudi3.txt",sep="\t",col.names=NA)
write.table(Y_con,file="otucon3.txt",sep="\t",col.names=NA)
YTotal = rbind(Y_dis, Y_con)
colnames(YTotal) = paste0("OTU", 1:dim(YTotal)[2])
write.table(YTotal,file="otutotal3.txt",sep="\t",col.names=NA) 




p=read.table(file="p3.txt",header = TRUE)
n=read.table(file="n3.txt",sep="\t",header = TRUE)
Y_dis = read.table(file="otudi3.txt",header = TRUE)
Y_con = read.table(file="otucon3.txt",header = TRUE)
YTotal = read.table(file="otutotal3.txt",header = TRUE)
YTotal = YTotal[1:578,2:128]
#normFactor = normFactors(as.data.frame(Y_dis))
disease = as.data.frame.matrix(t(YTotal))
Y<-newMRexperiment(disease,libSize=colSums(data.matrix(disease)),normFactors=1)
#lungp = cumNormStat(Y, pFlag = TRUE, main = "Data")
#Y = cumNorm(Y, p = lungp)
metadata<- c(rep("a",289), rep("b", 289))
#metadata <- disease[1,1:127]
metadata<-data.matrix(metadata) 

X<-model.matrix(~ metadata, data= disease)
#disease = newMRexperiment(Y_dis) 
#disease <- cumNorm(disease, p=.5)
#pdataDisease= pData(disease) 
#dataframe = as.data.frame(Y_dis)
#mod = model.matrix( ~ V1 , data=pdataDisease)
settings = zigControl(maxit = 1000, verbose = TRUE)
#fit = fitZig(obj = disease, mod = mod, useCSSoffset = FALSE,
#             control = settings)
#settings = zigControl(per_feature_zeroModel=TRUE)
fit = fitZig(obj = Y,mod=X,control=settings,useCSSoffset = FALSE)
final=MRcoefs(fit, group = 3)
fit2  = fitFeatureModel(Y, X)
final2=MRcoefs(fit2, group = 3)
