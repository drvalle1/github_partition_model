rm(list=ls(all=TRUE))
set.seed(6)

setwd('U:\\independent studies\\partition models\\no covariates\\github_partition_model')
source('gibbs sampler functions.R')

setwd('U:\\independent studies\\partition models\\no covariates')
dat=read.csv('fake data.csv',as.is=T)

max.regions=100
n.regions=50
uni.loc=unique(dat[,c('loc.id','LATNUM','LONGNUM')])
uni.loc=uni.loc[order(uni.loc$loc.id),]
nloc=nrow(uni.loc)

dist.mat=data.matrix(dist(uni.loc[,c('LONGNUM','LATNUM')]))
# x2=(uni.loc$LONGNUM[2]-uni.loc$LONGNUM[1:5])^2
# y2=(uni.loc$LATNUM[2] -uni.loc$LATNUM[1:5])^2
# sqrt(x2+y2)
# dist.mat[2,1:5]

centroid=sort(sample(1:nloc,size=n.regions))
tmp=assign.centroid(centroid,dist.mat)
region.map=data.frame(loc.id=1:nloc,centroid=tmp)
#---------------------------------------
#priors

phi.a=phi.b=1
#---------------------------------------
ngibbs=10000
vec.logl=matrix(NA,ngibbs,1)
vec.n=matrix(NA,ngibbs,1)
vec.centroids=matrix(NA,ngibbs,max.regions)
vec.phi=matrix(NA,ngibbs,max.regions)

param=list(region.map=region.map)
max.logl=-Inf

for (i in 1:ngibbs){
  #how many distinct regions?
  n=length(centroid)
  print(c(i,n))
  
  tmp=update.regions(param)
  param$region.map=tmp$region
  
  vec.logl[i]=tmp$prob
  vec.n[i]=n
  
  ind=1:length(tmp$phi)
  centroid=sort(unique(param$region.map$centroid))
  vec.centroids[i,ind]=centroid
  vec.phi[i,ind]=tmp$phi
  
  if (tmp$prob>max.logl) {
    max.logl=tmp$prob
    fim=param$region.map
  }
}

seq1=2000:ngibbs
plot(vec.logl[seq1],type='l')

setwd('U:\\independent studies\\partition models\\no covariates\\results')
write.csv(vec.centroids[seq1,],'centroids.csv',row.names=F)
write.csv(vec.phi[seq1,],'phi.csv',row.names=F)
