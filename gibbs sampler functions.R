assign.centroid=function(centroid,dist.mat){
  dist1=dist.mat[,centroid]
  res=rep(NA,nloc)
  for (i in 1:nloc){
    target=dist1[i,]
    ind=which(target==min(target))
    res[i]=centroid[ind]
  }
  res
}
#--------------------
get.likelihood=function(dat,region){
  dat1=merge(dat,region,all=T); dim(dat); dim(dat1)
  dat1$n=1
  dat2=aggregate(cbind(microsc1,n)~centroid,data=dat1,sum)
  p1=lgamma(dat2$microsc1+phi.a)+
     lgamma(dat2$n-dat2$microsc1+phi.b)-
     lgamma(phi.a+dat2$n+phi.b)
  phi=rbeta(nrow(dat2),dat2$microsc1+phi.a,dat2$n-dat2$microsc1+phi.b)
  list(logl=sum(p1),phi=phi)
}
#--------------------
update.regions=function(param){
  region.o=param$region.map
  centroid.o=sort(unique(region.o$centroid))
  nreg=length(centroid.o)

  if (nreg==1){
    centroid.n=birth(centroid.o)
    pdeath=(1/3)
    pbirth=1
    lpjump=log(pdeath/pbirth)
  }
  if (nreg==max.regions){
    centroid.n=death(centroid.o)
    pdeath=1
    pbirth=(1/3)
    lpjump=log(pbirth/pdeath)
  }
  if (nreg!=1 & nreg!=max.regions){
    rv=runif(1)
    if (rv < 1/3){
      centroid.n=birth(centroid.o)
      pdeath=(1/3)
      if (nreg+1==max.regions) pdeath=1
      pbirth=(1/3)
      lpjump=log(pdeath/pbirth)
    }
    if (rv > 1/3 & rv < 2/3){
      centroid.n=death(centroid.o)
      pdeath=(1/3)
      pbirth=1/3
      if (nreg-1==1) pbirth=1
      lpjump=log(pbirth/pdeath)
    }
    if (rv > 2/3){
      centroid.n=swap(centroid.o)
      lpjump=0
    }
  }
  tmp=assign.centroid(centroid.n,dist.mat)
  region.n=data.frame(loc.id=1:nloc,centroid=tmp)
  
  tmp=get.likelihood(dat,region.o)
  pold=tmp$logl
  phi.old=tmp$phi
  
  tmp=get.likelihood(dat,region.n)
  pnew=tmp$logl
  phi.new=tmp$phi
  
  prob=exp(pnew+lpjump-pold)
  rv=runif(1)
  if (rv<prob) return(list(region=region.n,prob=pnew,phi=phi.new))
  return(list(region=region.o,prob=pold,phi=phi.old))
}
#----------------------------------------
death=function(indinz){
  sort(sample(indinz,size=length(indinz)-1))
}
#---------------------------------------------------------------------------------------------------
swap=function(indinz){
  tmp=1:nloc
  indout=tmp[!tmp%in%indinz]
  
  keep=sample(indinz,size=length(indinz)-1)  
  include=sample(indout,size=1)
  sort(c(keep,include))
}
#---------------------------------------------------------------------------------------------------
birth=function(indinz){
  tmp=1:nloc
  indout=tmp[!tmp%in%indinz]
  include=sample(indout,size=1)
  sort(c(indinz,include))
}