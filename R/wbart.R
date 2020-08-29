
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

wbart=function(
x.train, y.train, x.test=matrix(0.0,0,0),
sparse=FALSE, theta=0, omega=1,
a=0.5, b=1, augment=FALSE, rho=NULL,
xinfo=matrix(0.0,0,0), usequants=FALSE,
cont=FALSE, rm.const=TRUE,
sigest=NA, sigdf=3, sigquant=.90,
k=2.0, power=2.0, base=.95,
sigmaf=NA, lambda=NA,
fmean=mean(y.train),
w=rep(1,length(y.train)),
ntree=200L, numcut=100L,
ndpost=1000L, nskip=100L, keepevery=1L,
nkeeptrain=ndpost, nkeeptest=ndpost,
nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
printevery=100L, transposed=FALSE,
# new arguments
location = NULL,
location_test = NULL,
range_select_sd=2.0, smoothness_select_sd=1.0,
sigma2_prior_a=10.0,
sigma2_prior_b=1.0,
tau2_prior_a=2.0,
tau2_prior_b=1.0,
range=4.6,
smoothness=0.0,
tau2=2.0,
range_prior_mean =5.0,
range_prior_sd = 1.5,
smoothness_prior_mean=0.0,
smoothness_prior_sd=0.5,
coordinates='ll'
)
{
#--------------------------------------------------
# locations in training set
names(location)[1:2] = c('x', 'y')
if (is.null(location_test)){} else {names(location_test)[1:2] = c('x', 'y')}
unique.location = as.data.frame(unique(location))
names(unique.location) = c('x', 'y')
unique.location$id = 1:nrow(unique.location)
location = merge(location, unique.location, by = c("x","y"), sort=FALSE)
Z = fastDummies::dummy_columns(location$id)
z = as.matrix(Z[,2:ncol(Z)])
if (coordinates=='ground') distance = as.matrix (dist( unique.location[c("x",'y')] ))
if (coordinates=='ll') distance = fields::rdist.earth(x1=unique.location[c("x",'y')] ,miles = FALSE)/100
nn = nrow(unique.location)
location <- as.matrix(location)
z1 = z


# locations in test set
if (is.null(location_test)){
new_loc = FALSE
distance_all = distance
z2 = NULL
location_all = unique.location
} else {
if (nrow(unique(location)) == nrow(unique(rbind(location[,1:2],location_test)))){
  new_loc = FALSE
  location_all = unique.location
  distance_all = distance
  ttest.location = dplyr::left_join(location_test, unique.location, by = c("x","y"),sort=FALSE)
  ttest.location <- rbind(ttest.location, unique.location)
  ttest.Z = fastDummies::dummy_columns(ttest.location$id)
  ttest.Z = as.matrix(ttest.Z[,2:ncol(ttest.Z)])
  ttest.Z = ttest.Z[1:nrow(location_test),]
  z2 = ttest.Z
} else {
  new_loc = TRUE
  unique.location.all = rbind(dplyr::anti_join(unique(location_test[,1:2]),unique.location[,1:2],by=c('x','y'),sort=FALSE),unique.location[,1:2])
  location_all = unique.location.all
  if (coordinates=='ground') distance_all = as.matrix (dist( unique.location.all ))
  if (coordinates=='ll') distance_all = fields::rdist.earth(x1=unique.location.all ,miles = FALSE)/100

  unique.location.all$id = 1:nrow(unique.location.all)
  ttest.location = left_join(location_test, unique.location.all, by = c("x","y"), sort=FALSE)
  ttest.location.2 <- rbind(ttest.location,data.frame(x=rep(0,nrow(unique.location.all)),
                                      y=rep(0,nrow(unique.location.all)),id=1:nrow(unique.location.all)))
  ttest.Z.2 = fastDummies::dummy_columns(ttest.location.2$id)
  ttest.Z.2 = as.matrix(ttest.Z.2[,2:ncol(ttest.Z.2)])
  z2 = ttest.Z.2[1:nrow(ttest.location),]
}
}
#Z2 = ttest.Z.2








#--------------------------------------------------
#data
n = length(y.train)

if(!transposed) {
    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if(length(x.test)>0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[ , temp$rm.const])
    }
    rm.const <- temp$rm.const
    grp <- temp$grp
    rm(temp)
}
else {
    rm.const <- NULL
    grp <- NULL
}

if(n!=ncol(x.train))
    stop('The length of y.train and the number of rows in x.train must be identical')

p = nrow(x.train)
np = ncol(x.test)
if(length(rho)==0) rho=p
if(length(rm.const)==0) rm.const <- 1:p
if(length(grp)==0) grp <- 1:p

##if(p>1 & length(numcut)==1) numcut=rep(numcut, p)

y.train = y.train-fmean
#--------------------------------------------------
#set nkeeps for thinning
if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
   nkeeptrain=ndpost
   cat('*****nkeeptrain set to ndpost\n')
}
if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
   nkeeptest=ndpost
   cat('*****nkeeptest set to ndpost\n')
}
if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
   nkeeptestmean=ndpost
   cat('*****nkeeptestmean set to ndpost\n')
}
if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
   nkeeptreedraws=ndpost
   cat('*****nkeeptreedraws set to ndpost\n')
}
#--------------------------------------------------
#prior
nu=sigdf
if(is.na(lambda)) {
   if(is.na(sigest)) {
      if(p < n) {
         df = data.frame(t(x.train),y.train)
         lmf = lm(y.train~.,df)
         sigest = summary(lmf)$sigma
      } else {
         sigest = sd(y.train)
      }
   }
   qchi = qchisq(1.0-sigquant,nu)
   lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
} else {
   sigest=sqrt(lambda)
}

if(is.na(sigmaf)) {
   tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
} else {
   tau = sigmaf/sqrt(ntree)
}
#--------------------------------------------------
ptm <- proc.time()
#call
res = .Call("cwbart",PACKAGE="BART",
            n,  #number of observations in training data
            p,  #dimension of x
            np, #number of observations in test data
            x.train,   #pxn training data x
            y.train,   #pxn training data x
            x.test,   #p*np test data x
            ntree,
            numcut,
            ndpost*keepevery,
            nskip,
            power,
            base,
            tau,
            nu,
            lambda,
            sigest,
            w,
            sparse,
            theta,
            omega,
            grp,
            a,
            b,
            rho,
            augment,
            nkeeptrain,
            nkeeptest,
            nkeeptestmean,
            nkeeptreedraws,
            printevery,
            
            z,
            distance,
            distance_all,
            new_loc,
            
            xinfo,
            
            range_select_sd,
            smoothness_select_sd,
            sigma2_prior_a,
            sigma2_prior_b,
            tau2_prior_a,
            tau2_prior_b,
            range,
            smoothness,
            tau2,
            range_prior_mean,
            range_prior_sd,
            smoothness_prior_mean,
            smoothness_prior_sd
)

res$proc.time <- proc.time()-ptm

res$z1 = z1
res$z2 = z2
res$locations_train = unique.location
res$locations_test = dplyr::anti_join(unique(location_test[,1:2]),unique.location[,1:2],by=c('x','y'),sort=FALSE)
#res$location_all = ttest.location.2

res$mu = fmean
res$yhat.train.mean = res$yhat.train.mean+fmean
res$yhat.train = res$yhat.train+fmean
res$yhat.test.mean = res$yhat.test.mean+fmean
res$yhat.test = res$yhat.test+fmean
if(nkeeptreedraws>0)
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
##res$nkeeptreedraws=nkeeptreedraws
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$rm.const <- rm.const
attr(res, 'class') <- 'wbart'
return(res)
}
