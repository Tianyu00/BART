
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2019 Rodney Sparapani

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


mc.gbmm <- function(
                     x.train, y.train,
                     x.test=matrix(0,0,0), type='wbart',
                     u.train=NULL, B=NULL,
                     ntype=as.integer(
                         factor(type,
                                levels=c('wbart', 'pbart'))),
                                ##levels=c('wbart', 'pbart', 'lbart'))),
                     sparse=FALSE, theta=0, omega=1,
                     a=0.5, b=1, augment=FALSE, rho=NULL,
                     xinfo=matrix(0,0,0), usequants=FALSE,
                     rm.const=TRUE,
                     sigest=NA, sigdf=3, sigquant=0.90,
                     k=2, power=2, base=0.95,
                     ##sigmaf=NA,
                     lambda=NA, tau.num=c(NA, 3, 6)[ntype],
                     ##tau.interval=0.9973,
                     offset=NULL, ##w=rep(1, length(y.train)),
                     ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
                     ndpost=1000L, nskip=100L,
                     keepevery=c(1L, 10L, 10L)[ntype],
                     printevery=100L, transposed=FALSE,
                     hostname=FALSE,
                     mc.cores = 2L, nice = 19L, seed = 99L
                     )
{
    if(is.na(ntype))
        stop("type argument must be set to either 'wbart' or 'pbart'")
        ##stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    if(length(u.train)==0)
        stop("the random effects indices must be provided")

    check <- unique(sort(y.train))

    if(length(check)==2) {
        if(!all(check==0:1))
            stop('Binary y.train must be coded as 0 and 1')
        if(type=='wbart')
            stop("The outcome is binary so set type to 'pbart'")
            ##stop("The outcome is binary so set type to 'pbart' or 'lbart'")
    }

    if(.Platform$OS.type!='unix')
        stop('parallel::mcparallel/mccollect do not exist on windows')

    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    if(!transposed) {
        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        ## if(length(x.test)>0)
        ##     x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        if(length(x.test)>0) {
            x.test = bartModelMatrix(x.test)
            x.test = t(x.test[ , temp$rm.const])
        }
        rm.const <- temp$rm.const
        rm(temp)
    }

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected

    mc.ndpost <- ceiling(ndpost/mc.cores)

    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
            gbmm(x.train=x.train, y.train=y.train,
                  x.test=x.test, type=type, ntype=ntype,
                  u.train=u.train, B=B,
                  sparse=sparse, theta=theta, omega=omega,
                  a=a, b=b, augment=augment, rho=rho,
                  xinfo=xinfo, usequants=usequants,
                  rm.const=rm.const,
                  sigest=sigest, sigdf=sigdf, sigquant=sigquant,
                  k=k, power=power, base=base,
                  ##sigmaf=sigmaf,
                  lambda=lambda, tau.num=tau.num,
                  ##tau.interval=tau.interval,
                  offset=offset, ##w=w,
                  ntree=ntree, numcut=numcut,
                  ndpost=mc.ndpost, nskip=nskip,
                  keepevery=keepevery, printevery=printevery,
                  transposed=TRUE, hostname=hostname)},
            silent=(i!=1))
        ## to avoid duplication of output
        ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]

    if(mc.cores==1 | attr(post, 'class')!=type) return(post)
    else {
        if(class(rm.const)[1]!='logical') post$rm.const <- rm.const

        post$ndpost <- mc.cores*mc.ndpost

        p <- nrow(x.train[post$rm.const, ])

        old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree),
                           ' ', as.character(p))
        old.stop <- nchar(old.text)

        post$treedraws$trees <- sub(old.text,
                                    paste0(as.character(post$ndpost), ' ',
                                           as.character(ntree), ' ',
                                           as.character(p)),
                                    post$treedraws$trees)

        keeptest <- length(x.test)>0

        for(i in 2:mc.cores) {
            post$hostname[i] <- post.list[[i]]$hostname
            
            post$yhat.train <- rbind(post$yhat.train,
                                     post.list[[i]]$yhat.train)

            if(keeptest) post$yhat.test <- rbind(post$yhat.test,
                                                 post.list[[i]]$yhat.test)

            if(type=='wbart')
                post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)

            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
            post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)

            post$sd.u <- cbind(post$sd.u, post.list[[i]]$sd.u)

            post$treedraws$trees <- paste0(post$treedraws$trees,
                                           substr(post.list[[i]]$treedraws$trees, old.stop+2,
                                                  nchar(post.list[[i]]$treedraws$trees)))

            post$proc.time['elapsed'] <- max(post$proc.time['elapsed'],
                                             post.list[[i]]$proc.time['elapsed'])
            for(j in 1:5)
                if(j!=3)
                    post$proc.time[j] <- post$proc.time[j]+post.list[[i]]$proc.time[j]
        }

        post$an.train <- apply(post$yhat.train, 1, function(x) !any(is.nan(x)))

        if(type=='wbart') {
            post$yhat.train.mean <-
                apply(post$yhat.train[post$an.train, ], 2, mean)

            if(keeptest)
                post$yhat.test.mean <-
                    apply(post$yhat.test[post$an.train, ], 2, mean)
        } else {
            post$prob.train.mean <-
                apply(post$prob.train[post$an.train, ], 2, mean)

            if(keeptest)
                post$prob.test.mean <-
                    apply(post$prob.test[post$an.train, ], 2, mean)
        }

        post$varcount.mean <- apply(post$varcount, 2, mean)
        post$varprob.mean <- apply(post$varprob, 2, mean)

        if(all(post$an.train)) post$an.train <- NULL ## no NAs
        
        attr(post, 'class') <- type

        return(post)
    }
}
