##author : Witold Wolski wewolski@gmail.com

## functions for transition group finding.

### Peak group finding scoring and quantification:
## Extract features given transition group
.getFeatures <- function( precursor , con , mzerror=0.1, rterror=NULL )
{
    parmass=precursor$mz
    parrt=precursor$rt
    rtrange=NULL
    if(length(rterror)>0)
        {
            rtrange = c(parrt-rterror,parrt+rterror)
        }
    
    ms1mapId <- getMS1id(con)
    massrange=c(parmass-mzerror,parmass+mzerror)
    parents <- getFeaturesFlex(con,ms1mapId,massrange=massrange,rtrange=rtrange,features="ms1view")

    swathid = getSwath( con , parmass )
    prodmass <- precursor$products[,"mz"]
    if(length(swathid)==0)
        {
            stop("no swath id given")
        }
    products <- NULL
    for(j in 1:length(prodmass) )
        {
            products <- rbind(products,
                              getFeaturesFlex(con, swathid ,
                                              c(prodmass[j]-mzerror , prodmass[j]+mzerror ),
                                              rtrange=rtrange ) )
        }
    return( list(parents= parents ,products= products) )
}

## 
## group features according RT
##
.groupFeatures<-function(product,parent,rterror=1.5,minsize=2)
{
    ## order products by mass
    product = product[order(product$Volume,decreasing=T),]
    # created feature groups
    groups <- list()
    count = 1
    while(dim(product)[1]>1)
        {
            largest<-product[1,]
            ## find products and parents given RT
            ## this will find also parent itself
            w <- which( (product$RT < (largest$RT + rterror)) & (product$RT > (largest$RT-rterror)))
            wpar <- which( (parent$RT < (largest$RT + rterror)) & (parent$RT > (largest$RT-rterror)))

            parf<-NULL
            group<-NULL
            if( length(wpar) > 0 ){
                parf <- cbind( rep(1,length(wpar)), parent[wpar,] )
                parent <- parent[-wpar,]#remove parent
            }
            if( length(w) > 0){
                group <- cbind( rep(2,length(w)), product[w,])
                product <- product[-w,]#remove group                
            }
            if(dim(group)[1] >= minsize){
                groups[[count]] <- list( parent = parf , prod = group )
                count = count+1                
            }
        }
    return(groups)
}

##
## 
##
.matchPeaks <- function(xmass,ymass,error=0.05){
    res <- NULL
    for(i in 1:length(xmass)){
        for(j in 1:length(ymass)){
                if(abs(ymass[j] - xmass[i])<error){
                    res <- rbind(res,c(i,j))
                }
            }
    }
    return(res)
}

### do scoring ###
scoreGroup <- function(trans,group){
    parent <- group$parent
    product <- group$prod
    RTmean <- mean(product$RT)
    RTmaxlocmean <- mean(product$maximumLocationRT)
    prodVol <- sum(product$Volume)
    prodInt <- sum(product$Max)
    RTcmass <- mean(product$centerOfMassRT)

    nrtexpected <- dim(trans$product)[1]
    nrtfound <- dim(product)[1]

    ## compute mean average mass deviation
    tmp <- abs(diff((sort(c(product$MZ,trans$product[,"mz"])))))
    idx <- which(tmp<0.1)
    massdev <- mean(tmp[idx])

    ## compute rank correlation
    idx <- .matchPeaks(product$MZ, trans$product[,"mz"])
    spearman <- cor( product$Volume[ idx[,1] ] ,  trans$product[,"intensity"][ idx[,2] ],method="spearman" )

    
    res <- c(transname=trans$name,
             sequence=trans$sequence,
             rt=trans$rt,
             pmass=trans$mz,
             decoy=trans$decoy,
             RTmean=RTmean,
             RTmaxlocmean=RTmaxlocmean,
             prodVol=prodVol,
             prodInt=prodInt,
             RTcmass=RTcmass,
             nrtexpected=nrtexpected,
             nrtfound=nrtfound,
             massdev=massdev,
             spearman = spearman
             )
    return(res)
}

## irt calibration ##
findTransitionsForIRT <- function(irttrans2
                                  ,con
                                  ,mzerror=0.03# for feature extraction
                                  ,rterror=2.5## for peak grouping
                                  ,relax=0
                          ){
    irts <- NULL
    for(j in 1:length(irttrans2))
        {
            trans<-irttrans2[[j]]
            ## retrieve features matching transition mass
            features = .getFeatures( trans , con , mzerror=mzerror)
            ## group features and select best group
            xx = .groupFeatures(features$products,features$parents,rterror=rterror,minsize=dim(trans$products)[1]-relax)
            cat("nr groups", length(xx), "\n")
            if(length(xx)>0)
                {
                                        # get retention time from first group.
                    time <- mean(c(xx[[1]]$parent[,"RT"],xx[[1]]$prod[,"RT"]))
                    irts <- rbind(irts,c(j,time)) # get experimental retention time
                }
        }

    reft <- unlist( lapply( irttrans2 , function(x){ return(x$rt) } ) )
    return(list(refrt=reft,irts=irts))
}



## find reference irTS
calibrateRTTransitions <- function(con, irttrans2, alltrans ,
                                   mzerror=0.03 ,
                                   rterror=2.5
                                   ,relax=0){
    library(MASS)
    tmp <-  findTransitionsForIRT(irttrans2,con,mzerror=mzerror,rterror=rterror
                                  ,relax=relax)
    ## the modelling part
    reft <- tmp$refrt[tmp$irts[,1]]
    irts <- tmp$irts[,2]
    if(length(reft)==length(irts)){
        plot(reft, irts)
        model <- rlm(irts~reft)
        abline(model)
        newtrans<-mapRT(alltrans,model)
        return(newtrans)
    }
    return(NULL)
}

getIRTcalibrationmodel <- function (con, 
                                    irttrans2, 
                                    mzerror = 0.03,
                                    rterror = 2.5,
                                    relax = 0
) 
{
  library(MASS)
  tmp <- findTransitionsForIRT(irttrans2, con, mzerror = mzerror, 
                               rterror = rterror, relax = relax)
  reft <- tmp$refrt[tmp$irts[, 1]]
  irts <- tmp$irts[, 2]
  if (length(reft) == length(irts)) {
    plot(reft, irts)
    model <- rlm(irts ~ reft)
    return(model)
  }
  return(NULL)
}



##
findTransitions <-  function(newtrans,con
                             ,mzerror=0.03# for feature extraction
                             ,extRTerr=90
                             ,rterror=2.5## for peak grouping
                             ,relax=0 ## 
                             ){
    res <- vector(length=length(newtrans),mode="list")
    
    score <- matrix("",nrow=length(newtrans),ncol=14)
    groups <- vector(length(newtrans),mode="list")
    count <- 1
    for(j in 1:length(newtrans))
        {
            trans<-newtrans[[j]]
            ##trans
            ## retrieve features matching transition mass
            features = .getFeatures( trans , con , mzerror=0.03, rterror=extRTerr )
            ## group features and select best group
            featuregroups=.groupFeatures(features$products,features$parents,rterror=rterror,minsize=dim(trans$products)[1]-relax)
            if(length(featuregroups)>0){
                for(i in 1:length(featuregroups)){
                    dd <- scoreGroup(trans,featuregroups[[i]])
                    score[count,] <- dd
                    groups[[count]] <- list(trans=trans, featuregroups= featuregroups[[i]])
                    count <- count + 1
                }
            }
        }
    print(count)

    if(count > 1){
        colnames(score) <- names(dd)
        groups <- groups[1:(count-1)]
        score <- score[1:(count-1),,drop=F]
        score <- as.data.frame(score,stringsAsFactors=FALSE)
        print(names(score))
        for(i in 3:dim(score)[2])
            {
                score[[i]] <- as.numeric(score[[i]])
            }
        return( list(score=score, groups=groups) )
    }
    return(NULL)
}#end findFeatures


#
#
visTransitions <- function(con,groups,main="",sleep=1){
    if(length(groups)>0){
        for(i in 1:length(groups)){
            plotGroup(con,groups[[i]]$trans,groups[[i]]$featuregroups,main=main)
            Sys.sleep(sleep)
        }
    }
}

##
##export help
##

toSimpleFormat <- function(sc){
    len <- dim(sc)[1]
    tmp <- data.frame(dataID=rep(dataset,len),
                      dataDilution=rep(dil,len),
                      pepseq=sc$sequence,
                      marea=sc$prodVol,
                      mcarea=sc$prodVol,
                      mrt=sc$rt,
                      mback=rep(0,len),stringsAsFactors=FALSE
                      )
    return(tmp)
}
