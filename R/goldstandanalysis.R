.computePairRatio <- function(tr1x,tr2x){
    nam <- c(tr1x$pepseq,tr1x$dataID , tr2x$dataID)
    dil <- c(tr1x$dataDilution, tr2x$dataDilution )
    rarea <- c(tr1x$marea,tr2x$marea)
    rcarea <- c(tr1x$mcarea,tr2x$mcarea)
    mrt <- c(tr1x$mrt,tr2x$mrt)
    mback <- c(tr1x$mback,tr2x$mback)
    tmp <- c(dil, rarea, rcarea, mrt, mback)
    return(list(nam=nam,nr=tmp))
}# end compute pairs
## dilution of tr1 > tr2
.allpossibleratiosFor2Dilutions <- function(tr1,tr2){
    nam <- NULL
    nr <- NULL
    for(i in 1:dim(tr1)[1])
        {
            for(j in 1:dim(tr2)[1])
                {
                    tr1x <- tr1[i,]
                    tr2x <- tr2[j,]
                    if(tr1x$dataID != tr2x$dataID){#don't compute for same dataset
                        x <- .computePairRatio(tr1x,tr2x)
                        nam <- rbind( nam , x$nam )
                        nr <- rbind( nr, x$nr )
                    }
                }
        }
    return (list(nam=nam, nr=nr))
}
.ratiosforTransition=function(trw){
    x<-sort(unique(trw$dataDilution),decreasing=TRUE)
    nam <- NULL
    nr <- NULL
    for(i in x)
        {
            for(j in x){
                if(j <= i){
                    tr1 <- trw[trw$dataDilution==i,]
                    tr2 <- trw[trw$dataDilution==j,]
                    res <- .allpossibleratiosFor2Dilutions(tr1,tr2)
                    nam <- rbind(nam, res$nam)
                    nr <- rbind(nr, res$nr)
                }
            }
        }
    return(list(nam=nam,nr=nr))
}
#' extract features from database
#' @param notna matrix with target data
#' @param con connection
#' @export
#' 
extractfeatures = function(notna, con , mzerror = 0.03, rterror=10)
{
  allfeat <- NULL
  for(i in 1:dim(notna)[1]){
    target <- notna[i,]
    massrange=c(target$ProductMz-mzerror,target$ProductMz+mzerror)
    rtrange=c(target$iRTs-rterror,target$iRTs+rterror)
    mapSummary(con)
    swathid = getSwath( con , target$PrecursorMz )
    target$PrecursorMz
    swathid
    feat <- getMZRTVol(con,swathid,mzrange=massrange,rtrange=rtrange)
    if(dim(feat)[1] > 0){
      for(j in dim(feat)[1]){
        featr <- c(target,feat[j,])
        allfeat<-rbind(allfeat,featr)
      }
    }else{
      featr <- c(target,rep(NA,length(feat)))
      allfeat <- rbind(allfeat,featr)
    }
  }
  return(allfeat)
}
#' ratios of transitions
#' @export
ratiosOfTransitionList=function(tr){
    pepseq <- (unique(tr$pepseq))
    res <- vector(length(pepseq),mode="list")
    length(res)
    for(ipep in 1:length(pepseq)){
        pep <- pepseq[ipep]
        cat("i",ipep,"pep",pep,"\n")
        trw <- tr[tr$pepseq == pep,]
        res[[ipep]] <- .ratiosforTransition(trw)
    }
    
    nam <- NULL
    nr <- NULL
    for(i in 1:length(res))
        {
            nam <- rbind(nam,res[[i]]$nam)
            nr <- rbind(nr,res[[i]]$nr)
        }

    res <- data.frame(nam,nr,stringsAsFactors=FALSE)
    colnames(res) <- c("pepseq","dataID1","dataID2","dil1","dil2",
                       "area1","area2","carea1","carea2","rt1","rt2","back1","mback2")
    return(res)
}
#' do data preparation
#' @export
prepareData <- function(tmp){
    replicA <- tmp$ReplicateName
    
    replic <- strsplit(as.character(replicA),"_")
    namDil <- c(unlist(replic))
    dsID <- namDil[ seq(2, length(namDil), 4) ]
    dsDil <- as.numeric(namDil[ seq(3, length(namDil), 4) ])

    
    iddil <- data.frame(dataID=dsID,dataDilution=dsDil,stringsAsFactors=FALSE)
    golds <- data.frame(iddil,tmp,stringsAsFactors=FALSE)
    #golds <- golds[,-5]
    
    indic <- (golds$RetentionTime=="#N/A")
    golds[indic,c(10:13)] <- as.numeric(NA)
    golds$RetentionTime <- as.numeric(golds$RetentionTime)
    golds$Area <- as.numeric(golds$Area)
    golds$Background <- as.numeric(golds$Background)
    golds$PeakRank <- as.numeric(golds$PeakRank)

    return(golds)
}
###############################
###############################
.getNewID <- function(x){
    pepse <- as.character(x[1])
    dataid <- as.character(as.numeric(x[2]))
    datadilution <- as.character(x[3])
    newid <- paste(pepse,dataid,datadilution,sep="_")
    return(newid)
}
.getKeys <- function(x){
    pepse <- as.character(x[1])
    dataid <- as.character(as.numeric(x[2]))
    datadilution <- as.character(x[3])
    newid <- c(pepse,dataid,datadilution)
    return(newid)
}
### compute transition stats for gold standard
.computeTransStat <- function(tran){
    tran <- data.frame(tran,stringsAsFactors=FALSE)
    area <- sum(tran$Area,na.rm=TRUE)
    areac <- sum(tran$Area-tran$Background,na.rm=TRUE)
    rt <- mean(tran$RetentionTime,na.rm=TRUE)
    back <- sum(tran$Background,na.rm=TRUE)
    return(c( area=area , carea=areac , rt=rt , back=back ))
}
#' compute summaries for parents from transitions
#'
#' @export
computeTransitionStats <- function(golds)
{
    newid <- apply(golds,1,.getNewID)
    unewid <- unique(newid)
    transres <- matrix(0,length(unewid),7)
    for(id in 1:length(unewid)){
        key <- unewid[id]
        indic <- which(key==newid)
        tran <- golds[ indic , ]
        res <- c(unlist(strsplit(key,"_")),.computeTransStat(tran))
        transres[id,]<-res
    }
    transresd <- data.frame("dataID"=as.character(transres[,1])
                        ,dataDilution=as.numeric(transres[,2])
                        ,pepseq=as.character(transres[,3])
                        ,marea=as.numeric(transres[,4])
                        ,mcarea=as.numeric(transres[,5])
                        ,mrt = as.numeric(transres[,6])
                        ,mback=as.numeric(transres[,7])
                        ,stringsAsFactors=FALSE)
    return(transresd)
}

    
