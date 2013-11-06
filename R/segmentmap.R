#
# functions for retention time aligment of swath data.
#


## function for creating the bins in mz and rt
## returns a list
.prepareBins <- function(con,mzbinw=1,rtbinw=3){
    swathinfo <- mapSummary(con)
    swathinfo <- swathinfo[swathinfo$mslevel==2,]
    mzwindow <- swathinfo[c("minMZsw","maxMZsw")]
    mzrange <- c(floor(min(swathinfo$minMZ)*10)/10-1,ceiling(max(swathinfo$maxMZ)*10)/10+1)
    rtrange <- c(floor(min(swathinfo$minRT)*10)/10-1,ceiling(max(swathinfo$maxRT)*10)/10+1)

    breaksmz <- seq(mzrange[1],mzrange[2],mzbinw)
    nBmz <- length(breaksmz)
    midsmz <- 0.5 * (breaksmz[-1L] + breaksmz[-nBmz])

    breaksrt <- seq(rtrange[1],rtrange[2],rtbinw)
    nBrt <- length(breaksrt)
    midsrt <- 0.5 * (breaksrt[-1L] + breaksrt[-nBrt])

    return(list( breaksrt=breaksrt , midsrt=midsrt  , breaksmz=breaksmz ,  midsmz=midsmz ))
}

##
##
.createMapImage <- function( con , idswath , bandm, rtshift=FALSE )
{
    swathinfo <- mapSummary(con)
    swathinfo <- swathinfo[swathinfo$id==idswath,]
    mzwindow <- swathinfo[c("minMZsw","maxMZsw")]
   
    data<-getMZRTVol(con, idswath)
    ##remove swath window
    datrm <- data[!(data$MZ > mzwindow$minMZsw & data$MZ < mzwindow$maxMZsw),]
    
    if(rtshift){
        shiftV <- (bandm$midsrt[1]- bandm$breaksrt[1])
    }else{
        shiftV <- 0
    }
    print(shiftV)
    categrt <- cut(datrm$RT,breaks=bandm$breaksrt+shiftV)
    categmz <- cut(datrm$MZ,breaks=bandm$breaksmz)
    tmp <- tapply(datrm$Volume,list(categrt,categmz),sum)
    rownames(tmp) <- bandm$midsrt+shiftV
    colnames(tmp) <- bandm$midsmz
    return(tmp)
}

## colnames are the masses, rownames are rts
.prepare2Compare <- function(img){
    img[is.na(img)] <- 0
    imgtr <- apply(img,1,function(x,y){x*y^2},as.numeric(colnames(img)))
    imgtr <- t(imgtr)
    return(imgtr)
}

##ensure that vector has more than minelems
.dropShortVecs <- function(imgtr,MARGIN=1,minelems=10){
   imgtr <- apply(imgtr,MARGIN,function( x )
                  {
                      if(sum(x>0) > minelems){
                          return(x)
                      }else{
                          return(rep(0,length(x)))
                      }
                  }
       )
    return(t(imgtr))
}


.vnormRows <- function(imgtr,MARGIN=1){
    imgtr <- apply(imgtr,MARGIN,function(x){ v <- sqrt( sum(x*x) ) ;if(v>0){x <- x / v};return(x) })
    return(t(imgtr))
}

createMapImage <- function(con,idswaths,mzbinw=0.9,rtbinw=3,
                           minelems=30,rtshift=FALSE){
    bandm <- .prepareBins(con,mzbinw,rtbinw)
    img <- NULL
    for(i in idswaths){
        tmp <- .createMapImage(con,i,bandm,rtshift)
        tmp <- .prepare2Compare(tmp)
        img <- cbind(img,tmp)
    }
    img <- .dropShortVecs(img,minelems=minelems)
    img <- .vnormRows(img)
    return(img)
}


## img1, img2 , threshold of dotprod
getMatches <- function(img1,img2,threshold=0.5){
    library(Matrix)
    iMgtr1 <- Matrix(img1)
    iMgtr2 <- Matrix(img2)
    X <- (iMgtr1) %*% t(iMgtr2)
    idx <- which(X>threshold,arr.ind=T)
    V <- X[X>threshold]
    x <- as.numeric(rownames(img1))[idx[,1]]
    y <- as.numeric(rownames(img1))[idx[,2]]
    return(list( x=x , y=y , V=V , X=X))
}


# xyV = result of getMatches
computeAlignment <- function(xyV, rtdevmax = 180, k.runmed = 19,f.lowess=0.05){
    x <- NULL
    y <- NULL
    V <- NULL
    for(i in xyV)
        {
            print(length(i$x))
            x <- c(x,i$x)
            y <- c(y,i$y)
            V <- c(V,i$V)
        }
    print(length(x))
    diffrt <- x-y
    meanrt <- (x+y)/2

    select <- diffrt < rtdevmax & diffrt >- rtdevmax
    
    diffrt <- diffrt[select] 
    meanrt <- meanrt[select]
    V <- V[select]

    rtorder <- order(meanrt)
    meanrt <- meanrt[rtorder]
    diffrt <- diffrt[rtorder]
    V <- V[rtorder]
    print(length(diffrt))
    diffrtmedian <- runmed(diffrt,k.runmed)
    return(list( meanrt=meanrt,
                diffrt=diffrt,
                diffrtmedian=diffrtmedian,
                V=V,
                lowess=lowess(meanrt,diffrtmedian,V,f=f.lowess))
           )
    
}

visAlignmentMatrix <- function(xyV,range=NULL )
{
    smoothGauss <- function(x,width=5){g <- dnorm(-20:20,0,width); x <- filter(x,g,circular=TRUE);return(x)}
    print(length(range))
   
    if(length(range)==2){
        print("test")
        XF <- xyV$X[range[1]:range[2],range[1]:range[2]]
    }
    else{
        XF <- xyV$X
    }
    XF <- apply(XF,1,smoothGauss)
    XF <- apply(XF,1,smoothGauss)
    image(XF,main="smoothed dotproduct matrix")
    abline(c(0,1),col="gray")
}
