#' creates a chromatogram for each feature
#' value - matrix with nr row equal to dim products
#' @export 
feature2Chrom = function(products,rtaxis){
  res <- NULL    
  if(dim(products)[1]==0){
    return(res)
  }
  for(i in 1:dim(products)[1]){
    y <- rep(0,length(rtaxis))
    product<-products[i,]    
    rtProjection=readBin(product$rtProjection[[1]],"numeric",size=4,n=100)
    minRTIndex<-product$minRTIndex
    y[minRTIndex:(minRTIndex+length(rtProjection)-1)] <- rtProjection
    res <- rbind(res,y)
  }
  return(res)
}
#' creates one chromatogram for all features (i.e. extracted at the same mass)
#' @export
features2Chrom = function( products , rtaxis ){
  res <- NULL
  if(dim(products)[1]==0){
    return(res)
  }
  y <- rep(0,length(rtaxis))
  for(i in 1:dim(products)[1]){
    product<-products[i,]    
    rtProjection=readBin(product$rtProjection[[1]],"numeric",size=4,n=100)
    minRTIndex<-product$minRTIndex
    w <- rep(0,length(rtaxis))
    w[minRTIndex:(minRTIndex+length(rtProjection)-1)] <- rtProjection
    y = y+w #add it in case of overlapping features.
  }
  return(y)
}
#' plot Products
#' @export
plotProducts=function(rtaxis,chroms,xlim=NULL){
  if(length(xlim)==0)
    xlim=range(rtaxis)
  ylim = c(0,(max(chroms)))
  plot(rtaxis,rep(0,length(rtaxis)),type="n",xlim=xlim,ylim=ylim)
  for(i in 1:dim(chroms)[1]){
    lines(rtaxis,(chroms[i,]),col=i)
  }
  
}

.plotBottom=function( rtaxis , chroms, xlim=NULL ){
  ##plot on bottom
  if(length(chroms)==0)
    return(NULL)
  if(length(xlim)==0)
    xlim=range(rtaxis)
  
  if(!is.matrix(chroms)){
    chroms <- matrix(chroms,nrow=1)
  }
  #set the range
  ylim = c(-(max(chroms)),(max(chroms)))
  
  plot(rtaxis,rep(0,length(rtaxis)) , type="n" , ylim=ylim , xlim=xlim , axes=FALSE, ylab="",xlab="")
  bbb=axTicks(4)
  axis(4,at=bbb,labels=-bbb)
  for(i in 1:dim(chroms)[1]){
    lines(rtaxis,-(chroms[i,]),col=i)
  }
}

.plotTop=function(rtaxis, chroms,main="", xlim=NULL){
  if(length(chroms)==0){
    return(NULL)
  }
  if(length(xlim)==0)
    xlim=range(rtaxis)
  if(!is.matrix(chroms)){
    chroms <- matrix(chroms,nrow=1)
  }
  ylim = c(-(max(chroms)),(max(chroms)))
  plot(rtaxis,rep(0,length(rtaxis)),type="n",main=main,ylim=ylim,xlim=xlim,ylab="MS2 ion intensities",xlab="rt")
  for(i in 1:dim(chroms)[1]){
    lines(rtaxis, (chroms[i,]),col=i)
  }
}
.plotFeatures = function(
  rtaxis, chroms,#plot on top
  rtaxis2, chroms2# plot on bottom
  , main="",xlim=NULL
)
{
  .plotTop(rtaxis,chroms,main=main,xlim)
  par(new=TRUE)
  .plotBottom(rtaxis2,chroms2,xlim)
  par(new=FALSE)
}
#' plots precursor
#' @param precursor - precursor mass
#' @export
plotprecursor=function(precursor,con,error=0.1,xlim=NULL){
  parmass<-precursor$mz
  ms1mapId <- getMS1id(con)
  
  parents <- getFeaturesFlex(con,ms1mapId,massrange=c(parmass-error,parmass+error))
  rtaxis1<-getRTAxis(con,ms1mapId)
  precurs <- features2Chrom(parents,rtaxis1)
  
  length(precurs)
  length(rtaxis1)
  
  swathid = getSwath(con,parmass)
  rtaxis<-getRTAxis(con,swathid)
  prodmass <- precursor$products[,"mz"]
  #print("test")
  xx<-NULL
  for(j in 1:length(prodmass)){
    products <- getFeaturesFlex(con, swathid ,c(prodmass[j]-error,prodmass[j]+error))
    xx <- rbind(xx,features2Chrom(products,rtaxis))
    print(dim(products))
  }
  .plotFeatures(rtaxis,xx,rtaxis1,precurs,xlim)
  return(invisible(list(rt1=rtaxis,chroms1=xx,rt2=rtaxis1,chroms2=precurs)))
}
#' plots a group
#' @param plots an transition group
#' @export
plotGroup = function(
  con,transition,
  group,main=""
  ,xlim=NULL
  ){
  parent <- group$parent
  product <- group$prod
  
  minIdx <- 10000
  maxIdx <- 0
  chromsParent <-NULL
  if(length(parent)[1]>0){
    ms1mapId <- getMS1id(con)
    rtaxis1 <-getRTAxis(con,ms1mapId)
    chromsParent <- features2Chrom(parent,rtaxis1)
    minIdx <- parent$minRTIndex
    maxIdx <- minIdx+parent$rtExtend
  }
  swathid = getSwath(con,transition$mz)

  rtaxis <- getRTAxis(con,swathid)
  xx <- NULL
  for(i in 1:dim(product)[1])
  {
    chromsProd <- features2Chrom( product[i,] , rtaxis )
    minIdx <- min(c(minIdx,product[i,]$minRTIndex))
    maxIdx <- max(c(maxIdx,product[i,]$minRTIndex+product[i,]$rtExtend))
    xx<-rbind(xx,chromsProd)
  }
  
  plotFeatures(rtaxis[minIdx:maxIdx], xx[,minIdx:maxIdx], rtaxis1[minIdx:maxIdx], chromsParent[minIdx:maxIdx],main=main,xlim)
  legend("bottomright",c(transition$sequence,transition$name))
}
#' peak group detection visualization:
#' @export
histOfMassFG = function(tmp){
  par(mfrow=c(2,1))
  parent = tmp$parents$RT
  product = tmp$products$RT
  ##use one RT axis
  br=seq( min(c(parent,product))-1 , max(c(parent,product))+1 , by=2)
  br2 = br+1
  pax1 = hist(parent , breaks=br)$counts
  pax2 = hist(parent , breaks=br2, add=T )$counts
  
  prx1=hist(product,breaks=br)$counts
  br2 = br+0.5
  prx2=hist(product,breaks=br2,add=T)$counts
  return(list(prod1=prx1,prod2=prx2,par1=pax1,par2=pax2,br1=br,br2=br2))
}
