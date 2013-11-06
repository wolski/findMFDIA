###### Extract peaks from ms1 and attempt deisotoping
# Sort all ms1 features according to volume
# start with most intense feature
# select all MS1 features in the vicinity of that feature from the ms1 map:
# within RT +-rterror sec and MZ+2.2 mz.
# at the moment all these features are considered being istopes.
# use the masses of the features to determine the charge state
# remove the feature and the isotope from the list of unprocessed features
# add the feature with determined charge to the resultset.
######

##
## selects RT MZ and Volume from database
##

.getPseudoSpec = function(
  con,
  peak, # vector c(RT, MZ)
  errorRT = 1.5, # error of RT in [s]
  features = "bestfeatures"
)
{
  idswath <-  getSwath(con,peak["MZ"])
  rtmin <- peak["RT"]-errorRT
  rtmax <- peak["RT"]+errorRT
  coord <- getMZRTVol(con,idswath,rtrange=c(rtmin,rtmax),features=features)
  return(coord)
}

## write mgf files
.writeMGF = function(fcon,filename,parpeak,peaks,specid=1){
  parpeak<-parpeak
  cat("BEGIN IONS\n", file=fcon)
  cat("TITLE=",filename,".",specid,"\n",file=fcon,sep="")
  cat("PEPMASS=",parpeak["MZ"]," ",parpeak["Volume"],"\n",file=fcon,sep="")
  cat("RTINSECONDS=",parpeak["RT"],"\n",file=fcon,sep="")
  if(!is.nan(parpeak["Charge"])){
    #cat("CHARGE=",parpeak["Charge"],"+\n",file=con,sep="")
  }
  for(i in 1:dim(peaks)[1])
  {
    cat(peaks[i,"MZ"],peaks[i,"Volume"],"\n",file=fcon,sep=" ")
  }
  cat("END IONS\n\n",file=fcon,sep="")
}

writeSpecList = function(fcon, filename, parpeak, peaklists )
{
  for(i in 1:length(peaklists))
  {
    .writeMGF(fcon,filename,parpeak,peaklists[[i]],specid = i)
  }
}

## generate Pseudospectrum given ms1 peaks
.generatePseudospec=function(
  dbcon,
  filecon,
  filename,
  peaks,
  errorRT=1.5, # error of RT in [s]
  minlength = 5
)
{
  generSummary <-list(found=0,notfound=0,length=NULL)
  for(i in 1:dim(peaks)[1])
  {
    tmp = .getPseudoSpec(dbcon,peaks[i,],errorRT=errorRT)
    if(dim(tmp)[1]>minlength){
      .writeMGF(filecon,filename,unlist(peaks[i,]),tmp,i)
      generSummary$found = generSummary$found + 1
      generSummary$length <- c(generSummary$length, dim(tmp)[1])
    }
    else{
      #cat("no fragments found for ms1 :" , i, "\n")
      generSummary$notfound = generSummary$notfound + 1
    }
  }
  return(generSummary)
}

# facade function
# fname - file to read
# outdir - directory to write to - filename will be assempled.
# TODO separate MS1 selection function:

getMS1seeds = function( con , rtext=3 , mzext=3 , ms1vol = 500, rtrange=NULL )
{
  ### Confine search to best features
  createBestFeaturesView(con,rtext,mzext,ms1vol)
  ms1id <- getMS1id(con)
  mzrange <- swathMS1Coverage( con )
  res <- getMZRTVol( con , ms1id , mzrange = mzrange , rtrange = rtrange)
  return( res )
}

generatePseudospecMS1 = function(fname , coordinates , outdir,  rtext=5 , mzext=4 ,  ms2vol = 500, errorRT = 2.5,
                                 minlength = 5)
{
  con=createConnection(fname)
  addIndexRT2Features(con) # index required for fast query
  createBestFeaturesView(con,rtext,mzext,ms2vol)
  ## generate pseudospectra and write to mgf file.
  fstem=sub("_0.sqlite$","",fname)
  fbase=basename(fstem)
  fbase=paste(fbase,"_erRT",errorRT,"_rtext",rtext,".mgf",sep="")
  fout=file.path(outdir,fbase)
  print(fout)
  print(fbase)
  filecon = file(fout,"w+")
  
  summ <- .generatePseudospec(con,filecon,fbase,coordinates,errorRT=errorRT, minlength=minlength)
  close(filecon)
  tryCatch(
    dbClearResult(dbListResults(con)[[1]]), finally=print("no db result to clear")
  )
  return(list(filename=fbase,summ=summ))
}
##
##
##



#
# generate Pseudospectra from extraction windows in ms2
# iterates over all maps, selects all features in extraction window
# starts to creating pseudospectra
#
pseudospecMS2ExtractWindow = function(
  fname,  
  outdir,
  rtext=5,
  mzext=4,
  ms2volup = 1000, ##parameters for view creation.
  ms2vollow = 200,
  errorRT = 2.5,
  rtrange = NULL,
  minlength=5
){
  con=createConnection(fname)
  addIndexRT2Features(con)
  ## confine search to best features
  ## select features from ms1 map and attempt deisotoping
  maps <- mapSummary(con)
  maps <- maps[maps[,"mslevel"]==2,]
  maps <- maps[order(maps[,"minMZsw"]),]
  maps <- maps[,c("minMZsw","maxMZsw","id")] #make sure it is ordered the same way
  
  fstem=sub("_0.sqlite$","",fname)
  fbase=basename(fstem)
  fbase=paste(fbase,"_erRT",errorRT,"_rtext",rtext,".mgf",sep="")
  fout=file.path(outdir,fbase)
  cat("outfile : ",fout,"\n")
  cat("fbase : " , fbase , "\n")
  filecon = file(fout,"w+")
  
  generSummary <-list(found=0,notfound=0,length=NULL)
  for(i in 1:dim(maps)[1])
  {
    x <- unlist( maps[i,c("id","minMZsw","maxMZsw")])
    ## TODO ADD DEISOTOPING
    createBestFeaturesView(con,rtext,mzext,ms2volup)
    res <- getMZRTVol( con , x["id"] , mzrange = x[2:3] , rtrange = rtrange)
    if(dim(res)[1]>0)
    {
      createBestFeaturesView(con,rtext,mzext,ms2vollow)
      summ <- .generatePseudospec(con, filecon, filename, res, errorRT=errorRT, minlength=minlength )
      generSummary$found = generSummary$found + summ$found
      generSummary$notfound = generSummary$notfound + summ$notfound
      generSummary$length = c(generSummary$length, summ$length)
    }
  }
  close(filecon)
  tryCatch(
    dbClearResult(dbListResults(con)[[1]]), finally=print("no db result to clear")
  )
  #dbDisconnect(con)
  return(generSummary)
}


###
generatePseudopsec = function( res , thrsm = 2 , thrlar = 4, minlength=5 )
{
  result = list()
  count = 1
  while( dim(res)[1] > 0 )
  {
    idxh=which(res$Volume == max(res$Volume))
    idx <- idxh[1]
    highest <- res[idx,] # 
    rt = highest$RT
    fp = highest
    res <- res[-idx,] # remove column 
    
    rtmin = rt - thrsm
    rtmax = rt + thrsm
    idxclose  = which( res$RT > rtmin & res$RT <  rtmax) 
    if(length(idxclose) > 0){
      fp1 = res[idxclose,] 
      fp = rbind(fp, fp1 )
      ## remove from dataset
      res = res[-idxclose,]
    }
    rtmin = rt - thrlar
    rtmax = rt + thrlar
    #cat("rtmin", length(rtmin),"rtmax", length(rtmax),"\n")
    idxfar = which(res$RT> rtmin & res$RT <  rtmax) 
    if(length(idxfar)>0){
      fpfar = res[idxfar,]
      fp = rbind(fp,fpfar)
    }
    if(dim(fp)[1]>minlength){
      result[[count]]=fp
      count = count + 1
    }
  }
  return(result)  
}

pseudospecMS2 = function(
  fname,  
  outdir,
  rtext=5,
  mzext=4,
  ms2vol = 200,
  errorRT = 2.5,
  rtrange = NULL,
  minlength=5
){
  con = createConnection( fname )
  addIndexRT2Features( con )
  ## confine search to best features
  ## select features from ms1 map and attempt deisotoping
  maps <- mapSummary(con)
  maps <- maps[maps[,"mslevel"]==2,]
  maps <- maps[order(maps[,"minMZsw"]),]
  maps <- maps[,c("minMZsw","maxMZsw","id")] #make sure it is ordered the same way
  
  fstem=sub("_0.sqlite$","",fname)
  fbase=basename(fstem)
  fbase=paste(fbase,"_erRT",errorRT,"_rtext",rtext,".mgf",sep="")
  fout=file.path(outdir,fbase)
  cat("outfile : ",fout,"\n")
  cat("fbase : " , fbase , "\n")
  filecon = file(fout,"w+")
  
  generSummary <-list(found=0,notfound=0,length=NULL)
  for(i in 1:dim(maps)[1])
  {
    x <- unlist( maps[i,c("id","minMZsw","maxMZsw")])
    ## TODO ADD DEISOTOPING
    createBestFeaturesView( con , rtext , mzext , ms2vol )
    res <- getMZRTVol( con , x["id"] , rtrange = rtrange )
    if(dim(res)[1]>0)
    {
      spec = generatePseudopsec(res,thrsm=errorRT, thrlar=2*errorRT)
    }
  }
  close(filecon)
  tryCatch(
    dbClearResult(dbListResults(con)[[1]]), finally=print("no db result to clear")
  )
  #dbDisconnect(con)
  return(generSummary)
}
