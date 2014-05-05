###
### functions to generate pseudospectra from MS2 maps only without precursor information.
###
## write mgf files ##
.writeMGFSwath = function(fcon,filename,parpeak,peaks,specid=1){
  parpeak<-parpeak
  cat("BEGIN IONS\n", file=fcon)
  cat("TITLE=",filename,".",specid,"\n",file=fcon,sep="")
  cat("PEPMASS=",parpeak["MZl"]," ",parpeak["MZh"]," ",parpeak["Volume"],"\n",file=fcon,sep="")
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
###  ...
.writeSpecList = function( fcon , filename , parmz , peaklists )
{
  if(length(peaklists)<1)
  {
    return(0)
  }
  parpeak=(c(id = 0, RT = 0, MZl =0, MZh = 0, Volume=0))
  for(i in 1:length(peaklists))
  {
    pp = peaklists[[i]][1,] # take most intense peak
    parpeak["MZl"] = parmz[1] # replace with swath extraction window center
    parpeak["MZh"] = parmz[2] # replace with swath extraction window center
    parpeak["id"] = pp["id"]
    parpeak["RT"] = pp["RT"]
    parpeak["Volume"] = pp["Volume"]
    parpeak <- unlist(parpeak)
    #print(class(parpeak))
    #print(parpeak)
    pl <- peaklists[[i]]
    pl <- pl[order(pl["MZ"]),] # sort by MZ
    .writeMGFSwath(fcon,filename,unlist(parpeak),pl,specid = i)
  }
}
#' Generate pseudospec from map
#'
#' @param res a data frame, must contain columns Volume and RT 
#' @param thrsm lower rt threshold
#' @param thrlar upper rt threshold
#' @param minlength minimum spectrum length
#' @return list of spectra 
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
#' generate ms2 pseudospectra
#' 
#' @export
pseudospecMS2 = function(
  fname,  
  outdir,
  rtext=5,
  mzext=4,
  ms2vol = 200,
  errorRT = 2.5, # allowed RT error
  rtrange = NULL, # rtrange of peaklist
  minlength=5 # minimum length of peaklist
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
  createBestFeaturesView( con , rtext , mzext , ms2vol )

  generSummary <-list(found=0,notfound=0,length=NULL)
  for(i in 1:dim(maps)[1])
  {
    # get the swath window
    map <- unlist( maps[i,c("id","minMZsw","maxMZsw")])
    ## TODO ADD DEISOTOPING
    res <- getMZRTVol( con , map["id"] , rtrange = rtrange )
    print(dim(res))
    if( dim(res)[1]>0 ) # if there are features in the map ...
    {
      spec = generatePseudopsec( res , thrsm = errorRT , thrlar=2 * errorRT, minlength=minlength )
      writeSpecList(filecon,filename="test.txt",parmz=map[2:3],peaklists=spec)
    }
  }
  close(filecon)
  tryCatch(
    dbClearResult(dbListResults(con)[[1]]), finally=print("no db result to clear")
  )
  #dbDisconnect(con)
  return(generSummary)
}

