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

## Write mgf files
.writeMGF = function(fcon,filename,parpeak,peaks,specid=1){
  parpeak = parpeak
  peaks = peaks[order(peaks["MZ"]),] # sort by MZ
  
  cat("BEGIN IONS\n", file=fcon)
  cat("TITLE=",filename,".",specid,"\n",file=fcon,sep="")
  cat("PEPMASS=",parpeak["MZ"]," ",parpeak["Volume"],"\n",file=fcon,sep="")
  cat("RTINSECONDS=",parpeak["RT"],"\n",file=fcon,sep="")
  if(!is.nan(parpeak["Charge"])){
    #cat("CHARGE=",parpeak["Charge"],"+\n",file=con,sep="")
  }
  for(i in 1:dim(peaks)[1])
  {
    cat(peaks[i,"MZ"], peaks[i,"Volume"], "\n", file=fcon, sep=" ")
  }
  cat("END IONS\n\n",file=fcon,sep="")
}

## generate Pseudospectrum given ms1 peaks
.generatePseudospec = function(
  dbcon ,
  filecon ,
  filename ,
  peaks ,
  errorRT=1.5 , # error of RT in [s]
  minlength = 5
)
{
  generSummary <-list(found=0,notfound=0,length=NULL)
  for(i in 1:dim(peaks)[1])
  {
    tmp = .getPseudoSpec(dbcon,peaks[i,] , errorRT = errorRT )
    if(dim(tmp)[1]>minlength)
    {
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

## use this function to retrieve ms1 features
getMS1seeds = function( con , rtext=3 , mzext=3 , ms1vol = 500, rtrange=NULL )
{
  ### Confine search to best features
  createBestFeaturesView(con,rtext,mzext,ms1vol)
  ms1id <- getMS1id(con)
  mzrange <- swathMS1Coverage( con )
  res <- getMZRTVol( con , ms1id , mzrange = mzrange , rtrange = rtrange)
  return( res )
}


#' Generate pseudospectra using MS1 seeds
#'
#' @param  coordinates dataframe with peaks to use to generate pseudospectra from
#' @param fname sql file name
#' @param outdir directory to write to
#' @param rtext feature extend in rt
#' @param mzext feature extend in mz
#' @param ms2vol minimum volume of feature
#' @param errorRT retention time error 
#' @param minlength minimum peaklist length
#' 
generatePseudospecMS1 = function(fname ,
                                 coordinates , 
                                 outdir ,
                                 rtext=5 ,
                                 mzext=4 ,
                                 ms2vol = 500 ,
                                 errorRT = 2.5 ,
                                 minlength = 5
)
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
  
  summ <- .generatePseudospec(con,filecon,fbase,peaks=coordinates,errorRT=errorRT, minlength=minlength)
  close(filecon)
  tryCatch(
    dbClearResult(dbListResults(con)[[1]]), finally=print("no db result to clear")
  )
  return( list(filename=fbase,summ=summ) )
}
##
##
##




