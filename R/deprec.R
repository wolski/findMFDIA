#
# generate Pseudospectra from extraction windows in ms2 ( using mass available there ) 
# iterates over all maps, selects all features in extraction window
# starts to creating pseudospectra
#
# mark as deprecated - wasn't working.
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
