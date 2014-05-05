library("RSQLite")
library("DBI")
#' creates database connection
#' @export
createConnection=function(dbstring){
  drv <-dbDriver("SQLite")
  con <- dbConnect(drv,dbstring,cache_size=10000)
}
#' preselects features for querying
#' @export
createBestFeaturesView=function(con, rtext=2,mzext=2,Volume=1000,name="bestfeatures"){
  dbSendQuery(con,paste("Drop VIEW if exists ", name, sep=""))
  query<-paste("Create VIEW ", name, " AS select * from features where rtExtend > ",
               rtext, " and  mzExtend> ", mzext, " and Volume> ", Volume ,sep="")
  print(query)
  dbSendQuery(con,query)
}
#' adds MZ RT index to features table
#' @export
addIndexMZRT2Features = function(con){
  dbSendQuery(con,"Drop index if exists idswathMZRTfeatures;")
  dbSendQuery(con,"Create index if not exists idswathMZRTONfeatures ON features (idswath,MZ,RT);")
}
#' add RT index to features
#' @export
addIndexRT2Features = function(con){
  ## dbSendQuery(con,"Drop index if exists idswathMZRTfeatures;")
  dbSendQuery(con,"Create index if not exists idswathRTONfeatures ON features (idswath,RT);")
}
#' drop MZ RT index
#' @export
dropIndexMZRTFeatures = function(con){
  dbSendQuery(con,"Drop index if exists idswathMZRTfeatures;")
}
#' newwindows matrix or df
#' @param con database connection
#' @param newwindows must have 2 columns named minMZsw and maxMZsw
updateSwathWindows=function(con, newwindows){
  # sort newwindows by minmass
  # ret
  maps <- mapSummary(con)
  maps <- maps[maps[,"mslevel"]==2,]
  maps <- maps[order(maps[,"minMZsw"]),]
  #make sure it is ordered the same way
  newwindows <- newwindows[newwindows[,"mslevel"]==2,]
  newwindows <- newwindows[order(newwindows[,"minMZsw"]),]
  if(dim(maps)[1] != dim(newwindows)[1]){
    stop("dims do not agree")
  }
  
  dbBeginTransaction(con)
  for(i in 1:dim(maps)[1]){
    query=paste("update mapinfo set minMZsw=" ,
                newwindows[i,"minMZsw"], " ,maxMZsw=" , newwindows[i,"maxMZsw"],
                " where id = " , maps[i,"id"],
                sep=""
    )
    print(query)
    rs <- dbSendQuery(con,query)
  }
  dbCommit(con)
  
  return(rs);
}
