#' query the swath of parentmass
#' @export
getSwath=function(con,parentmass){
    query <- paste("select id from mapinfo where maxMZsw > " , parentmass , " and minMZsw <= " , parentmass)
    res<-dbSendQuery(con,query)
    coord = fetch(res,n=-1)
    return(coord$id)
}
#' retrieving features
#' @export
getFeaturesFlex = function( con, idswath, mzrange =NULL, rtrange=NULL , features="bestfeatures" )
{
    query=paste( "select * from ", features , " where idswath = ", idswath , sep="")
    if(length(mzrange)==2){
        query=paste( query, " and  MZ >= ", mzrange[1], " and  MZ < " , mzrange[2],  sep="")
    }
    if(length(rtrange)==2){
        query = paste(query, " and RT >= ", rtrange[1], " and RT < " , rtrange[2] ,sep="" )
    }
    res<-dbSendQuery(con,query)
    coord <- fetch(res,n=-1)
    return(coord)
}
#' Xtract features MZ, RT and Volume from swath and sort by intensity
#' @export
getMZRTVol = function(con, idswath, features="bestfeatures" , rtrange=NULL, mzrange=NULL)
{
    test <- features
    query <- paste("select id, RT, MZ, Volume from ",features , " where idswath = " , idswath )
    if(length(rtrange)==2){
         query <- paste(query," and RT > " , rtrange[1] , " and RT < ", rtrange[2] )
    }
    if(length(mzrange)==2){
        query <- paste(query," and MZ > " , mzrange[1] , " and MZ < ", mzrange[2] ,
        " order by Volume DESC")
    }
    #print(query)
    qres=dbSendQuery(con, query)
    res = fetch(qres, n=-1)
    if(dim(res)[1] == 0){
      res = data.frame("id" = NA,"RT" = NA,"MZ" = NA,"Volume"=NA)
    }
    return(res)
}

#' get a single column from map
#' @export
getFromMap = function(con, idswath, column = "Volume" , features = "bestfeatures" ){
    query=paste( "select ", column, " from ", features , " where idswath = ", idswath ,sep="" )
    res<-dbSendQuery(con,query)
    coord <- fetch(res,n=-1)
    return(coord)
}
#' getID of ms1 map
#' @export
getMS1id= function(con){
    id <- fetch(dbSendQuery(con, "select id from mapinfo where mslevel=1"),n=-1)
    return(id[[1]])
}
#' retrieves the RT axis of an lcms map
#' @export
getRTAxis=function(con,idmap){
    quer <-  paste("select rts from mapinfo where id = ", idmap, sep="")
    #print(quer)
    axis <- fetch(dbSendQuery(con,quer ), n=-1)
    axis <- readBin(axis$rts[[1]],"numeric",size=8,n=20000)
    return(axis)
}
#' summary methods - computes some summaries over the db.
#' @export
summaryDB=function(con)
{
    featurespswath <- fetch(dbSendQuery(con,"select features.idswath , mapinfo.mslevel, mapinfo.minMZsw  , mapinfo.maxMZsw , count(*) as count, sum(features.Volume) as volume from features, mapinfo where features.idswath = mapinfo.id group by features.idswath order by mapinfo.minMZsw;"),n=-1)
    return( featurespswath )
}
#' produces a map summary
#' @export
mapSummary=function(con){
  query=paste("select id,  mslevel, minMZsw, maxMZsw, minMZ, maxMZ, minRT, maxRT from mapinfo order by minMZsw ")
  res<-dbSendQuery(con,query)
  coord <- fetch(res,n=-1)
  return(coord)
}
#' summarize what mass range is covered by swath
#' @export
swathMS1Coverage=function(con){
  query=paste("select min(minMZsw),max( maxMZsw) from mapinfo where mslevel=2")
  res<-dbSendQuery(con,query)
  coord <- fetch(res,n=-1)
  return(coord)
}
#' extract ion chromatogram from map
#' @export
getXIC=function(con,idmap, mzmin,mzmax){
    #print(idmap)
    rtaxis<-getRTAxis(con,idmap)
    tmp = getFeaturesFlex(con,idmap,c(mzmin,mzmax))
    return(cbind(rtaxis,features2Chrom(tmp,rtaxis)))
}
 
