## author : Witold Wolski wewolski@gmail.com
## reading transition csv
## handling transitions

.getTransitionsCSV <- function(filename){
  return(read.table(filename,sep="\t",header=T))
}

createProduct = function(name, mz, intensity, annotation){
    tmp <- cbind(name,mz,intensity,annotation)
    class(tmp)="products"
    return(tmp)
}

createTransition <- function(name, mz, rt, sequence, proteinname, charge, decoy,products){
    transition <- list(
        name=name,
        mz=mz,
        rt=rt,
        sequence=sequence,
        decoy=decoy,
        products=products
        )
    class(transition)<-"precursor"
    return(transition)
}

getTransitions <- function(filename){
    dfcsv <- .getTransitionsCSV(filename)
    res<-.dataframe2transitions(dfcsv)
    res <- sortTransitionsbyRT(res)
    return(res)
}
    

##sort transitions by RT
sortTransitionsbyRT <- function(trans){
    rts <- NULL
    for(i in 1:length(trans))
        {
            rts=c(rts,trans[[i]]$rt)
        }
    trans <- trans[order(rts)]
    class(trans) <- "transitionlist"
    return(trans)
}

.checkReqNames <- function(req,nam){
    req=tolower(req)
    nam=tolower(nam)
    x = req %in% nam
    if(  sum(x) != length(req)){
        stop(paste("required column : ", req[!x]))
    }
}

##remove low mass transitions
filterProdByMass <- function(trans, minmass){
    for( i in 1:length(trans) )
        {
            tmp <- trans[[i]]$products
            trans[[i]]$products <- tmp[tmp[,"mz"]>minmass,]
        }
    return(trans)
}

##maps the retention time to dataset
## model must be model <- lm(irts~reft)
mapRT <- function(trans,model){
for( i in 1:length(trans) )
        {
            tmp <- trans[[i]]$rt
            trans[[i]]$rt <- predict(model,data.frame(reft=tmp))[[1]]
        }
    return(trans)
}

##converts dataframe to transitions
.dataframe2transitions <- function(frame){
    trid <- unique(frame$transition_group_id)
    res <- vector(mode="list",length(trid))
    count <- 1

    nam <- tolower(names(frame))
    req <- c("tr_calibrated","precursormz","peptidesequence","decoy")
    .checkReqNames(req,nam)
    req <- c("transition_name","productmz","libraryintensity","annotation")
    .checkReqNames(req,nam)

    for(i in trid){
        trans<-frame[frame$transition_group_id==i,]
        prod=createProduct(
            name = trans$transition_name ,
            mz = trans$ProductMz,
            intensity = trans$LibraryIntensity,
            annotation = trans$Annotation
            )
        
        res[[count]] = createTransition(name = as.character(trans$ProteinName[1]),
               mz = trans$PrecursorMz[1],
               rt = trans$Tr_calibrated[1] ,
               sequence = as.character(trans$PeptideSequence[1]),
               decoy = trans$decoy[1],
               products = prod )
        count = count +1
    }
    return(res)
    
}

#todo expand also the products
as.data.frame.transitionlist = function(precursor){
    res <- matrix(0, nrow=length(precursor),ncol=5)
    res <- as.data.frame(res)
    for(i in 1:length(precursor) ){
       row <- precursor[[i]][1:5]
       nam <- (names(row))
       res[i,]<-row
    }
    colnames(res) <- nam
    return(res)
}
