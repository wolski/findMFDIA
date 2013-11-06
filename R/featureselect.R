## author : Witold Wolski wewolski@gmail.com

##
## This files contains functions for RT peak similarity scoring.
##


###making KURTOSIS less volume dependent and scaling
transformKURT <- function(yk,volume)
    {
        yk <- yk-min(yk)+1
        yk2 <- log(yk)/(log(volume)+5)
        yk2 <- scale(yk2)
        return(yk2)
    }

transformSD <- function(sdk,rtext,volume)
    {
        res <- (sdk/rtext)*(log(volume)+5)
        return(res)
    }


computeScore4Group=function(group, transition)
{
    parent <- group$parent
    prod <- group$prod
    ## [1] "rep(2, length(w))" "id"                "idmap"            
    ## [4] "idswath"           "centerOfMassMZ"    "centerOfMassRT"   
    ## [7] "MZ"                "RT"                "MZSD"             
    ## [10] "RTSD"              "MZSKEW"            "RTSKEW"           
    ## [13] "MZKURT"            "RTKURT"            "Max"              
    ## [16] "Count"             "Volume"            "maximumLocationMZ"
    ## [19] "maximumLocationRT" "minMZIndex"        "mzExtend"         
    ## [22] "minRTIndex"        "rtExtend"          "mzProjection"     
    ## [25] "rtProjection"     
    volumeMS2 <- sum(prod$Volume)
    volumeMS1 <- parent$Volume
    sdkurt <- sd(prod$RTKURT)
    sdsd <- sd(prod$RTSD)
    sdskw <- sd(prod$RTSKEW)
    sdRT <- sd(prod$RT)
    sdcomRT <- sd(prod$centerOfMassRT)
}
