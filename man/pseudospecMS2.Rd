\name{pseudospecMS2}
\alias{pseudospecMS2}
\alias{writeSpecList}
\alias{generatePseudopsec}
\alias{pseudospecMS2}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating pseudospectra without precursor
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pseudospecMS2(fname, outdir, rtext = 5, mzext = 4, ms2vol = 200, errorRT = 2.5, rtrange = NULL, minlength = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
%%     ~~Describe \code{fname} here~~
}
  \item{outdir}{
%%     ~~Describe \code{outdir} here~~
}
  \item{rtext}{
%%     ~~Describe \code{rtext} here~~
}
  \item{mzext}{
%%     ~~Describe \code{mzext} here~~
}
  \item{ms2vol}{
%%     ~~Describe \code{ms2vol} here~~
}
  \item{errorRT}{
%%     ~~Describe \code{errorRT} here~~
}
  \item{rtrange}{
%%     ~~Describe \code{rtrange} here~~
}
  \item{minlength}{
%%     ~~Describe \code{minlength} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fname, outdir, rtext = 5, mzext = 4, ms2vol = 200, 
    errorRT = 2.5, rtrange = NULL, minlength = 5) 
{
    con = createConnection(fname)
    addIndexRT2Features(con)
    maps <- mapSummary(con)
    maps <- maps[maps[, "mslevel"] == 2, ]
    maps <- maps[order(maps[, "minMZsw"]), ]
    maps <- maps[, c("minMZsw", "maxMZsw", "id")]
    fstem = sub("_0.sqlite$", "", fname)
    fbase = basename(fstem)
    fbase = paste(fbase, "_erRT", errorRT, "_rtext", rtext, ".mgf", 
        sep = "")
    fout = file.path(outdir, fbase)
    cat("outfile : ", fout, "\n")
    cat("fbase : ", fbase, "\n")
    filecon = file(fout, "w+")
    generSummary <- list(found = 0, notfound = 0, length = NULL)
    for (i in 1:dim(maps)[1]) {
        map <- unlist(maps[i, c("id", "minMZsw", "maxMZsw")])
        createBestFeaturesView(con, rtext, mzext, ms2vol)
        res <- getMZRTVol(con, map["id"], rtrange = rtrange)
        if (dim(res)[1] > 0) {
            spec = generatePseudopsec(res, thrsm = errorRT, thrlar = 2 * 
                errorRT, minlength = minlength)
            writeSpecList(filecon, filename = "test.txt", parmz = map[2:3], 
                peaklists = spec)
        }
    }
    close(filecon)
    tryCatch(dbClearResult(dbListResults(con)[[1]]), finally = print("no db result to clear"))
    return(generSummary)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
