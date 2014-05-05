#' file stem
#' @export
file.stem <- function(x){
        bn <- basename(x)
        gsub("\\..*$", "", bn)
}
