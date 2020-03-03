## Makes all permuations of a vector and returns a matrix // Based on the 'multicool' package
## of Prof. James M. Curran
#' @importFrom multicool initMC allPerm
permAll <- function(x) {
  if (length(x) == 1) 
    return(x)
  xx = multicool::initMC(x)
  multicool::allPerm(xx)
}
