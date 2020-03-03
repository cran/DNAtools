#' Estimate theta and the fraction of comparisons between close relatives
#' 
#' Estimates the fraction of comparisons between pairs of close relatives while
#' fitting the theta parameter minimising the object function. The function
#' makes use of the R-package 'Rsolnp' which is an implementation of an solver
#' for non-linear minimisation problems with parameter constraints.
#' 
#' Computes the proportion of comparisons between close relatives in a database
#' matching exercise for each theta value under investigation.
#' 
#' @param obs The matrix or vector of observed matches/partial-matches as
#' returned by the dbCompare()-function
#' @param theta0 The left value of the interval in which a bisection-like
#' search is performed for theta
#' @param theta1 Right value of interval (see theta0)
#' @param theta.tol A stopping criterion for the search. If the search narrows
#' within theta.tol the function terminates
#' @param theta.step Default is NULL. If not a grid search will be performed on
#' seq(from = theta0, to = theta1, by = theta.step)
#' @param max.bisect The maximum number of bisectional iterations perform prior
#' to termination
#' @param probs List of vectors with allele probabilities for each locus
#' @param var.list A named list of components for computing variances, see
#' dbVariance. The names of the elements are the associated theta-values, and
#' each component is a list of (V1,V2,V3) - see dbVariance with n=1
#' @param init.alpha Initial values for alpha, where the order is
#' (First-cousins, Avuncular, Parent-child, Full-siblings). The value for
#' Unrelated is computed as 1-sum(init.alpha)
#' @param init.keep Whether the initial values should be used in successive
#' steps for the current optimum should be used.
#' @param objFunction Which of the five different object functions should be
#' used to compare observed and expected
#' @param collapse Not yet implemented
#' @param trace Should iteration steps and other process indicators be printed
#' @param solnp.ctrl See solnp for details
#' @return Returns a list of three components: value, solution and var.list.
#' The first element, value, is a dataframe with the value of the objection
#' function for each of the theta values investigated. Solution is the
#' estimated alpha-vector where the objection function was minimised. Finally,
#' var.list is a names list of components for computing variances. May be
#' reused in later computations for increased speed in some iterations.
#' @author James Curran and Torben Tvedebrink
#' @references T Tvedebrink, PS Eriksen, J Curran, HS Mogensen, N Morling.
#' 'Analysis of matches and partial-matches in Danish DNA reference profile
#' database'. Forensic Science International: Genetics, 2011.
#' @examples
#' 
#'   \dontrun{
#'   ## Simulate some allele frequencies:
#'   freqs <-  replicate(10, { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
#'               simplify=FALSE)
#'   ## Load the sample database:
#'   data(dbExample)
#'   obs <- dbCompare(dbExample,trace=FALSE)$m
#'   C3 <- optim.relatedness(obs,theta0=0.0,theta1=0.03,probs=freqs,
#'           objFunction='C3',max.bisect=30,trace=TRUE)
#'   }
#' 
#' @importFrom Rsolnp solnp
#' 
#' @export optim.relatedness
optim.relatedness <- function(obs, theta0 = 0, theta1 = 0.03, theta.tol = 10^(-7), theta.step = NULL, 
  max.bisect = 15, probs, var.list = NULL, init.alpha = 10^c(-4, -6, -8, -10), init.keep = FALSE, 
  objFunction = c("T2", "T1", "C3", "C2", "C1"), collapse = FALSE, trace = FALSE, solnp.ctrl = list(tol = 10^(-9), 
    rho = 10, delta = min(init.alpha) * 0.01, trace = FALSE)) {
  if (!any(objFunction == c("T2", "T1", "C3", "C2", "C1"))) 
    stop("Wrong objFunction is supplied. Use any of: \"T2\", \"T1\", \"C3\", \"C2\" or \"C1\")")
  if (sum(init.alpha) > 1 | any(init.alpha < 0) | any(init.alpha > 1)) 
    stop("init.alpha is a probability vector, i.e. 0<=alpha[i]<=1 and sum(alpha)<=1")
  objFunction <- objFunction[1]
  n <- (1 + sqrt(1 + 8 * sum(obs, na.rm = TRUE)))/2
  if (is.matrix(obs)) 
    obs <- t(obs)[up.tri(obs)]
  if (collapse) {
    ## Collapse the M_{m/p} matrix into matrix of matching alleles
    warning("attempt to 'collapse' observations into mathcing alleles ignored. The options is not yet implemented.")
  }
  
  ## Declaration of the different objective functions
  C1 <- function(x, expected, observed, variance = NULL) {
    EE <- x[1] * expected[[1]]
    for (r in 2:length(expected)) EE <- EE + x[r] * expected[[r]]
    sum(((observed - EE)^2/EE), na.rm = TRUE)
  }
  C2 <- function(x, expected, observed, variance = NULL) {
    EE <- x[1] * expected[[1]]
    for (r in 2:length(expected)) EE <- EE + x[r] * expected[[r]]
    sqrt(sum((observed - EE)^2, na.rm = TRUE))
  }
  C3 <- function(x, expected, observed, variance = NULL) {
    EE <- x[1] * expected[[1]]
    for (r in 2:length(expected)) EE <- EE + x[r] * expected[[r]]
    sum(abs(observed - EE)/EE, na.rm = TRUE)
  }
  T1 <- function(x, expected, observed, variance = NULL) {
    EE <- x[1] * expected[[1]]
    for (r in 2:length(expected)) EE <- EE + x[r] * expected[[r]]
    sum((EE - observed)^2/diag(variance))
  }
  T2 <- function(x, expected, observed, variance = NULL) {
    EE <- x[1] * expected[[1]]
    for (r in 2:length(expected)) EE <- EE + x[r] * expected[[r]]
    gginv <- eigen(variance, T, F, T)
    ivariance <- gginv$vec %*% diag(c(1/gginv$val[-length(gginv$val)], 0)) %*% t(gginv$vec)
    as.numeric(t(EE - observed) %*% ivariance %*% (EE - observed))
  }
  objFun <- get(objFunction)
  ## Now the objFun contains the function declaration for the object function needed in the
  ## optimisation steps
  
  tgrid <- seq(from = theta0, to = theta1, len = 3)  ## Initial theta grid
  val.df <- data.frame(theta = NA, value = NA)[0, ]  ## The output containing the theta values and 
  min.t <- rep(tgrid[2], 2)
  
  ## GRID SEARCH
  if (!is.null(theta.step)) {
    tgrid <- seq(from = theta0, to = theta1, by = theta.step)
    grid.search <- TRUE
  } else grid.search <- FALSE
  
  bb <- 0
  min.val <- Inf
  min.id <- 1
  if (sum(init.alpha) != 1) 
    init.alpha <- c(1 - sum(init.alpha), init.alpha)
  min.res <- init.alpha
  if (trace) {
    if (grid.search) 
      cat("Grid search... ") else cat("Bisectional search...\n")
  }
  while (bb < max.bisect) {
    bb <- bb + 1  ## number of bisection searches performed
    if (trace & !grid.search) 
      cat(paste("Iteration: ", format(bb, width = 2), "\n", sep = ""))
    if (diff(range(tgrid)) < theta.tol & !grid.search) {
      if (trace) {
        if (solnp.ctrl$trace) 
          cat("\n")
        cat("Interval converged\n")
      }
      break
    }
    if (length(min.t) > 10) {
      if (abs(diff(min.t[length(min.t) - c(0, 5)])) < (theta.tol)^2 & !grid.search) {
        if (trace) 
          cat("No change in theta for several iterations\n")
        break
      }
    }
    if (grepl("T", objFunction)) {
      ## Object function is either T2 or T1 - variances are needed..
      if (is.null(var.list)) {
        if (trace) 
          cat("Variances are being computed... Please wait")
        var.list <- dbVariance(probs = probs, theta = tgrid, n = 1)
        names(var.list) <- paste(tgrid)
        if (trace) 
          cat("  Done..!\n")
      } else {
        if (all(is.element(paste(tgrid), names(var.list))) & trace) 
          cat("All needed variances are provided in the input...\n") else {
          if (trace) 
          cat(paste("Missing variances (", sum(!is.element(paste(tgrid), names(var.list))), 
            ") are being computed... Please wait", sep = ""))
          ttgrid <- tgrid[!is.element(paste(tgrid), names(var.list))]
          vvar.list <- dbVariance(probs, theta = ttgrid, n = 1)
          if (length(ttgrid) == 1) 
          vvar.list <- list(vvar.list)
          names(vvar.list) <- ttgrid
          var.list <- c(var.list, vvar.list)
          var.list <- var.list[sort.list(names(var.list))]
          if (trace) 
          cat("  Done..!\n")
        }
      }
      variances <- lapply(var.list, function(x, n) choose(n, 2) * x$V1 + 6 * choose(n, 
        3) * x$V2 + 6 * choose(n, 4) * x$V3, n = n)
      variances <- variances[paste(tgrid)]
    } else {
      ## Object function is either C1, C2 or C3
      var.list <- NULL
      variances <- replicate(length(tgrid), NULL, simplify = FALSE)
    }
    expects <- lapply(tgrid, function(t, n) list(UN = dbExpect(probs = probs, theta = t, 
      n = n, vector = TRUE, k = c(0, 0, 1)), FC = dbExpect(probs = probs, theta = t, n = n, 
      vector = TRUE, k = c(0, 1, 3)/4), AV = dbExpect(probs = probs, theta = t, n = n, 
      vector = TRUE, k = c(0, 1, 1)/2), PC = dbExpect(probs = probs, theta = t, n = n, 
      vector = TRUE, k = c(0, 1, 0)), FS = dbExpect(probs = probs, theta = t, n = n, vector = TRUE, 
      k = c(1, 2, 1)/4)), n = n)
    ## Part of the actual computations for each value of theta in the bisection or grid search
    for (i in 1:length(tgrid)) {
      t <- tgrid[i]  ## current theta value under consideration
      solnpObjFun <- function(x) objFun(x, expected = expects[[i]], observed = obs, variance = variances[[i]])
      alpha <- min.res
      if (init.keep) 
        alpha <- init.alpha
      est <- try(Rsolnp::solnp(pars = alpha, fun = solnpObjFun, eqfun = sum, eqB = 1, LB = rep(0, 
        5), UB = rep(1, 5), control = solnp.ctrl), silent = TRUE)
      if (length(est) == 1) {
        ## 'est' contains the error message produced by the solnp function
        est <- list(pars = init.alpha, values = NA)
        val.df <- rbind(val.df, data.frame(theta = t, value = NA))
        warn.message <- paste("NAs were returned for theta =", format(t, digits = 5), 
          "by the solnp procedure. This could indicate numerical problems with the identified solution.", 
          sep = " ")
        if (!grid.search) 
          warn.message <- paste(warn.message, "\n  An attempt to overcome or approximate/solve the problem is to use a grid search instead using 'theta.step' to set step size.")
        warning(warn.message)
      } else {
        step.val <- est$values[length(est$values)]
        if (step.val <= min.val) {
          min.val <- step.val
          min.res <- est$pars
          min.id <- i
          min.t <- c(min.t, t)
        }
        val.df <- rbind(val.df, data.frame(theta = t, value = step.val))
      }
    }
    ## 
    if (grid.search) {
      if (trace) 
        cat("\n")
      break  ## breaks out of bisectional search loop
    }
    ## The mid point t1 is the minima of [t0,t1,t2].. Update to search in [t0+0.5*l,t2-0.5*l],
    ## where l = t1-t0 = t2-t1
    if (min.id == 1) 
      tgrid <- seq(from = tgrid[1], to = tgrid[2], len = 3) else if (min.id == length(tgrid)) 
      tgrid <- seq(from = tgrid[min.id - 1], to = tgrid[min.id], len = 3) else {
      tgrid.length <- c(tgrid[min.id] - min(tgrid), max(tgrid) - tgrid[min.id])
      tgrid <- sort(c(tgrid[1], tgrid[min.id], tgrid[length(tgrid)], tgrid[min.id] + 0.25^bb * 
        tgrid.length[1], tgrid[min.id] - 0.25^bb * tgrid.length[2]))
    }
  }
  names(min.res) <- c("Unrelated", "First-Cousins", "Avuncular", "Parent-child", "Full-siblings")
  val.df <- val.df[!duplicated(val.df$theta), ]
  res <- list(value = val.df[order(val.df$theta), ], solution = c(theta = tgrid[min.id], min.res), 
    var.list = var.list)
  attributes(res)$objFun <- objFunction
  attributes(res)$class <- "dbOptim"
  res
}



#' Plots the fitted object function for estimated familial relationships in the
#' database and theta.
#' 
#' Plots the minimised object function for included values of theta
#' 
#' Plots the object function
#' 
#' @aliases plot.dbOptim points.dbOptim lines.dbOptim
#' @param x Object returned by optim.relatedness
#' @param type The type of plot character ('l'=line, 'p'=points, ...), see
#' 'par' for more details
#' @param ... Other plot options
#' @return A plot of the object function
#' @author James Curran and Torben Tvedebrink
#' @seealso optim.relatedness
#' @examples
#' 
#'   \dontrun{
#'   ## Simulate some allele frequencies:
#'   freqs <-  replicate(10, { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
#'               simplify=FALSE)
#'   ## Load the sample database:
#'   data(dbExample)
#'   obs <- dbCompare(dbExample,trace=FALSE)$m
#'   C3 <- optim.relatedness(obs,theta0=0.0,theta1=0.03,probs=freqs,
#'           objFunction='C3',max.bisect=30,trace=TRUE)
#'   plot(C3)
#'   }
#' 
#' @rawNamespace S3method(plot, dbOptim)
#' @importFrom graphics plot
#' @export plot.dbOptim
plot.dbOptim <- function(x, type = "l", ...) {
  objFun <- attributes(x)$objFun
  ylabel <- switch(objFun, C1 = expression(C[1](theta)), C2 = expression(C[2](theta)), C3 = expression(C[3](theta)), 
    T1 = expression(T[1](theta)), T2 = expression(T[2](theta)))
  graphics::plot(value ~ theta, x$value, xlab = expression(theta), ylab = ylabel, type = type, ...)
}

#' @importFrom graphics points
points.dbOptim <- function(x, type = "p", ...) {
  graphics::points(value ~ theta, x$value, type = type, ...)
}

#' @importFrom graphics points
lines.dbOptim <- function(x, type = "l", ...) {
  graphics::points(value ~ theta, x$value, type = type, ...)
}



#' Prints the results from optim.relatedness()
#' 
#' Prints the evaluated functions for the object function, best estimate of
#' alpha and possibly list of variances.
#' 
#' Prints the summary details of the fit
#' 
#' @param x Object returned by optim.relatedness()
#' @param var.list Logical. Whether the (long) list of variance components
#' should be printed to the screen.
#' @param ... ...
#' @return A dataframe with [theta,value] and a vector of fitted alpha
#' parameters
#' @author James Curran and Torben Tvedebrink
#' @seealso optim.relatedness
#' @examples
#' 
#'   \dontrun{
#'   ## Simulate some allele frequencies:
#'   freqs <-  replicate(10, { g = rgamma(n=10,scale=4,shape=3); g/sum(g)},
#'               simplify=FALSE)
#'   ## Load the sample database:
#'   data(dbExample)
#'   obs <- dbCompare(dbExample,trace=FALSE)$m
#'   C3 <- optim.relatedness(obs,theta0=0.0,theta1=0.03,probs=freqs,
#'           objFunction='C3',max.bisect=30,trace=TRUE)
#'   print(C3)
#'   }
#' @rawNamespace S3method(print, dbOptim)
#' @export print.dbOptim
print.dbOptim <- function(x, var.list = FALSE, ...) {
  if (var.list) 
    print(x[c("value", "solution", "var.list")]) else print(x[c("value", "solution")])
}
