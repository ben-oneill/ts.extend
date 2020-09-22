#' Generate random vectors from the stationary GARMA distribution
#'
#' This function generates random vectors from the stationary Gaussian auto-regressive moving-average (GARMA) distribution.  The user specifies
#' the number of vectors ```n``` and their dimension ```m``` and the function returns an n x m matrix of generated time-series from the GARMA
#' distribution with the specified parameters.  By default the function generates from the marginal GARMA distribution, but the user may give
#' conditional values in the ```condvals``` vector to generate from the associated conditional distribution (non-conditional values in this
#' vector are given as ```NA```).
#'
#' @param n Positive integer giving the number of random vectors to generate
#' @param m Positive integer giving the dimension of the random vectors to generate (i.e., the number of values in each time-series)
#' @param condvals Either a single value ```NA``` or a numeric vector with ```m``` elements; numeric entries are conditioning values for the generated vector
#' @param mean The mean parameter
#' @param errorvar The error variance parameter
#' @param ar Vector of auto-regressive coefficients (all roots of AR characteristic polynomial must be outside the unit circle)
#' @param ma Vector of moving-average coefficients
#'
#' @examples
#' #Set the model parameters
#' AR <- c(0.8, -0.2)
#' MA <- c(0.6,  0.3)
#' #Generate random time-series from the GARMA distribution
#' SERIES <- rGARMA(n = 16, m = 30, ar = AR, ma = MA)
#'
#' #Set the conditional values
#' CONDVALS     <- rep(NA, 30)
#' CONDVALS[1]  <- -4
#' CONDVALS[12] <-  0
#' CONDVALS[30] <-  4
#'
#' #Generate and plot random time-series from the GARMA distribution
#' SERIES.COND <- rGARMA(n = 16, m = 30, ar = AR, ma = MA, condvals = CONDVALS)

rGARMA <- function(n,
                   m,
                   condvals = as.numeric(NA),
                   mean     = 0,
                   errorvar = 1,
                   ar       = numeric(0),
                   ma       = numeric(0)) {

  #Check inputs
  if (!is.numeric(n))            { stop('Error: n should be a positive integer')  }
  if (length(n) != 1)            { stop('Error: n should be a positive integer')  }
  if (n != as.integer(n))        { stop('Error: n should be a positive integer')  }
  if (n < 1)                     { stop('Error: n cannot be less than one')  }
  if (!is.numeric(m))            { stop('Error: m should be a positive integer')  }
  if (length(m) != 1)            { stop('Error: m should be a positive integer')  }
  if (m != as.integer(m))        { stop('Error: m should be a positive integer')  }
  if (m < 1)                     { stop('Error: m cannot be less than one')  }
  if (!is.numeric(mean))         { stop('Error: mean must be numeric') }
  if (length(mean) != 1)         { stop('Error: mean must be a single numeric value') }
  if (!is.numeric(errorvar))     { stop('Error: errorvar must be numeric') }
  if (length(errorvar) != 1)     { stop('Error: errorvar must be a single numeric value') }
  if (errorvar < 0)              { stop('Error: errorvar cannot be negative') }

  #Check input condvals
  if (!is.vector(condvals))      { stop('Error: condvals must be a conditioning vector') }
  if (length(condvals) == 1) {
    if (is.na(condvals))       {
      condvals <- as.numeric(rep(NA, m)); } }
  mm   <- length(condvals);
  if (mm != m)                   { stop('Error: condvals must have length m') }
  cond <- !is.na(condvals);
  cc   <- sum(cond);

  #Deal with trivial case where all elements are conditioning elements
  if (cc == m) {
    OUT <- matrix(condvals, nrow = n, ncol = m, byrow = TRUE);
    class(OUT)    <- 'time.series';
    rownames(OUT) <- sprintf('Series[%s]', 1:n);
    colnames(OUT) <- sprintf('Time[%s]', 1:m);
    return(OUT); }

  #Compute mean vector and variance matrix
  MEAN <- rep(mean, m);
  VAR  <- errorvar*ARMA.var(n = m, ar = ar, ma = ma);

  #Compute conditional mean vector and conditional variance matrix
  if (cc == 0) {
    CMEAN  <- MEAN;
    CVAR   <- VAR; }
  if (cc >  0) {
    MEAN1  <- MEAN[!cond, drop = FALSE];
    MEAN2  <- MEAN[ cond, drop = FALSE];
    VAR11  <- VAR[!cond, !cond, drop = FALSE];
    VAR12  <- VAR[!cond,  cond, drop = FALSE];
    VAR21  <- VAR[ cond, !cond, drop = FALSE];
    VAR22  <- VAR[ cond,  cond, drop = FALSE];
    CMEAN  <- MEAN1 + VAR12 %*% solve(VAR22, condvals[cond] - MEAN2);
    CVAR   <- VAR11 - VAR12 %*% solve(VAR22, VAR21); }

  #Generate random vectors
  OUT  <- matrix(condvals, nrow = n, ncol = m, byrow = TRUE);
  if (cc <  m) {
    OUT[, !cond] <- mvtnorm::rmvnorm(n, mean = CMEAN, sigma = CVAR); }

  #Add class and labels for rows and columns
  class(OUT)    <- 'time.series';
  rownames(OUT) <- sprintf('Series[%s]', 1:n);
  colnames(OUT) <- sprintf('Time[%s]', 1:m);

  OUT; }
