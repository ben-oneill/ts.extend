#' Cumulative distribution function for the stationary GARMA distribution
#'
#' This function computes the cumulative distribution function from the stationary Gaussian auto-regressive moving-average (GARMA) distribution.
#' The user specifies a vector \code{x} giving a single time-series vector, or a matrix \code{x} giving one time-series vector in each row, and the
#' function returns the vector of cumulative probabilities corresponding to the input time-series vectors.  By default the function generates from
#' the marginal GARMA distribution, but the user may give conditioning indicators in the \code{cond} vector to compute the conditional density where
#' some of the elements in the input vectors are conditioning values.
#'
#' @param x A vector or matrix of time-series values (if a matrix, each time-series should be one row of the matrix)
#' @param cond Either a single logical value \code{FALSE} or a logical vector with the same number of elements; as each time-series vector; each logical value indicates whether the density is conditional on the associated time-series value in \code{x}.
#' @param mean The mean parameter
#' @param errorvar The error variance parameter
#' @param ar Vector of auto-regressive coefficients (all roots of AR characteristic polynomial must be outside the unit circle)
#' @param ma Vector of moving-average coefficients
#' @param log Logical; if \code{TRUE} the function returns the log-probability; if \code{FALSE} the function returns the probability
#'
#' data(garma)
#' AR <- c(0.8, -0.2)
#' MA <- c(0.6,  0.3)
#'
#' #Compute the cumulative probability of the GARMA output
#' (PROBS <- pGARMA(SERIES, ar = AR, ma = MA))
#'
pGARMA <- function(x,
                   cond     = FALSE,
                   mean     = 0,
                   errorvar = 1,
                   ar       = numeric(0),
                   ma       = numeric(0),
                   log      = FALSE) {

  #Check inputs
  if (!is.numeric(x))            { stop('Error: Time series input should be numeric') }
  if (is.array(x)) {
    if (length(dim(x)) > 2)      { stop('Error: Time series input should be a vector or matrix, not a larger array') } }
  if (!is.numeric(mean))         { stop('Error: mean must be numeric') }
  if (length(mean) != 1)         { stop('Error: mean must be a single numeric value') }
  if (!is.numeric(errorvar))     { stop('Error: errorvar must be numeric') }
  if (length(errorvar) != 1)     { stop('Error: errorvar must be a single numeric value') }
  if (errorvar < 0)              { stop('Error: errorvar cannot be negative') }
  if (!is.logical(log))          { stop('Error: log must be a logical value') }
  if (length(log) != 1)          { stop('Error: log must be a single logical value') }

  #Set dimensions (converting vector to matrix if required)
  if (!is.array(x)) { x <- matrix(x, nrow = 1) }
  n <- nrow(x);
  m <- ncol(x);

  #Check input cond
  if (!is.vector(cond))          { stop('Error: cond should be FALSE or a logical vector') }
  if (length(cond) == 1)         {
    if (cond == FALSE)           { cond <- rep(FALSE, m); } }
  mm   <- length(cond);
  cc   <- sum(cond);
  if (mm != m)                   { stop('Error: cond must have same length as input time-series') }

  #Check if there are missing conditional values
  MISS <- sum(is.na(x[, cond]));
  if (MISS != 0)                 { stop('Error: Missing conditioning values --- there are NA values in the conditioning elements') }

  #Deal with trivial case where all elements are conditioning elements
  if (cc == m)                   {
    OUT <- rep(0, n);
    names(OUT) <- sprintf('Series[%s]', 1:n);
    warning('All values are set as conditional values --- output probability is set equal to one by convention');
    if (log) { return(OUT) } else { return(exp(OUT)) } }

  #Compute mean vector and variance matrix
  MEAN <- rep(mean, m);
  VAR  <- errorvar*ARMA.var(n = m, ar = ar, ma = ma);

  #Decompose mean vector and conditional variance matrix and compute conditional variance matrix
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
    CVAR   <- VAR11 - VAR12 %*% solve(VAR22, VAR21); }

  #Compute output probabilities
  OUT <- rep(NA, n);
  for (i in 1:n) {
    XX <- x[i, !cond];
    YY <- x[i,  cond];
    MM <- !is.na(XX);
    if (cc == 0) {
      CMEAN <- MEAN; } else {
      CMEAN <- MEAN1 + VAR12 %*% solve(VAR22, YY - MEAN2); }
    if (sum(MM) == 0) {
      OUT[i] <- 0;
      warning('All non-conditional elements in row ', i, ' of input are missing (marginal) values --- output probability is set equal to one by convention'); }
    if (sum(MM) == 1) {
      OUT[i] <- pnorm(q     = XX[MM],
                      mean  = CMEAN[MM],
                      sd    = CVAR[MM,MM],
                      log.p = TRUE); }
    if (sum(MM) >  1) {
      stopifnot("mvtnorm packaged required for multivariate features."=requireNamespace('mvtnorm', quietly=TRUE))
      OUT[i] <- mvtnorm::pmvnorm(lower  = rep(-Inf, sum(MM)),
                                 upper  = XX[MM],
                                 mean   = CMEAN[MM],
                                 sigma  = CVAR[MM, MM],
                                 abseps = 10^(-6),
                                 algorithm = mvtnorm::GenzBretz());
      OUT[i] <- log(OUT[i]); } }

  #Add labels
  names(OUT) <- sprintf('Series[%s]', 1:n);

  if (log) { OUT } else { exp(OUT) } }
