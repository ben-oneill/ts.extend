#' Auto-covariance/auto-correlation function for the stationary ARMA model
#'
#' This function computes a vector of output values from the auto-covariance/auto-correlation function for a stationary auto-regressive
#' moving-average (ARMA) model.  The user specifies the vector size \code{n} and the function returns a vector of auto-covariance/
#' auto-correlation values at all lags \code{Lag[0], ... , Lag[n-1]}.  The function requires the model to be stationary, which means that
#' the vector of auto-regression coefficients must give an auto-regressive characteristic polynomial with roots outside the unit circle.
#'
#' @param n Positive integer giving the number of consecutive values in the time-series (output is a vector of length \code{n})
#' @param ar Vector of auto-regressive coefficients (all roots of AR characteristic polynomial must be outside the unit circle)
#' @param ma Vector of moving-average coefficients
#' @param corr Logical; if \code{TRUE} the function returns the auto-correlation function; if \code{FALSE} the function returns the auto-covariance function
#'
#' @examples
#'
#' data(garma)
#'
#' AR <- c(0.8, -0.2)
#' MA <- c(0.6,  0.3)
#' #Compute the auto-correlation function
#' ARMA.autocov(n = 6, ar = AR, ma = MA, corr = TRUE)

ARMA.autocov <- function(n,
                         ar      = numeric(0),
                         ma      = numeric(0),
                         corr    = FALSE) {

  #Check inputs
  if (!is.numeric(n))            { stop('Error: n should be a positive integer') }
  if (length(n) != 1)            { stop('Error: n should be a single positive integer') }
  if (as.integer(n) != n)        { stop('Error: n should be a positive integer') }
  if (n <= 0)                    { stop('Error: n must be at least one') }
  if (!is.numeric(ar))           { stop('Error: ar should be a numeric vector') }
  if (!is.numeric(ma))           { stop('Error: ma should be a numeric vector') }
  if (!is.logical(corr))         { stop('Error: corr should be a logical value') }
  if (length(corr) != 1)         { stop('Error: corr should be a single logical value') }

  #Compute model orders
  p <- length(ar);
  q <- length(ma);
  r <- max(p, q+1);

  #Compute the autocorrelation function and leading variance
  #Case where there is no AR or MA
  if ((p == 0) && (q == 0)) {
    ACF <- c(1, rep(0, n-1));
    VAR <- ACF[1]; }

  #Compute the autocorrelation function and leading variance
  #Case where there is no AR but there is MA
  if (p == 0) {
    if (q > 0)  {
      x      <- c(1, ma);
      FILTER <- stats::filter(c(x, rep(0, q)), rev(x), sides = 1);
      ACF    <- FILTER[-(1L:q)];
      if (n > q+1) { ACF <- c(ACF, rep(0, n-q-1)); }
      VAR <- ACF[1L]
      ACF <- ACF/VAR; } }

  #Compute the autocorrelation function and leading variance
  #Case where there is AR and MA
  if (p > 0) {
    if (r > 1) {

      #Adjust ar vector and order p
      if (r > p) { ar <- c(ar, rep(0, r-p));
                   p  <- r; }

      #Construct linear equations for ACF (A %*% ACF = b)
      AA <- c(1, -ar);
      A  <- matrix(0, nrow = p+1L, ncol = 2*p+1L);
      for (i in 1L:nrow(A)) {
      for (k in 0L:p)       {
          A[i,i+k] <- AA[k+1]; } }
      A[, 1L:p] <- A[, 1L:p] + A[, (2*p+1):(p+2)];
      A <- A[rev(1:(p+1)), rev(1:(p+1))];
      b <- c(1, rep(0, p));

      if (q > 0) {
        psi   <- c(1, ARMAtoMA(ar, ma, q));
        theta <- c(1, ma, rep(0, q+1));
        for (k in 1:(q+1)) { b[k] <- sum(psi*theta[k:(q+k)]); } }

      #Compute ACF and VAR
      ACF <- solve(A, b);
      VAR <- ACF[1];
      ACF <- ACF[-1L]/VAR; } else {
        VAR <- 1/(1-ar^2);
        ACF <- ar; }

    #Extend to required length
    if (n > p+1) {
      FILTER <- stats::filter(rep(0, n-p-1), ar, "recursive", init = rev(ACF));
      ACF    <- c(ACF, FILTER); }
    ACF <- c(1, ACF[1L:(n-1)]); }

  #Add names and output
  names(ACF) <- sprintf('Lag[%s]', 0:(n-1));
  if (corr) { ACF; } else { VAR*ACF; } }
