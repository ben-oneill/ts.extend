#' Covariance/correlation matrix for the stationary ARMA model
#'
#' This function computes the covariance/correlation matrix for a stationary auto-regressive moving-average (ARMA) model.  The user specifies
#' the matrix size ```n``` and the function returns a matrix of covariance/correlation values at all times ```Time[1], ... , Time[n]``` (in the
#' case where conditioning values are specified using the ```condvals``` argument, only the time values for non-conditional values are included).
#' The function requires the model to be stationary, which means that the vector of auto-regression coefficients must give an auto-regressive
#' characteristic polynomial with roots outside the unit circle.
#'
#' @param n Positive integer giving the number of values in the time-series (output variance matrix is an n x n matrix)
#' @param condvals Either a single value ```NA``` or a numeric vector with ```n``` elements; numeric entries are conditioning values for the generated vector
#' @param ar Vector of auto-regressive coefficients (all roots of AR characteristic polynomial must be outside the unit circle)
#' @param ma Vector of moving-average coefficients
#' @param corr Logical; if ```TRUE``` the function returns the correlation matrix; if ```FALSE``` the function returns the covariance matrix

ARMA.var <- function(n,
                     condvals = as.numeric(NA),
                     ar       = numeric(0),
                     ma       = numeric(0),
                     corr     = FALSE) {

  #Check inputs
  if (!is.numeric(n))            { stop('Error: n should be a positive integer') }
  if (length(n) != 1)            { stop('Error: n should be a single positive integer') }
  if (as.integer(n) != n)        { stop('Error: n should be a positive integer') }
  if (n <= 0)                    { stop('Error: n must be at least one') }
  if (!is.numeric(ar))           { stop('Error: ar should be a numeric vector') }
  if (!is.numeric(ma))           { stop('Error: ma should be a numeric vector') }
  if (!is.logical(corr))         { stop('Error: corr should be a logical value') }
  if (length(corr) != 1)         { stop('Error: corr should be a single logical value') }

  #Check input condvals
  if (!is.vector(condvals))      { stop('Error: condvals must be a conditioning vector') }
  if (length(condvals) == 1) {
    if (is.na(condvals))       {
      condvals <- as.numeric(rep(NA, n)); } }
  mm   <- length(condvals);
  if (mm != n)                   { stop('Error: condvals must have length n') }
  cond <- !is.na(condvals);
  cc   <- sum(cond);

  #Compute marginal variance matrix
  ACV <- ARMA.autocov(n = n, ar, ma, corr = corr);
  VAR <- matrix(0, nrow = n, ncol = n);
  for (i in 1:n) {
  for (j in 1:n) {
    VAR[i,j] <- ACV[abs(i-j)+1]; } }
  rownames(VAR) <- sprintf('Time[%s]', 1:n);
  colnames(VAR) <- sprintf('Time[%s]', 1:n);

  #Compute conditional variance matrix
  if (cc == 0) {
    CVAR   <- VAR; }
  if ((cc > 0) & (cc < n)) {
    VAR11  <- VAR[!cond, !cond, drop = FALSE];
    VAR12  <- VAR[!cond,  cond, drop = FALSE];
    VAR21  <- VAR[ cond, !cond, drop = FALSE];
    VAR22  <- VAR[ cond,  cond, drop = FALSE];
    CVAR   <- VAR11 - VAR12 %*% solve(VAR22, VAR21); }

  CVAR; }
