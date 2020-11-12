#' Permutation-spectrum test for time-series data
#'
#' Computes the permutation-spectrum test to detect a periodic signal in a real or complex time-series vector.  The null hypothesis for
#' the test is that the values in the time-series vector are independent and identically distributed with no signal, and the alternative
#' hypothesis is that the time-series has at least one remaining periodic signal.  The test statistic is the maximum scaled intensity
#' of the observed time-series vector.  The p-value for the test is the probability of observing a maximum scaled intensity at least as
#' large as the observed value under the null distribution of exchangeability.  The null distribution for the test is simulated by applying
#' random permutations to the observed time-series.  The number of simulations is controlled by the \code{sims} parameter.
#'
#' @param x A vector of time-series values (must have at least two data points)
#' @param sims Positive integer for the number of simulations to perform in the test
#' @param progress Logical; if \code{TRUE} the function uses a progress bar to track its simulations
#'
#' @examples
#'
#' data(garma)
#'
#' #Show the intensity of a time-series vector
#' spectrum.test(SERIES1, sims = 100)
spectrum.test <- function(x = NULL, sims = 10^6, progress = TRUE) {

  #Check data inputs
  if (!is.numeric(x)) {
  if (!is.complex(x))           { stop('Error: time-series x should be numeric or complex') } }
  if (!is.vector(x))            { stop('Error: time-series x should be a vector') }
  if (length(x) <= 1)           { stop('Error: time-series x must have at least two data points') }
  if (!is.numeric(sims))        { stop('Error: sims should be numeric') }
  if (length(sims) != 1)        { stop('Error: sims should be a single positive integer '); }
  if (as.integer(sims) != sims) { stop('Error: sims should be a positive integer'); }
  if (sims <= 0)                { stop('Error: sims should be a positive integer '); }
  if (sims <= 10^5)             { warning('low sims may give poor results for the spectrum test'); }
  if (!is.logical(progress))    { stop('Error: progress must be a logical value') }
  if (length(progress) != 1)    { stop('Error: progress must be a logical value'); }

  #Compute maximum intensity
  n   <- length(x);
  DFT <- fft(x - mean(x))/sqrt(n);
  SD  <- sqrt(sum(Mod(x - mean(x))^2/(n-1)));
  INT <- Mod(DFT)/SD;
  MAXINT <- max(INT);

  #Set description of test and data
  method    <- 'Permutation-Spectrum Test';
  x.type    <- ifelse(is.complex(x),
                      ifelse(!identical(Im(x), rep(0,n)), 'complex', 'real'), 'real');
  data.name <- paste0(x.type, ' time-series vector ', deparse(substitute(x)), ' with ', n, ' values');

  #Set null and alternative hypotheses
  null.value  <- NULL
  alternative <- "distribution of time-series vector is not exchangeable (at least one periodic signal is present)";

  #Determine maximum intensity of each noise vector
  MAXINT.SIM <- rep(0, sims);
  if (progress) { PROGRESS <- txtProgressBar(min = 0, max = sims, initial = 0, char = '='); }
  for (i in 1:sims) {
    S   <- sample(x);
    DFT <- fft(S - mean(S))/sqrt(n);
    INT <- Mod(DFT)/SD;
    MAXINT.SIM[i] <- max(INT);
    if (progress) { setTxtProgressBar(PROGRESS, value = i); } }
  if (progress) { close(PROGRESS); }

  #Calculate test statistics
  estimate    <- NULL;
  statistic   <- MAXINT;
  attr(statistic, "names") <- "maximum scaled intensity";

  #Calculate p-value
  p.value     <- sum(MAXINT.SIM >= MAXINT)/sims;
  attr(p.value, "names") <- NULL;

  #Create htest object
  TEST        <- list(method = method, data.name = data.name, x = x,
                      null.value = null.value, alternative = alternative,
                      sample.size = n, sims = sims, maxint.sim = MAXINT.SIM,
                      estimate = estimate, statistic = statistic, p.value = p.value);
  class(TEST) <- c("spectrum.test", "htest");
  TEST; }
