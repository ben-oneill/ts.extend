#' Compute the spectral intensity of a time-series vector/matrix
#'
#' This function computes the spectral intensity of a time series vector/matrix.  The user inputs the time-series vector or matrix and specifies
#' whether it is to be centered and/or scaled.  Centering subtracts the sample mean of the vector prior to conversion into frequency space; this
#' sets the intensity of the signal to zero at the zero frequency.  Scaling scales the intensity so the that norm of the intensity vector is equal
#' to the number of values in the time-series.
#'
#' @param x A vector of time-series values
#' @param centered Logical; if \code{TRUE} the time-series vector is centred (mean removed) before computing its intensity
#' @param centred An alternative name for the \code{centered} input
#' @param scaled Logical; if \code{TRUE} the intensity measure is scaled so that its norm is equal to the number of values in the time-series
#' @param nyquist Logical; if \code{TRUE} the intensity vector is reduces for real time-series vectors to limit it to frequencies in the Nyquist range
#'
#' @examples
#'
#' data(garma)
#' intensity(SERIES1)

intensity <- function(x, centered = TRUE, centred = centered, scaled = TRUE, nyquist = TRUE) {

  #Check inputs
  if (!is.numeric(x))        {
    if (!is.complex(x))        { stop('Error: x must either be numeric or complex') } }
  if (is.array(x)) {
    if (length(dim(x)) > 2)     { stop('Error: Time series input should be a vector or matrix, not a larger array') } }
  if (length(x) == 0)        { stop('Error: x must have at least one value') }
  if (length(x) == 1)        {
    if (scaled)                { stop('Error: x must have at least two values for scaled intensity') } }
  if (!is.logical(centred))  { stop('Error: centered/centred must be a logical value') }
  if (length(centred) != 1)  { stop('Error: centered/centred must be a single logical value') }
  if (!is.logical(scaled))   { stop('Error: scaled must be a logical value') }
  if (length(scaled) != 1)   { stop('Error: scaled must be a single logical value') }
  if (!is.logical(nyquist))  { stop('Error: nyquist must be a logical value') }
  if (length(nyquist) != 1)  { stop('Error: nyquist must be a single logical value') }

  #Set dimensions (converting vector to matrix if required)
  if (!is.array(x)) { x <- matrix(x, nrow = 1) }
  n <- nrow(x);
  m <- ncol(x);

  #Compute intensity
  INT <- matrix(0, nrow = n, ncol = m);
  rownames(INT) <- rownames(x);
  colnames(INT) <- sprintf('Freq[%s/%s]', ((1:m) - 1), m);
  for (i in 1:n) {
    DFT <- fft(x[i,] - centred*mean(x[i,]))/sqrt(m);
    SD  <- sqrt(sum(Mod(x[i,] - mean(x[i,]))^2/(m-1)));
    INT[i,] <- Mod(DFT);
    if (centred) { INT[i, 1] <- 0; } else { SD <- sqrt((m-1)/m)*SD; }
    if (scaled)  { INT[i,] <- INT[i,]/SD; } }

  #Adjust INT for real vectors
  x.type <- ifelse(is.complex(x),
                   ifelse(!identical(Im(x), rep(0,m)), 'complex', 'real'), 'real');
  if (nyquist) {
    if (x.type == 'real') {
      mm  <- ceiling((m+1)/2);
      INT <- INT[, 1:mm]; } }

  #Convert back to vector if n = 1
  class(INT) <- c('matrix', 'intensity');
  if (n == 1) {
    INT <- as.vector(INT);
    names(INT) <- sprintf('Freq[%s/%s]', ((1:mm) - 1), m);
    class(INT) <- 'intensity'; }

  #Add class and attributes
  attr(INT, 'scaled') <- scaled;
  attr(INT, 'length') <- m;

  #Output the intensity
  INT; }
