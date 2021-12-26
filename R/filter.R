#' Filter a time-series vector/matrix based on its spectral intensity
#'
#' This function filters a time-series based on its spectral intensity.  The user inputs the time-series vector or matrix and 
#' specifies whether it is to be centered and/or scaled.  Centering subtracts the sample mean of the vector prior to conversion 
#' into frequency space; this sets the intensity of the signal to zero at the zero frequency.  Scaling scales the intensity so 
#' the that norm of the intensity vector is equal to the number of values in the time-series.  The function then applies the 
#' specified highpass/lowpass filter and returns the filtered signals.  By default, both the highpass and lowpass filters remove
#' parts of the signal that do not meet the filter requirements; if preferred the filters can attenuate (rather than removing) 
#' these signals.
#'
#' @param x A vector of time-series values
#' @param centered Logical; if \code{TRUE} the time-series vector is centred (mean removed) before computing its intensity
#' @param centred An alternative name for the \code{centered} input
#' @param scaled Logical; if \code{TRUE} the intensity measure is scaled so that its norm is equal to the number of values in the time-series
#' @param highpass The minimum intensity value for the highpass filter (or NULL if there is no highpass filter)
#' @param highpass.remove Logical; if \code{TRUE} then the highpass filter removes (instead of attenuating) the low-frequency signal
#' @param lowpass The maximum intensity value for the lowpass filter (or NULL if there is no lowpass filter)
#' @param lowpass.remove Logical; if \code{TRUE} then the lowpass filter removes (instead of attenuating) the high-frequency signal

filter <- function(x, centered = TRUE, centred = centered, scaled = TRUE, 
                   highpass = NULL, highpass.remove = TRUE, 
                   lowpass  = NULL, lowpass.remove  = TRUE) {

  #Check input x
  if (!is.numeric(x))        {
    if (!is.complex(x))        { stop('Error: x must either be numeric or complex') } }
  if (is.array(x)) {
    if (length(dim(x)) > 2)     { stop('Error: Time series input should be a vector or matrix, not a larger array') } }
  if (!is.array(x)) { x <- matrix(x, nrow = 1) }
  n <- nrow(x)
  m <- ncol(x)
  if (m == 0)        { stop('Error: x must have at least one value') }
  if (m == 1)        {
    if (scaled)                { stop('Error: x must have at least two values for scaled intensity') } }
  
  #Check inputs centred/centered and scaled
  if (!is.logical(centred))  { stop('Error: centered/centred must be a logical value') }
  if (length(centred) != 1)  { stop('Error: centered/centred must be a single logical value') }
  if (!is.logical(scaled))   { stop('Error: scaled must be a logical value') }
  if (length(scaled) != 1)   { stop('Error: scaled must be a single logical value') }
  
  #Check inputs highpass and lowpass
  if (!is.null(highpass)) {
    if (!is.numeric(highpass)) stop('Error: highpass must be a positive number (or NULL)')
    if (length(highpass) != 1) stop('Error: highpass must be a single number (or NULL)')
    if (highpass <= 0)         stop('Error: highpass must be a positive number (or NULL)') }
  if (!is.null(lowpass)) {
    if (!is.numeric(lowpass))  stop('Error: lowpass must be a positive number (or NULL)')
    if (length(lowpass) != 1)  stop('Error: lowpass must be a single number (or NULL)')
    if (lowpass <= 0)          stop('Error: lowpass must be a positive number (or NULL)') }
  if ((!is.null(highpass))&(!is.null(lowpass))) {
    if (lowpass < highpass)    stop('You have lowpass < highpass --- this means that all data is filtered out') }
  if ((is.null(highpass))&(is.null(lowpass))) {
    warning('You have not specified highpass or lowpass --- no filter applied')
    return(x) }
 
  #Compute DFT and intensity
  DFT <- matrix(complex(1), nrow = n, ncol = m)
  rownames(DFT) <- rownames(x)
  colnames(DFT) <- sprintf('Freq[%s/%s]', ((1:m) - 1), m)
  INT <- matrix(0, nrow = n, ncol = m)
  rownames(INT) <- rownames(x)
  colnames(INT) <- sprintf('Freq[%s/%s]', ((1:m) - 1), m)
  MEAN <- rep(0, n)
  SD   <- rep(0, n)
  for (i in 1:n) {
    MEAN[i] <- mean(x[i,])
    SD[i]   <- sqrt(sum((Mod(x[i,] - MEAN[i])^2)/(m-1)))
    DFT[i,] <- fft(x[i,] - centred*MEAN[i])/sqrt(m)
    INT[i,] <- Mod(DFT[i,])
    if (centred) { 
      DFT[i, 1] <- 0
      INT[i, 1] <- 0 } else { 
      SD[i] <- sqrt((m-1)/m)*SD[i] }
    if (scaled)  { 
      DFT[i,] <- DFT[i,]/SD[i]
      INT[i,] <- INT[i,]/SD[i] } }
  
  #Filter the signal using highpass/lowpass filter
  DFT.FILT <- DFT
  INT.FILT <- INT
  x.FILT   <- matrix(is.complex(0), nrow = n, ncol = m)
  for (i in 1:n) {
    if (!is.null(lowpass))  { 
      for (k in 1:m) {
        if (INT[i, k] > lowpass)  { 
          DFT.FILT[i, k] <- ifelse(lowpass.remove,  0, (DFT[i, k]/INT[i, k])*lowpass)
          INT.FILT[i, k] <- Mod(DFT.FILT[i, k]) } } }
    if (!is.null(highpass)) { 
      for (k in 1:m) {
        if (INT[i, k] < highpass) { 
          DFT.FILT[i, k] <- ifelse(highpass.remove, 0, (DFT[i, k]/INT[i, k])*highpass)
          INT.FILT[i, k] <- Mod(DFT.FILT[i, k]) } } }
    x.FILT[i, ]   <- fft(DFT.FILT[i, ], inverse = TRUE)/sqrt(m)
    if (scaled)  { x.FILT[i, ] <- x.FILT[i, ]*SD[i] }
    if (centred) { x.FILT[i, ] <- MEAN[i] + x.FILT[i, ] } }
  if (n == 1) { x.FILT <- as.vector(x.FILT) }
  if (is.numeric(x)) { x.FILT <- Re(x.FILT) }

  #Output the filtered vector/matrix
  x.FILT }
