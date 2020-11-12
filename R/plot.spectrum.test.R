#' Plot of the Permutation-Spectrum Test
#'
#' This function generates a dual plot showing the results of the permutation-spectrum testusing either \code{ggplot} or \code{base} graphics.
#' The user must input a test object produced by the \code{spectrum.test} function.  The function produces a dual plot showing the scaled intensity
#' of the time-series vector and the simulated null distribution of the maximum scaled intensity under the null hypothesis of an IID vector.
#' The plots also report the value of the maximum scaled intensity and the resulting p-value for the test.  This dual plot forms a useful
#' companion to the permutation-spectrum test; it allows the user to visualise the simulated null distribution and test statistic.
#'
#' @param x A \code{spectrum.test} object produced by the \code{spectrum.test} function
#' @param ggplot Logical; if \code{TRUE} the scatterplot is a \code{ggplot} object; if \code{FALSE} it is a \code{base} plot object
#' @param print Logical; if \code{TRUE} the scatterplot is printed
#' @param ...   unused
#' @examples
#'
#' data(garma)
#'
#' #Show the intensity of a time-series vector
#' TEST <- spectrum.test(SERIES1, sims = 100)
#' plot(TEST)

plot.spectrum.test <- function(x, ggplot = TRUE, print = TRUE, ...) {

  test <- x

  #Check test input
  if (is.null(x[["x"]]))       { stop('Error: test object should contain a time-series vector x') } else {
    x <- x[["x"]] }
  if (is.null(test[["maxint.sim"]])) { stop('Error: test object should contain simulated intensity values in maxint.sim') } else {
    maxint.sim <- test[["maxint.sim"]] }

  #Check other inputs
  if (!is.numeric(x)) {
    if (!is.complex(x))           { stop('Error: time-series x should be numeric or complex') } }
  if (!is.vector(x))            { stop('Error: time series x should be a vector') }
  if (!is.numeric(maxint.sim))  { stop('Error: maxint.sim input should be numeric') }
  if (!is.logical(ggplot))      { stop('Error: ggplot should be a logical value') }
  if (length(ggplot) != 1)      { stop('Error: ggplot should be a single logical value') }
  if (!is.logical(print))       { stop('Error: print should be a logical value') }
  if (length(print) != 1)       { stop('Error: print should be a single logical value') }

  #Find Intensity of time-series vector
  sims <- length(maxint.sim);
  n    <- length(x);
  DFT  <- fft(x - mean(x))/sqrt(n);
  SD   <- sqrt(sum(Mod(x - mean(x))^2/(n-1)));
  INT  <- Mod(DFT)/SD;
  FREQ <- (0:(n-1))/n;

  #Adjust INT and FREQ for real vectors
  x.type <- ifelse(is.complex(x),
                   ifelse(!identical(Im(x), rep(0,n)), 'complex', 'real'), 'real');
  if (x.type == 'real') { nn <- ceiling((n+1)/2) } else { nn <- n }
  INT  <- INT[1:nn];
  FREQ <- FREQ[1:nn];

  #Find maximum intensity and frequency
  MAXFREQ <- FREQ[which.max(INT)][1];
  MAXINT  <- max(INT);
  MAXY    <- max(c(maxint.sim, MAXINT));
  p.value <- sum(maxint.sim >= MAXINT)/sims;

  #Set plot title and label
  TITLE <- paste0('Permutation Spectrum Plot \n (Maximum Scaled Intensity = ', round(MAXINT, 4),
                  ', p-value = ', sprintf('%.4f', p.value), ')');
  LABEL <- paste0('Simulated Density (',
                  formatC(sims, format = 'f', big.mark = ',', digits = 0), ' simulations)');

  #If ggplot2 is selected and installed we create the plot in this package
  if (ggplot) {
    if (requireNamespace('ggplot2',   quietly = TRUE)) {
      if (requireNamespace('grid',      quietly = TRUE)) {
        if (requireNamespace('gridExtra', quietly = TRUE)) {

          #Set theme
          THEME <- ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                                  plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold'),
                                  axis.title.y  = ggplot2::element_text(margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 0)));

          #Generate plots
          FIGURE1 <- ggplot2::ggplot(ggplot2::aes_string(x = "Frequency", y = "Intensity"),
                                     data = data.frame(Frequency = FREQ, Intensity = INT)) +
            ggplot2::geom_bar(stat = 'identity', fill = 'DarkRed') +
            ggplot2::geom_point(ggplot2::aes_string(x = "Frequency", y = "Intensity"),
                                data = data.frame(Frequency = MAXFREQ, Intensity = MAXINT),
                                colour = 'Red', size = 4) +
            ggplot2::expand_limits(y = c(0, MAXY)) + THEME +
            ggplot2::ylab('Scaled Intensity');
          FIGURE2 <- ggplot2::ggplot(ggplot2::aes_(x = ~0, y =~MaxIntensity),
                                     data = data.frame(MaxIntensity = maxint.sim)) +
            ggplot2::geom_violin(fill = 'orange', draw_quantiles = c(0.25, 0.5, 0.75)) +
            ggplot2::annotate('point', x = 0, y = MAXINT, colour = 'red', size = 4) +
            ggplot2::expand_limits(y = c(0, MAXY)) + THEME +
            ggplot2::theme(axis.text.x = ggplot2::element_text(colour = 'white')) +
            ggplot2::scale_y_continuous(position = 'right') +
            ggplot2::xlab(LABEL) + ggplot2::ylab('');
          TOP     <- grid::textGrob(TITLE, gp = grid::gpar(fontface = 'bold'));
          PLOT    <- gridExtra::grid.arrange(FIGURE1, FIGURE2, nrow = 1, ncol = 2, top = TOP); } else {

            warning('Package gridExtra is not installed - switching to base graphics instead');  } } else {
              warning('Package grid is not installed - switching to base graphics instead');       } } else {
                warning('Package ggplot2 is not installed - switching to base graphics instead');    } }

  #If ggplot2 is not selected or not installed we create the plot using base graphics
  if ((!ggplot)|(!requireNamespace('ggplot2', quietly = TRUE))) {

    #Restore par settings on exit
    OLDPAR <- par('mfrow', 'mar', 'mgp', 'las');
    on.exit(par(OLDPAR));

    #Enable the device
    dev.new();
    dev.control('enable');

    #Create layout for plot
    LAYOUT <- matrix(c(rep(1, 2), 2:3), nrow = 2, ncol = 2, byrow = TRUE);
    layout(LAYOUT, heights = c(0.1, 1));

    #Create and save plot matrix
    par(mar = c(0.5, 2.1, 0.5, 2.1), mgp = c(3, 1, 0), las = 0);
    plot.new();
    text(0.5,0.5, TITLE, cex = 1.2, font = 2);
    par(mar = c(5.1, 4.1, 2.1, 2.1), mgp = c(3, 1, 0), las = 0);

    #Generate plots
    DENS <- density(maxint.sim);
    YMAX <- max(DENS$x);
    barplot(INT, names.arg = FREQ, ylim = c(0, YMAX), xaxt = 'n', xlab = 'Frequency', ylab = 'Scaled Intensity');
    axis(1, at = 0.2*(0:5)*length(INT)*(6/5), labels = 0.1*(0:5));
    points(x = which.max(INT)*(6/5), y = MAXINT, col = 'red', pch = 19, cex = 2);
    matplot(x = cbind(-DENS$y, DENS$y), y = DENS$x, xaxt = 'n', yaxs = 'i',
            type = c('l', 'l'), lty = c(1, 1),
            col = c('black', 'black'),
            xlim = c(-max(DENS$y), max(DENS$y)),
            ylim = c(0, YMAX),
            xlab = LABEL, ylab = '')
    polygon(x =  DENS$y, y = DENS$x, col = 'orange', border = FALSE);
    polygon(x = -DENS$y, y = DENS$x, col = 'orange', border = FALSE);
    matplot(x = cbind(-DENS$y, DENS$y), y = DENS$x, yaxs = 'i',
            type = c('l', 'l'), lty = c(1, 1),
            col = c('black', 'black'), add = TRUE)
    points(x = 0, y = MAXINT, col = 'red', pch = 19, cex = 2);
    mtext(TITLE, outer = TRUE, cex = 1.5);

    #Record plot and turn off device
    PLOT <- recordPlot();
    dev.off(); }

  if (print) { print(PLOT); }

  PLOT; }
