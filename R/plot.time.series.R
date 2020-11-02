#' Plot scatterplot matrix of time-series vectors
#'
#' This is a custom plot function that operates on objects of class \code{time.series} (which is the output generated from the \code{rGARMA}
#' function).  This plot function generates a scatterplot matrix of time-series vectors using either \code{ggplot} or \code{base} graphics.
#' The user must input either a single time-series vector \code{x} or a matrix \code{x} where each row is one time-series vector.  The function
#' generates a scatterplot matrix showing each of the time-series.  If the plot is generated using \code{ggplot} graphics then the user has an
#' option to include a background showing all points from all the series in grey.  The user may choose to print or assign the object, or both.
#' Since the function generates a scatterplot of time-series plots, there are certain limits in the output.  If the user attempts to generate
#' the plot for a time-series matrix with more than 36 time-series, the user will be prompted to continue and prompted about whether they want
#' to include the background.  The prompts can be removed in the arguments of the function.
#'
#' @param x A vector or matrix of time-series values (if a matrix, each time-series should be one row of the matrix)
#' @param ggplot Logical; if \code{TRUE} the scatterplot is a \code{ggplot} object; if \code{FALSE} it is a \code{base} plot object
#' @param background Logical; if \code{TRUE} then each \code{ggplot} scatterplot will include background points for all the time-series vectors
#' @param print Logical; if \code{TRUE} the scatterplot is printed
#' @param user.prompt Logical; if \code{TRUE} the user will be prompted for choices when the number if time-series is large
#' @param ...   unused
#'
#' #@examples
#'
#' data(garma)
#' plot(SERIES)
#'
plot.time.series <- function(x, ggplot = TRUE, background = TRUE, print = TRUE, user.prompt = TRUE, ...) {

  #Check inputs
  if (!is.numeric(x))           { stop('Error: Time series input should be numeric') }
  if (is.array(x)) {
    if (length(dim(x)) > 2)     { stop('Error: Time series input should be a vector or matrix, not a larger array') } }
  if (!is.logical(ggplot))      { stop('Error: ggplot should be a logical value') }
  if (length(ggplot) != 1)      { stop('Error: ggplot should be a single logical value') }
  if (!is.logical(background))  { stop('Error: background should be a logical value') }
  if (length(background) != 1)  { stop('Error: background should be a single logical value') }
  if (!is.logical(print))       { stop('Error: print should be a logical value') }
  if (length(print) != 1)       { stop('Error: print should be a single logical value') }
  if (!is.logical(user.prompt)) { stop('Error: user.prompt should be a logical value') }
  if (length(user.prompt) != 1) { stop('Error: user.prompt should be a single logical value') }

  #Set dimensions (converting vector to matrix if required)
  if (!is.array(x)) { x <- matrix(x, nrow = 1) }
  n <- nrow(x);
  m <- ncol(x);
  B <- background;

  #Check number of plots and prompt user if needed
  #Applies if there are more than 36 plots (6 x 6 plot matrix)
  if (n > 36) {
    if (user.prompt) {
      MSG  <- paste0('You propose to create a scatterplot of ', n, ' plots - This may be slow and it may look bad \nAre you sure you want to continue? Y/N ');
      INV  <- paste0('Invalid response --- Are you sure you want to continue? Y/N ');
      CONT <- readline(MSG);
      while (!(CONT %in% c('Y','N'))) { CONT <- readline(INV); } }
    if (CONT == 'N') { stop('Error: Too many time-series --- plot aborted') }
    if (all(c(ggplot, background, user.prompt))) {
      MSG  <- paste0('We recommend removing the background \nDo you want to remove the background? Y/N ');
      INV  <- paste0('Invalid response --- Do you want to remove the background? Y/N ');
      BACK <- readline(MSG);
      while (!(BACK %in% c('Y','N'))) { BACK <- readline(INV); }
      B    <- (BACK == 'N'); } }

  #If ggplot2 is selected and installed we create the plot in this package
  if (ggplot) {
    if (requireNamespace('ggplot2', quietly = TRUE)) {

      #Create data-frame for plotting
      PLOTDATA <- expand.grid(1:m, 1:n)[,2:1];
      colnames(PLOTDATA) = c('n', 'Time');
      PLOTDATA$Value <- 0;
      for (i in 1:nrow(PLOTDATA)) {
        PLOTDATA$Value[i] <- x[PLOTDATA$n[i], PLOTDATA$Time[i]] }
      PLOTDATA <- as.data.frame(PLOTDATA);

      #Create labels for plots
      if (is.null(rownames(x))) { LABELS <- 1:n; } else { LABELS <- rownames(x); }
      names(LABELS) <- 1:n;

      #Create plot in ggplot2
      if (B) {
        PLOT <- ggplot2::ggplot(ggplot2::aes_string(x = "Time", y = "Value"), data = PLOTDATA) +
          ggplot2::geom_point(ggplot2::aes_string(x = "Time", y = "Value"), data = PLOTDATA[,2:3],
                              size = 1, colour = 'grey') +
          ggplot2::geom_line(size = 1, colour = 'blue') +
          ggplot2::facet_wrap(~ n, labeller = ggplot2::as_labeller(LABELS)) +
          ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold')) +
          ggplot2::ggtitle('Time Series Plots'); } else {
            PLOT <- ggplot2::ggplot(ggplot2::aes_string(x = "Time", y = "Value"), data = PLOTDATA) +
              ggplot2::geom_line(size = 1, colour = 'blue') +
              ggplot2::facet_wrap(~ n, labeller = ggplot2::as_labeller(LABELS)) +
              ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                             plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold')) +
              ggplot2::ggtitle('Time Series Plots'); } } else {
                warning('Package ggplot2 is not installed - switching to base graphics instead'); } }

  #If ggplot2 is not selected or not installed we create the plot using base graphics
  if ((!ggplot)|(!requireNamespace('ggplot2', quietly = TRUE))) {

    #Restore par settings on exit
    OLDPAR <- par('mfrow', 'mar', 'mgp', 'las');
    on.exit(par(OLDPAR));

    #Enable the device
    dev.new();
    dev.control('enable');

    #Create layout for plot
    COLS   <- max(1, ceiling(sqrt(n)));
    ROWS   <- ceiling(n/COLS);
    LAYOUT <- matrix(c(rep(1, COLS), 2:(ROWS*COLS+1)), nrow = ROWS+1, ncol = COLS, byrow = TRUE);
    layout(LAYOUT, heights = c(0.3, rep(1, ROWS)));

    #Create and save plot matrix
    par(mar = c(0.5, 2.1, 0.5, 2.1), mgp = c(3, 1, 0), las = 0);
    plot.new();
    text(0.5,0.5, 'Time Series Plots', cex = 2, font = 2);
    par(mar = c(5.1, 4.1, 2.1, 2.1), mgp = c(3, 1, 0), las = 0);
    YMIN <- min(x);
    YMAX <- max(x);
    for (i in 1:n) { plot(x[i,], type = 'l',
                          ylim = c(YMIN, YMAX),
                          xlab = 'Time', ylab = 'Value') }

    #Record plot and turn off device
    PLOT <- recordPlot();
    dev.off(); }

  if (print) { print(PLOT); }

  PLOT; }
