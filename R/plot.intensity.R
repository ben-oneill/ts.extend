#' Plot scatterplot matrix of intensity vectors
#'
#' This is a custom plot function that operates on objects of class ```intensity``` (which is the output generated from the ```intensity```
#' function).  This plot function generates a scatterplot matrix of intensity vectors using either ```ggplot``` or ```base``` graphics.
#' The user must input either a single intensity vector ```x``` or a matrix ```x``` where each row is one intensity vector.  The function
#' generates a scatterplot matrix showing each of the intensity barplots.  The user may choose to print or assign the object, or both.
#' Since the function generates a scatterplot of intensity plots, there are certain limits in the output.  If the user attempts to generate
#' the plot for a time-series matrix with more than 36 intensity vectors, the user will be prompted to continue.  The prompts can be removed
#' in the arguments of the function.
#'
#' @param x A vector or matrix of intensity values (if a matrix, each intensity vector should be one row of the matrix)
#' @param ggplot Logical; if ```TRUE``` the scatterplot is a ```ggplot``` object; if ```FALSE``` it is a ```base``` plot object
#' @param print Logical; if ```TRUE``` the scatterplot is printed
#' @param user.prompt Logical; if ```TRUE``` the user will be prompted for choices when the number if intensity vectors is large
#' @param ...   unused
#'
#' @examples
#'
#' data(garma)
#' INT     <- intensity(SERIES1)
#' plot(INT)
plot.intensity <- function(x, ggplot = TRUE, print = TRUE, user.prompt = TRUE, ...) {

  #Check inputs
  if (!is.numeric(x))           { stop('Error: Time series input should be numeric') }
  if (is.array(x)) {
    if (length(dim(x)) > 2)     { stop('Error: Time series input should be a vector or matrix, not a larger array') } }
  if (!is.logical(ggplot))      { stop('Error: ggplot should be a logical value') }
  if (length(ggplot) != 1)      { stop('Error: ggplot should be a single logical value') }
  if (!is.logical(print))       { stop('Error: print should be a logical value') }
  if (length(print) != 1)       { stop('Error: print should be a single logical value') }
  if (!is.logical(user.prompt)) { stop('Error: user.prompt should be a logical value') }
  if (length(user.prompt) != 1) { stop('Error: user.prompt should be a single logical value') }

  #Set dimensions (converting vector to matrix if required)
  ATTR <- attributes(x);
  if (!is.array(x)) { x <- matrix(x, nrow = 1); }
  n <- nrow(x);
  m <- ncol(x);
  LENGTH <- ifelse('length' %in% names(ATTR), ATTR$length, m);
  SCALED <- ifelse('scaled' %in% names(ATTR), ATTR$scaled, FALSE);
  YLAB   <- ifelse(SCALED, 'Scaled Intensity', 'Intensity')

  #Check number of plots and prompt user if needed
  #Applies if there are more than 36 plots (6 x 6 plot matrix)
  if (n > 36) {
    if (user.prompt) {
      MSG  <- paste0('You propose to create a scatterplot of ', n, ' plots - This may be slow and it may look bad \nAre you sure you want to continue? Y/N ');
      INV  <- paste0('Invalid response --- Are you sure you want to continue? Y/N ');
      CONT <- readline(MSG);
      while (!(CONT %in% c('Y','N'))) { CONT <- readline(INV); } }
    if (CONT == 'N') { stop('Error: Too many time-series --- plot aborted') } }

  #If ggplot2 is selected and installed we create the plot in this package
  if (ggplot) {
    if (requireNamespace('ggplot2', quietly = TRUE)) {

      #Create data-frame for plotting
      PLOTDATA <- expand.grid(1:m, 1:n)[,2:1];
      colnames(PLOTDATA) = c('n', 'm');
      PLOTDATA$ Intensity <- 0;
      for (i in 1:nrow(PLOTDATA)) {
        PLOTDATA$Frequency[i] <- (PLOTDATA$m[i]-1)/LENGTH;
        PLOTDATA$Intensity[i] <- x[PLOTDATA$n[i], PLOTDATA$m[i]]; }
      PLOTDATA <- as.data.frame(PLOTDATA);

      #Create labels for plots
      if (is.null(rownames(x))) { LABELS <- 1:n; } else { LABELS <- rownames(x); }
      names(LABELS) <- 1:n;

      #Create plot in ggplot2
      if (n == 1) {
        PLOT <- ggplot2::ggplot(ggplot2::aes_string(x = "Frequency", y = "Intensity"), data = PLOTDATA) +
          ggplot2::geom_bar(stat = 'identity', fill = 'blue') +
          ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                         plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold')) +
          ggplot2::ggtitle('Intensity Plot') +
          ggplot2::ylab(YLAB); } else {
            PLOT <- ggplot2::ggplot(ggplot2::aes_string(x = "Frequency", y = "Intensity"), data = PLOTDATA) +
              ggplot2::geom_bar(stat = 'identity', fill = 'blue') +
              ggplot2::facet_wrap(~ n, labeller = ggplot2::as_labeller(LABELS)) +
              ggplot2::theme(plot.title    = ggplot2::element_text(hjust = 0.5, size = 14, face = 'bold'),
                             plot.subtitle = ggplot2::element_text(hjust = 0.5, face = 'bold')) +
              ggplot2::ggtitle('Intensity Plots') +
              ggplot2::ylab(YLAB); } } else {
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
    text(0.5,0.5, ifelse(n == 1, 'Intensity Plot', 'Intensity Plots'), cex = 2, font = 2);
    par(mar = c(5.1, 4.1, 2.1, 2.1), mgp = c(3, 1, 0), las = 0);
    YMAX <- max(x);
    for (i in 1:n) {
      barplot(x[i,], col = 'blue',
              ylim = c(0, YMAX), xaxt = 'n',
              xlab = 'Frequency', ylab = YLAB);
      if (m < LENGTH) {
        axis(1, at = 0.2*(0:5)*length(x[i,])*(6/5), labels = 0.1*(0:5)); } else {
          axis(1, at = 0.2*(0:5)*length(x[i,])*(6/5), labels = 0.2*(0:5)); } }

    #Record plot and turn off device
    PLOT <- recordPlot();
    dev.off(); }

  if (print) { print(PLOT); }

  PLOT; }

