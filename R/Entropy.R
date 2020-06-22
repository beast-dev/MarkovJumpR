#' Compute entropy of MarkovJump paths at specific times
#'
#' @param paths  MarkovJumps paths
#' @param time   Vector of times at which to evaluate
#'
#' @importFrom entropy entropy
#'
#' @export
getEntropy <- function(paths, time) {

  execute <- function(paths, time) {
    weights <- paths %>% filter(.data$startTime >= time, .data$endTime <= time) %>%
      group_by(.data$location) %>% summarise(weight = n())
    return(entropy(weights$weight))
  }

  Vectorize(execute, vectorize.args = "time")(paths = paths, time)
}

#' Plots entropy for one taxon
#'
#' @param entropy  Vector of entropy values
#' @param time     Vector of times
#' @param minTime   Minimum time to plot
#' @param maxTime   Maximum time to plot
#' @param drawBox  Add box around plot
#' @param xNoAxis  Do not show x-axis
#' @param xAt Time axis \code{at} values
#' @param xLabels Time axis \code{labels} values
#' @param fillColor Color to fill in below entropy time-series
#' @param lineColor Color of entropy time-series
#' @param ylab Y-axis label text
#' @param cex.label Y-axis label size
#'
#' @return NULL
#'
#' @importFrom graphics lines polygon title
#'
#' @export
plotEntropy <- function(entropy,
                      time,
                      minTime = min(time),
                      maxTime = max(time),
                      drawBox = FALSE,
                      xNoAxis = TRUE,
                      xAt = NULL,
                      xLabels = TRUE,
                      fillColor = NULL,
                      lineColor = "black",
                      ylab = "Entropy",
                      cex.label = 0.75) {

  plot(time, entropy, xlim = c(minTime, maxTime),
       # ylim = c(0, max(entropy)),
       type = "n",
       ylab = "",
       xlab = "",
       axes = FALSE)

  lines(time, entropy,  col = lineColor)

  if (ylab != "") {
    title(ylab = ylab, line = 0, cex.lab = cex.label)
  }

  if (!xNoAxis) {
    axis(1, at = xAt, labels = xLabels)
  }

  if (!is.null(fillColor)) {
    polygon(time, entropy, col = fillColor, border = fillColor)
  }
  # axis(2)

  # lines(time, entropy)

  if (drawBox) {
    box()
  }
}
