#' MarkovJumpR: R plotting routines for MarkovJump reconstructions
#'
#' MarkovJumpR is a new toy in the BEAST arsensal to facilitate phylogenetic
#' reconstruction using travel history information.
#'
#' @docType package
#' @name MarkovJumpR
#' @import dplyr
#'
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics abline axis box plot segments
#' @importFrom stats rnorm
#' @importFrom utils read.csv
NULL

#' Loads MarkovJump paths from a CSV file
#'
#' @param fileName  CSV file name
#'
#' @return MarkovJumps path object
#'
#' @export
loadPaths <- function(fileName) {
  paths <- read.csv(fileName, header = TRUE, stringsAsFactors = FALSE)
  result <- list(
    paths = paths,
    locations = unique(paths$location),
    minTime = min(paths$startTime, paths$endTime),
    maxTime = max(paths$startTime, paths$endTime)
  )
  return(result)
}

#' Plots MarkovJump paths for one taxon
#'
#' @param paths  MarkovJumps paths
#' @param locationMap  A \code{locationMap} object
#' @param minTime   Minimum time to plot
#' @param maxTime   Maximum time to plot
#' @param defaultColor Default color for dwell times
#' @param otherColor Color for locations not provided in \code{locationMap}
#' @param verticalColor Color of transitions (default `NULL` uses origin location color)
#' @param alpha Alpha value to overlay dwell times
#' @param yJitterSd Standard deviation of jitter around location
#' @param cex.labels Location text label size
#' @param addLocationLine Plot dashed line across locations
#' @param mustDisplayAllLocations Check that \code{locationMap} contains all locations in \code{paths}
#' @param xAt Time axis \code{at} values
#' @param xLabels Time axis \code{labels} values
#'
#' @return MarkovJumps path object
#'
#' @import dplyr
#'
#' @export
plotPaths <- function(paths,
                      locationMap = NULL,
                      minTime = min(paths$startTime,
                                    paths$endTime),
                      maxTime = max(paths$startTime,
                                    paths$endTime),
                      defaultColor = "black",
                      otherColor = "black",
                      verticalColor = NULL,
                      alpha = 0.1,
                      yJitterSd = 0.1,
                      cex.labels = 0.5,
                      addLocationLine = FALSE,
                      mustDisplayAllLocations = TRUE,
                      xAt = NULL,
                      xLabels = TRUE) {

  if (length(unique(paths$taxonId)) > 1) {
    stop("Multiple taxa are not yet supported")
  }

  if (is.null(locationMap)) {
    uniqueLocations <- unique(paths$location)
    locationMap <- data.frame(
      location = uniqueLocations,
      position = c(1:length(uniqueLocations)),
      color = rep(defaultColor, length(uniqueLocations))
    )
  } else {
    if (is.null(locationMap$color)) {
      locationMap$color <- defaultColor
    }
  }

  allLocationsMatch <- all(unique(paths$location) %in% locationMap$location)
  if (!allLocationsMatch) {
    if (mustDisplayAllLocations) {
      stop("All locations in `paths` are not represented in `locationMap`")
    } else {
      paths[!(paths$location %in% locationMap$location), "location"] <- "Other"
      locationMap <- rbind(locationMap,
                           data.frame(location = "Other",
                                      position = nrow(locationMap) + 1,
                                      color = otherColor))
    }
  }

  yJitter <- rnorm(n = nrow(paths),
                   mean = 0,
                   sd = yJitterSd)

  vcol2rgb <- Vectorize(col2rgb)

  paths <- paths %>% left_join(locationMap, by = "location") %>%
    mutate(position = .data$position + yJitter,
           color = rgb(vcol2rgb(.data$color)[1,], vcol2rgb(.data$color)[2,], vcol2rgb(.data$color)[3,],
                       alpha = alpha * 255, maxColorValue = 255))

  colorRgb <- col2rgb(defaultColor)
  defaultColor <- rgb(colorRgb[1],
               colorRgb[2],
               colorRgb[3],
               alpha = alpha * 255,
               maxColorValue = 255)

  plot(0, 0, xlim = c(minTime, maxTime),
       ylim = c(min(locationMap$position) - 0.5,
                max(locationMap$position) + 0.5),
       type = "n",
       ylab = "",
       xlab = "",
       axes = FALSE)

  axis(1, at = xAt, labels = xLabels)
  axis(2, at = locationMap$position,
       labels = locationMap$location,
       cex.axis = cex.labels,
       line = -0.75,
       las = 1, tick = FALSE)

  transitions <- paths %>% inner_join(
    paths, by = c("treeId", "taxonId",
               "startTime" = "endTime")) %>%
    select(-.data$endTime, -.data$startTime.y) %>%
    rename(time = .data$startTime)

  if (!is.null(verticalColor)) {
    transitions$color.x <- verticalColor
  }

  if (addLocationLine) {
    abline(
      h = locationMap$position,
      col = locationMap$color,
      lty = 2
    )
  }

  segments(x0 = transitions$time,
           x1 = transitions$time,
           y0 = transitions$position.x,
           y1 = transitions$position.y,
           col = transitions$color.x)

  segments(x0 = paths$startTime,
           x1 = paths$endTime,
           y0 = paths$position,
           y1 = paths$position,
           col = paths$color)

  box()
}

#' Create a \code{locationMap} ordered by location weights
#'
#' @param paths  MarkovJump paths for one or more taxa
#' @param heaviestAtTop Order locations starting at top of plot
#' @param minRelativeWeight Minimum relative weight to include in resulting \code{locationMap}
#'
#' @return \code{locationMap} object
#'
#' @export
orderByWeight <- function(paths,
                          heaviestAtTop = FALSE,
                          minRelativeWeight = 0) {

  location <- paths %>% group_by(.data$location) %>%
    summarize(weight = sum(abs(.data$endTime - .data$startTime))) %>%
    mutate(relativeWeight = .data$weight / sum(.data$weight)) %>%
    filter(.data$relativeWeight > minRelativeWeight) %>%
    arrange(.data$weight) %>% select(.data$location) %>% pull()

  if (length(location) == 0) {
    stop("No locations meet `minRelativeWeight` restriction")
  }

  position = c(1:length(location))
  if (!heaviestAtTop) {
    position <- rev(position)
  }

  map <- data.frame(location = location,
                    position = position)
  map$location <- as.character(map$location)

  return(map)
}
