#' Opens default web browser and loads a web page for length estimation and
#' testing.
#'
#' Function \code{startweb.esttest()} opens a web browser and loads the page
#' for length estimation and hypothesis testing of data read from a picture.
#' Use this function if you would like to compute the length of an interval,
#' which is the domain of uniform distribution, but the data are contaminated
#' with additive error. Sequence of actions on the web page is as follows:
#' \enumerate{
#'  \item Load a picture in JPG format
#'  \item Click on start and end point of the line passing through the observed
#'    object
#'  \item Set data preparation parameters
#'  \item Click on generate data so the data are prepared
#'  \item \itemize{
#'    \item set estimation parameters (see \code{\link{lengthest}}) and click
#'      on Estimate, or
#'    \item set testing parameters (see \code{\link{lengthtest}}) and click on
#'      Test.
#'    }
#'  }
#'
#'  Parameters that can be set on the web page are as follows:
#'
#'  \strong{Data parameters}
#'  \describe{
#'    \item{Levels of grey}{Number of colors (shades of grey) used in analysis.}
#'    \item{Box size}{The algorithm takes each pixel of a picture and maps it
#'      to box_size * box_size matrix. It is done in a way that the brightness
#'      of the observed pixel dictates the quantity of \emph{dots} in mentioned
#'      matrix. Distribution of dots in matrix is uniform. Ultimately, length
#'      estimation or testing is done on the set of the resulting matrices.}
#'    \item{Line thickness}{Width of the line, i.e. the maximum length between
#'      surrounding pixel and the drawn line so that pixel is to be taken into
#'      account for length estimation or hypothesis testing. All surrounding
#'      pixels are orthogonally projected on the line.}
#'    \item{Observed object is...}{Sets whether observed object is bright or
#'      dark.}
#'  }
#'  \strong{Estimation parameters}
#'  \describe{
#'    \item{Error distribution}{The type of the error distribution. Can be
#'      Gauss, Laplace or Student.}
#'    \item{Error standard deviation}{Estimation method for the error
#'      standard deviation. Can be Maximum Likelihood (ML) or the Method of
#'      Moments. If one does not want to estimate the deviation but to
#'      explicitly enter it, he should choose "Enter value" and enter the
#'      deviation in the lower field.}
#'    \item{Confidence level}{Confidence level of the confidence interval.}
#'  }
#'  \strong{Testing parameters}
#'  \describe{
#'    \item{\eqn{latex}{H_0}}{Specified null value being tested (measured in
#'      pixel width or percentage of image width).}
#'    \item{Alternative}{The alternative hypothesis. Can be less, greater or
#'      two-sided.}
#'  }
#'  \strong{Results}
#'  \itemize{
#'    \item Length expressed in pixel width, as well as in percentage of the
#'      whole image's width
#'    \item Standard deviation
#'    \item Method - Asymptotic distribution of MLE or LR statistic
#'    \item Confidence interval
#'    \item p-value of the test and the value of the test statistic (if
#'      hypothesis testing has been performed)
#'  }
#'
#' @section Note:
#' In order to have quadratic pixels on the screen, please use proportional
#' screen resolution. In the case of modern LCD (LED) displays, these are
#' usually native screen resolutions. If your display has aspect ratio
#' width:height = 16:9, these resolutions are 1280x720, 1600x900, 1920x1080,
#' etc. In the case od 16:10 display, use 1280x800, 1440x900, 1920x1200, etc.
#' If you use nonproportional screen resolution, pixels on the screen will not
#' be quadratic, so estimated values measured in pixels may not be correct.
#'
#' @examples
#' # open the web page for length estimation and hypothesis testing of an
#' # object shown in the picture
#' \donttest{startweb.esttest()}
#'
#' @export
startweb.esttest <- function() {
  opencpu::ocpu_start_server(preload = "LeArEst", on_startup = function(server){
    utils::browseURL(paste0(server, "/library/LeArEst/www/index_esttest.html"))
  })
}

# Utility function for web interface call
# Not intended for user calls!
#
#' @importFrom stats runif
#' @importFrom graphics hist
load.image <- function(generatedData = NULL, dataCenter = NULL, testType,
  picture, dataBright = "bright", error, var, varEst,
  confLevel, levelsOfGray, boxSize, thickness, nulla, unit,
  alternative, startX, startY, endX, endY) {
  # 0 < intensity < levelsOfGray-1 must apply: boxSize^2 >= levelsOfGray-1 !!!

  var <- as.numeric(var)
  conf.level <- as.numeric(confLevel)
  levelsOfGray <- as.numeric(levelsOfGray)
  boxSize <- as.numeric(boxSize)
  thickness <- as.numeric(thickness)
  null.a <- as.numeric(nulla)
  if (alternative == "t")
    alt <- "two.sided" else if (alternative == "g")
      alt <- "greater" else if (alternative == "l")
        alt <- "less"

  if (!is.null(generatedData) && generatedData != 0) {
    data <- as.double(unlist(strsplit(generatedData, split = ", ")))
    dataCenter <- as.double(dataCenter)
  }

  img_temp <- jpeg::readJPEG(picture)

  # if the image is color, convert it to greyscale first
  if (!is.na(dim(img_temp)[3])) {
    img <- img_temp[, , 1] * 0.33 + img_temp[, , 2] * 0.33 +
      img_temp[, , 3] * 0.33
  } else {
    img <- img_temp
  }
  rm(img_temp)

  # if we watch the dark area, make the negative of the image
  if (dataBright == "dark") {
    img <- 1 - img
  }

  img <- floor(levelsOfGray * img)

  height <- dim(img)[1]
  width <- dim(img)[2]

  if (unit == "perc") {
    null.a <- (null.a / 100 * width)
  }

  startPoint <- c(as.numeric(startX), as.numeric(startY))
  endPoint <- c(as.numeric(endX), as.numeric(endY))

  # if the endPoint is 'more left' then the startPoint - switch them
  if (endPoint[1] < startPoint[1]) {
    temp <- startPoint
    startPoint <- endPoint
    endPoint <- temp
  }

  # correction of the points, if they are too close to the picture borders
  if (startPoint[1] < thickness/2)
    startPoint[1] <- ceiling(thickness/2)
  if (startPoint[2] < thickness/2)
    startPoint[2] <- ceiling(thickness/2)
  if (startPoint[1] > height - thickness/2)
    startPoint[1] <- height - ceiling(thickness/2)
  if (startPoint[2] > width - thickness/2)
    startPoint[2] <- width - ceiling(thickness/2)
  if (endPoint[1] < thickness/2)
    endPoint[1] <- ceiling(thickness/2)
  if (endPoint[2] < thickness/2)
    endPoint[2] <- ceiling(thickness/2)
  if (endPoint[1] > height - thickness/2)
    endPoint[1] <- height - ceiling(thickness/2)
  if (endPoint[2] > width - thickness/2)
    endPoint[2] <- width - ceiling(thickness/2)

  # the angle between the line and the abscissa.
  fiLine <- atan((endPoint[2] - startPoint[2])/(endPoint[1] - startPoint[1]))

  # if 'Prepare Data' is clicked, calculate the vector 'data' and 'dataCenter'
  if (testType == "prepare") {
    dataMatrix = matrix(FALSE, nrow = height * boxSize, ncol = width * boxSize)

    for (j in 1:width) {
      for (i in 1:height) {
        x <- floor(runif(img[i, j], max = boxSize))
        y <- floor(runif(img[i, j], max = boxSize))
        dataMatrix[(i - 1) * boxSize + x, (j - 1) * boxSize + y] <- TRUE
      }
    }

    # vector of projected points' distances to the left point
    distancesVector = vector(mode = "numeric", length = sum(dataMatrix))

    # convert the line to Ax + By + C = 0 in the points area
    A = boxSize * (endPoint[2] - startPoint[2])
    B = boxSize * (startPoint[1] - endPoint[1])
    C = boxSize^2 * (startPoint[2] * (endPoint[1] - startPoint[1]) -
        startPoint[1] * (endPoint[2] - startPoint[2]))

    XStart <- floor((startPoint[1] - thickness/2) * boxSize)
    XEnd <- floor((endPoint[1] + thickness/2) * boxSize)
    YStart <- floor((min(startPoint[2], endPoint[2]) - thickness/2) * boxSize)
    YEnd <- floor((max(startPoint[2], endPoint[2]) + thickness/2) * boxSize)

    # line segment on the ordinate (x = 0)
    lineSegment <- -C/B

    denom <- sqrt(A^2 + B^2)
    cnt <- 1
    for (j in YStart:YEnd) {
      for (i in XStart:XEnd) {
        if (dataMatrix[i, j] == FALSE)
          next

        # calculate distance between the point and line,
        # if it is larger then thickness/2 => next
        d <- abs(A * i + B * j + C)/denom
        if (d > thickness * boxSize/2)
          next

        fiPoint <- atan((j - lineSegment)/i)
        pointBelow <- if (fiPoint < fiLine)
          1 else -1

        projX <- i - pointBelow * d * sin(fiLine)
        if ((projX < startPoint[1] * boxSize) || (projX > endPoint[1] * boxSize))
          next
        projY <- j + pointBelow * d * cos(fiLine)

        distancesVector[cnt] <- sqrt((projX - startPoint[1] * boxSize)^2 +
            (projY - startPoint[2] * boxSize)^2)

        cnt <- cnt + 1
      }
    }
    data <- distancesVector[distancesVector != 0]
    dataCenter <- mean(data)
    data <- (data - dataCenter)/boxSize

    hist(data)

    output <- "Data successfully prepared.<br/>"

    return(c(output, 0, 0, 0, 0, toString(data), dataCenter))

  } else if (testType == "est") {
    if (varEst == "value") {
      estResult <- lengthest(data, error = error, var = var,
        conf.level = conf.level)
    } else {
      estResult <- lengthest(data, error = error, var.est = varEst,
        conf.level = conf.level)
    }

    centerPoint <- c(startPoint[1] + cos(fiLine) * dataCenter/boxSize,
      startPoint[2] + sin(fiLine) * dataCenter/boxSize)
    leftPoint <- centerPoint - c(cos(fiLine) * estResult$radius,
      sin(fiLine) * estResult$radius)
    rightPoint <- centerPoint + c(cos(fiLine) * estResult$radius,
      sin(fiLine) * estResult$radius)

    if (varEst == "MM") {
      varianceCalcd <- " (MM estimated) "
    } else if (varEst == "ML") {
      varianceCalcd <- " (ML estimated) "
    } else {
      varianceCalcd <- " (explicitly given) "
    }

    percLength <- 100 * (2 * estResult$radius / width)

    output <- paste("Levels of grey: ", levelsOfGray, ", Box size: ", boxSize,
      ", Line thickness: ", thickness, "<br>Error distribution: ", error,
      ", Error standard deviation: ", varEst, ", Confidence level: ", confLevel,
      "<br><br><strong>Length</strong>: ", round(2 * estResult$radius, 2),
      " pixel width (", round(percLength, 2), "% of the image width)",
      "<br><strong>Standard deviation</strong>: ",
      round(sqrt(estResult$var.error), 2), varianceCalcd,
      "<br><strong>Method</strong>: ", estResult$method,
      "<br><strong>Confidence interval</strong>: ("
      , round(2 * estResult$conf.int[1], 2), ", ",
      round(2 * estResult$conf.int[2], 2), ")", sep = "")

    return(c(output, leftPoint[1], leftPoint[2],
      rightPoint[1], rightPoint[2], 0, 0))

  } else if (testType == "test") {
    if (varEst == "value") {
      estResult <- lengthtest(x = data, error = error, alternative = alt,
        var = var, null.a = null.a/2,
        conf.level = conf.level)
    } else {
      estResult <- lengthtest(x = data, error = error, alternative = alt,
        null.a = null.a/2, var.est = varEst,
        conf.level = conf.level)
    }

    centerPoint <- c(startPoint[1] + cos(fiLine) * dataCenter/boxSize,
      startPoint[2] + sin(fiLine) * dataCenter/boxSize)
    leftPoint <- centerPoint - c(cos(fiLine) * estResult$radius,
      sin(fiLine) * estResult$radius)
    rightPoint <- centerPoint + c(cos(fiLine) * estResult$radius,
      sin(fiLine) * estResult$radius)

    if (varEst == "MM") {
      varianceCalcd <- " (MM estimated) "
    } else if (varEst == "ML") {
      varianceCalcd <- " (ML estimated) "
    } else {
      varianceCalcd <- " (explicitly given) "
    }

    percLength <- 100 * (2 * estResult$radius / width)

    output <- paste("Levels of grey: ", levelsOfGray, ", Box size: ", boxSize,
      ", Line thickness: ", thickness, "<br>Error distribution: ", error,
      ", Error standard deviation: ", varEst, ", Confidence level: ", confLevel,
      "<br><br><strong>Estimated length</strong>: ",
      round(2 * estResult$radius, 2),
      " pixel width (", round(percLength, 2), "% of the image width)",
      "<br><strong>p value</strong>: ", round(estResult$p.value, 4),
      "<br><strong>Alternative</strong>: ", estResult$alternative,
      "<br><strong>T</strong>: ", round(estResult$tstat, 4),
      "<br><strong>Method</strong>: ", estResult$method,
      "<br><strong>Standard deviation</strong>: ",
      round(sqrt(estResult$var.error), 2), varianceCalcd, sep = "")

    if (estResult$alternative == "two.sided") {
      output <- paste(output, "<br><strong>Confidence interval</strong>: (",
        round(2 * estResult$conf.int[1], 2), ", ",
        round(2 * estResult$conf.int[2], 2), ")", sep = "")
    }

    return(c(output, leftPoint[1], leftPoint[2],
      rightPoint[1], rightPoint[2], 0, 0))
  }
}
