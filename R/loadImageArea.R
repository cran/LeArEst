#' Opens default web browser and loads a web page for area estimation.
#'
#' Function \code{startweb.area()} opens a web browser and loads the page
#' for area estimation of object shown in a picture. Use this function if you
#' can isolate a part of the picture with uniformly distributed dots on an
#' elliptical domain with unclear borders. Sequence of actions on the web page
#' is as follows:
#' \enumerate{
#'  \item Load a picture in JPG format
#'  \item Click on upper left and lower right corner of a rectangle surrounding
#'   observed object so the rectangle is drawn
#'  \item Set data and estimation parameters
#'  \item Click on Estimate
#'  }
#'
#' The area estimation algorithm takes many horizontal and vertical (if
#' "horizontal + vertical" slicing is selected) or star-shaped (if "star"
#' slicing is selected) slices of the object. Length estimation procedure is
#' conducted on each slice and in that way set of edge points is obtained.
#' Lastly, ellipse or circle is fitted on that set of points by function
#' \code{\link[conicfit]{EllipseDirectFit}} or
#' \code{\link[conicfit]{CircleFitByPratt}} from the package \code{conicfit}
#' and area of that ellipse or circle is returned as the result. The area is
#' measured in pixels, as well as percentage of the whole image.
#'
#' Parameters that can be set on the web page are as follows:
#'
#' \strong{Data parameters}
#'  \describe{
#'    \item{Levels of grey}{Number of colors (shades of grey) used in
#'      analysis.}
#'    \item{Box size}{The algorithm takes each pixel of a picture and maps it
#'      to box_size * box_size matrix. It is done in a way that the brightness
#'      of the observed pixel dictates the quantity of dots in mentioned
#'      matrix. Distribution of dots in matrix is uniform. Ultimately, length
#'      estimation is done on the set of the resulting matrices.}
#'    \item{Line thickness}{Width of the slice, i.e. the maximum length between
#'      surrounding pixel and the drawn line so that pixel is to be taken into
#'      account for length estimation. All surrounding pixels are orthogonally
#'      projected on the central line.}
#'    \item{Number of slices}{Number of slices after cutting in one direction.
#'      Defaults to 10. Slices are equally thick in both directions. Smaller
#'      number of cuts will be automatically applied for smaller dimension if
#'      the chosen rectangle is not a square.}
#'    \item{Slicing}{Sets slicing method for the edge point estimation. Can be
#'      "horizontal + vertical" or "star".}
#'    \item{Parallelization}{Sets whether to distribute area estimation on
#'      multiple CPU cores. If set to On, total number of cores - 1 are used.}
#'    \item{Object brightness}{Sets whether observed object is bright or
#'      dark.}
#'    \item{Represent object as}{Represent estimated object as an ellipse or as
#'      a circle.}
#'  }
#'  \strong{Estimation parameters}
#'  \describe{
#'    \item{Error distribution}{Type of the error distribution. Can be Gauss,
#'      Laplace or Student.}
#'    \item{Error standard deviation}{Estimation method for the error
#'      standard deviation. Can be Maximum Likelihood (ML) or the Method of
#'      Moments. If one does not want to estimate the deviation but to
#'      explicitly enter it, he should choose "Enter value" and enter the
#'      deviation in the lower field.}
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
#' # open the web page for area estimation of an object shown in the picture
#' \donttest{startweb.area()}
#'
#' @export
startweb.area <- function() {
  opencpu::ocpu_start_server(preload = "LeArEst", on_startup = function(server){
    utils::browseURL(paste0(server, "/library/LeArEst/www/index_area.html"))
  })
}

# Utility function for web interface call
# Not intended for user calls!
#
#' @importFrom stats runif
#' @importFrom stats kmeans
load.image.area <- function(picture, dataBright = "bright", representation,
  error, var, varEst, confLevel, levelsOfGray, boxSize, thickness, parallel,
  slicing, nrSlices, startX, startY, endX, endY) {
  # 0 < intensity < levelsOfGray-1 must apply: boxSize^2 >= levelsOfGray-1 !!!

  var <- as.numeric(var)
  conf.level <- as.numeric(confLevel)
  levelsOfGray <- as.numeric(levelsOfGray)
  boxSize <- as.numeric(boxSize)
  thickness <- as.numeric(thickness)
  nrSlices <- as.numeric(nrSlices)
  parallel <- as.logical(as.numeric(parallel))

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

  # prepare dataMatrix for the whole image
  dataMatrix = matrix(FALSE, nrow = height * boxSize, ncol = width * boxSize)
  set.seed(123)

  # for (j in 1:width) {
  #   for (i in 1:height) {
  for (j in startPoint[2]:endPoint[2]) {
    for (i in startPoint[1]:endPoint[1]) {
      x <- floor(runif(img[i, j], max = boxSize))
      y <- floor(runif(img[i, j], max = boxSize))
      dataMatrix[(i - 1) * boxSize + x, (j - 1) * boxSize + y] <- TRUE
    }
  } # TODO: za manje zauzece memorije smanjiti dimenziju dataMatrixa

  # divide selected square into two clusters to obtain limits for chopping
  kMeans <- kmeans(c(img[startPoint[1]:endPoint[1], startPoint[2]:endPoint[2]]),
    2, nstart = 100)
  clusterMatrix <- matrix(kMeans$cluster, nrow = endPoint[1] - startPoint[1] + 1,
    ncol = endPoint[2] - startPoint[2] + 1)

  # remove the outliers from clusterMatrix 'consecutive' tells how many
  # consecutive rows or columns should be 'an object' to trigger clusterStartPoint
  consecutive <- 20

  # suppose there is no object on the selected square's borders.
  # Taking that into account, adjust object to be in cluster 2.
  borderClusterPoints <- c(clusterMatrix[1, ], clusterMatrix[, dim(clusterMatrix)[2]],
    clusterMatrix[dim(clusterMatrix)[1], ], clusterMatrix[, 1])

  if (sum(borderClusterPoints == 2) > sum(borderClusterPoints == 1)) {
    clusterMatrix <- clusterMatrix + 1
    clusterMatrix[clusterMatrix == 3] <- 1
  }

  # image(clusterMatrix)

  maxrows <- apply(clusterMatrix, 1, max)
  maxcols <- apply(clusterMatrix, 2, max)

  for (i in 1:(length(maxrows) - consecutive + 1)) {
    if (sum(maxrows[i:(i + consecutive - 1)] == 2) == consecutive) {
      clusterMatrix[1:(i - 1), ] <- 1
      break
    }
  }

  for (i in (length(maxrows) - consecutive + 1):1) {
    if (sum(maxrows[i:(i + consecutive - 1)] == 2) == consecutive) {
      clusterMatrix[(i + consecutive - 1):length(maxrows), ] <- 1
      break
    }
  }

  for (i in 1:(length(maxcols) - consecutive + 1)) {
    if (sum(maxcols[i:(i + consecutive - 1)] == 2) == consecutive) {
      clusterMatrix[, 1:(i - 1)] <- 1
      break
    }
  }

  for (i in (length(maxcols) - consecutive + 1):1) {
    if (sum(maxcols[i:(i + consecutive - 1)] == 2) == consecutive) {
      clusterMatrix[, (i + consecutive - 1):length(maxcols)] <- 1
      break
    }
  }

  watchedPoints <- which(clusterMatrix == 2, arr.ind = TRUE)
  clusterStartPoint <- c(min(watchedPoints[, 1]), min(watchedPoints[, 2])) +
    startPoint  # -1
  clusterEndPoint <- c(max(watchedPoints[, 1]), max(watchedPoints[, 2])) +
    startPoint  # - 1

  ellipsePoints <- matrix(NA, nrow = 1000, ncol = 2)
  pointNr <- 1

  biggerDim <- max((clusterEndPoint[1] - clusterStartPoint[1]),
    (clusterEndPoint[2] - clusterStartPoint[2]))
  lineDistance <- biggerDim %/% nrSlices

  rectHeight <- endPoint[1] - startPoint[1]
  rectWidth <- abs(endPoint[2] - startPoint[2])

  currentX <- floor((clusterStartPoint[1] + lineDistance/2) * boxSize)

  start.time <- Sys.time()

  if (slicing == "hv") {
    # horizontal chopping
    if (parallel) {
      no_cores <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(no_cores)
      doParallel::registerDoParallel(cl)
      ellipsePointsH <- foreach::foreach (currentX = seq(
            from = floor((clusterStartPoint[1] + lineDistance/2) * boxSize),
            to = clusterEndPoint[1] * boxSize,
            by = lineDistance * boxSize)) %dopar% {
        distancesVector = vector(mode = "numeric",
          length = thickness * boxSize * rectWidth)
        distFromStart <- 0
        cnt <- 1

        for (currentY in min(startPoint[2] * boxSize, endPoint[2] * boxSize):
            max(startPoint[2] * boxSize, endPoint[2] * boxSize)) {
          leftestPoint <- currentX - floor(thickness * boxSize / 2)
          rightestPoint <- currentX + floor(thickness * boxSize / 2)
          sumThisY <- sum(dataMatrix[leftestPoint:rightestPoint, currentY])

          if (sumThisY > 0) {
            distancesVector[cnt:(cnt + sumThisY - 1)] <- distFromStart
          }

          distFromStart <- distFromStart + 1
          cnt <- cnt + sumThisY
        }

        data <- distancesVector[distancesVector != 0]
        dataCenter <- mean(data)
        data <- (data - dataCenter) / boxSize

        if (varEst == "value") {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = varEst,
            conf.level = conf.level)
        }

        centerY <- (min(startPoint[2], endPoint[2]) + dataCenter / boxSize)
        leftPoint <- floor(c(currentX / boxSize,
                             centerY - unname(estResult$radius)))
        rightPoint <- floor(c(currentX / boxSize,
                              centerY + unname(estResult$radius)))

        # rm(distancesVector, data, dataCenter, data, estResult, centerY)

        c(leftPoint, rightPoint)
      }
      # razvrti tu listu u matricu:
      ellipsePoints <- matrix(unlist(ellipsePointsH), ncol = 2, byrow = TRUE)
    } else {
      while (currentX < clusterEndPoint[1] * boxSize) {
        distancesVector = vector(mode = "numeric",
          length = thickness * boxSize * rectWidth)
        distFromStart <- 0
        cnt <- 1

        for (currentY in min(startPoint[2] * boxSize, endPoint[2] * boxSize):
            max(startPoint[2] * boxSize, endPoint[2] * boxSize)) {
          leftestPoint <- currentX - floor(thickness * boxSize / 2)
          rightestPoint <- currentX + floor(thickness * boxSize / 2)
          sumThisY <- sum(dataMatrix[leftestPoint:rightestPoint, currentY])

          if (sumThisY > 0) {
            distancesVector[cnt:(cnt + sumThisY - 1)] <- distFromStart
          }

          distFromStart <- distFromStart + 1
          cnt <- cnt + sumThisY
        }

        data <- distancesVector[distancesVector != 0]
        dataCenter <- mean(data)
        data <- (data - dataCenter) / boxSize

        if (varEst == "value") {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = varEst,
            conf.level = conf.level)
        }

        centerY <- (min(startPoint[2], endPoint[2]) + dataCenter / boxSize)
        leftPoint <- floor(c(currentX / boxSize, centerY - estResult$radius))
        rightPoint <- floor(c(currentX / boxSize, centerY + estResult$radius))

        ellipsePoints[(2 * pointNr - 1), ] <- leftPoint
        ellipsePoints[(2 * pointNr), ] <- rightPoint

        currentX <- currentX + lineDistance * boxSize
        pointNr <- pointNr + 1
      }
    }

    currentY <- floor((clusterStartPoint[2] + lineDistance / 2) * boxSize)

    # vertical chopping
    if (parallel) {
      ellipsePointsV <- foreach::foreach (currentY = seq(
            from = floor((clusterStartPoint[2] + lineDistance / 2) * boxSize),
            to = clusterEndPoint[2] * boxSize,
            by = lineDistance * boxSize)) %dopar% {
        distancesVector = vector(mode = "numeric",
          length = thickness * boxSize * rectHeight)
        distFromStart <- 0
        cnt <- 1

        for (currentX in min(startPoint[1] * boxSize, endPoint[1] * boxSize):
            max(startPoint[1] * boxSize, endPoint[1] * boxSize)) {
          leftestPoint <- currentY - floor(thickness * boxSize / 2)
          rightestPoint <- currentY + floor(thickness * boxSize / 2)
          sumThisX <- sum(dataMatrix[currentX, leftestPoint:rightestPoint])

          if (sumThisX > 0) {
            distancesVector[cnt:(cnt + sumThisX - 1)] <- distFromStart
          }

          distFromStart <- distFromStart + 1
          cnt <- cnt + sumThisX
        }

        data <- distancesVector[distancesVector != 0]
        dataCenter <- mean(data)
        data <- (data - dataCenter) / boxSize

        if (varEst == "value") {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = varEst,
            conf.level = conf.level)
        }

        centerX <- (min(startPoint[1], endPoint[1]) + dataCenter / boxSize)
        upPoint <- floor(c(centerX - estResult$radius, currentY / boxSize))
        downPoint <- floor(c(centerX + estResult$radius, currentY / boxSize))

        # rm(distancesVector, data, dataCenter, data, estResult, centerY)

        c(upPoint, downPoint)
      }
      ellipsePointsV <- matrix(unlist(ellipsePointsV), ncol = 2, byrow = TRUE)
      ellipsePoints <- as.matrix(merge(ellipsePoints, ellipsePointsV,
                                       all = TRUE))
    } else {
      while (currentY < clusterEndPoint[2] * boxSize) {
        distancesVector = vector(mode = "numeric",
          length = thickness * boxSize * rectHeight)
        distFromStart <- 0
        cnt <- 1

        for (currentX in min(startPoint[1] * boxSize, endPoint[1] * boxSize):
            max(startPoint[1] * boxSize, endPoint[1] * boxSize)) {
          leftestPoint <- currentY - floor(thickness * boxSize / 2)
          rightestPoint <- currentY + floor(thickness * boxSize / 2)
          sumThisX <- sum(dataMatrix[currentX, leftestPoint:rightestPoint])

          if (sumThisX > 0) {
            distancesVector[cnt:(cnt + sumThisX - 1)] <- distFromStart
          }

          distFromStart <- distFromStart + 1
          cnt <- cnt + sumThisX
        }

        data <- distancesVector[distancesVector != 0]
        dataCenter <- mean(data)
        data <- (data - dataCenter) / boxSize

        if (varEst == "value") {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = varEst,
            conf.level = conf.level)
        }

        centerX <- (min(startPoint[1], endPoint[1]) + dataCenter / boxSize)
        upPoint <- floor(c(centerX - estResult$radius, currentY / boxSize))
        downPoint <- floor(c(centerX + estResult$radius, currentY / boxSize))

        ellipsePoints[(2 * pointNr - 1), ] <- upPoint
        ellipsePoints[(2 * pointNr), ] <- downPoint

        currentY <- currentY + lineDistance * boxSize
        pointNr <- pointNr + 1
      }
    }
  } else if (slicing == "star") {
    # common part (parallel and serial):

    # create matrix of "TRUE" locations in dataMatrix, filter it to data
    # cluster and translate in a way that clusterStartPoint * boxSize is the
    # origin of the coordinate system
    mydatafull <- as.data.frame(which(dataMatrix == TRUE, arr.ind=T))
    colnames(mydatafull) <- c("x", "y")

    mydata <- mydatafull[mydatafull$x > clusterStartPoint[1] * boxSize &
        mydatafull$y > clusterStartPoint[2] * boxSize &
        mydatafull$x < clusterEndPoint[1] * boxSize &
        mydatafull$y < clusterEndPoint[2] * boxSize, ]
    mydata$x <- mydata$x - (clusterStartPoint[1] * boxSize)
    mydata$y <- mydata$y - (clusterStartPoint[2] * boxSize)

    rm(mydatafull)

    # translate input points so the origin is data mean
    # centerx <- mean (mydata$x)
    # centery <- mean (mydata$y)

    # translate input points so the origin is center point of the data cluster
    centerx <- abs((clusterEndPoint[1] - clusterStartPoint[1]) * boxSize %/% 2)
    centery <- abs((clusterEndPoint[2] - clusterStartPoint[2]) * boxSize %/% 2)

    mydata <- mydata - c(centerx, centery)

    if (parallel) {
      no_cores <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(no_cores)
      doParallel::registerDoParallel(cl)

      ellipsePoints <- foreach::foreach (currAngle = seq(
        from = 0, to = pi - 0.001, by = pi / nrSlices)) %dopar% {

          mydataRot <- mydata
          if (currAngle != 0) {
            mydataRot$x <- mydata$x * cos(currAngle) - mydata$y * sin(currAngle)
            mydataRot$y <- mydata$y * cos(currAngle) + mydata$x * sin(currAngle)
          }

          minRoty <- min(mydataRot$y)
          linedata <- mydataRot[mydataRot$x > (-thickness * boxSize / 2) &
              mydataRot$x < (thickness * boxSize / 2), ]

          distancesVector <- sapply(linedata$y, '-', minRoty)
          dataCenter <- mean(distancesVector)
          data <- distancesVector - dataCenter

          if (varEst == "value") {
            estResult <- lengthest(data, error = error, var = var,
              conf.level = conf.level)
          } else {
            estResult <- lengthest(data, error = error, var.est = varEst,
              conf.level = conf.level)
          }

          downPointRot <- c(0, dataCenter + minRoty - unname(estResult$radius))
          upPointRot <- c(0, dataCenter + minRoty + unname(estResult$radius))

          # rotate calculated points back and undo translate:
          downPointX <- downPointRot[1] * cos(currAngle) +
                        downPointRot[2] * sin(currAngle) +
                        centerx + clusterStartPoint[1] * boxSize
          downPointY <- downPointRot[2] * cos(currAngle) -
                        downPointRot[1] * sin(currAngle) +
                        centery + clusterStartPoint[2] * boxSize
          upPointX <- upPointRot[1] * cos(currAngle) +
                      upPointRot[2] * sin(currAngle) +
                      centerx + clusterStartPoint[1] * boxSize
          upPointY <- upPointRot[2] * cos(currAngle) -
                      upPointRot[1] * sin(currAngle) +
                      centery + clusterStartPoint[2] * boxSize

          upPoint <- c(upPointX, upPointY) %/% boxSize
          downPoint <- c(downPointX, downPointY) %/% boxSize

          c(upPoint, downPoint)
        }
      # cobvert that list to a matrix:
      ellipsePoints <- matrix(unlist(ellipsePoints), ncol = 2, byrow = TRUE)
    } else {
      for (currAngle in seq(0, pi-0.001, by = pi / nrSlices)) {
        # rotate the points by currAngle
        mydataRot <- mydata
        if (currAngle != 0) {
          mydataRot$x <- mydata$x * cos(currAngle) - mydata$y * sin(currAngle)
          mydataRot$y <- mydata$y * cos(currAngle) + mydata$x * sin(currAngle)
        }

        minRoty <- min(mydataRot$y)
        linedata <- mydataRot[mydataRot$x > (-thickness * boxSize / 2) &
            mydataRot$x < (thickness * boxSize / 2), ]

        distancesVector <- sapply(linedata$y, '-', minRoty)
        dataCenter <- mean(distancesVector)
        data <- distancesVector - dataCenter

        if (varEst == "value") {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = varEst,
            conf.level = conf.level)
        }

        downPointRot <- c(0, dataCenter + minRoty - unname(estResult$radius))
        upPointRot <- c(0, dataCenter + minRoty + unname(estResult$radius))

        # rotate calculated points back and undo translate:
        downPointX <- downPointRot[1] * cos(currAngle) +
                      downPointRot[2] * sin(currAngle) +
                      centerx + clusterStartPoint[1] * boxSize
        downPointY <- downPointRot[2] * cos(currAngle) -
                      downPointRot[1] * sin(currAngle) +
                      centery + clusterStartPoint[2] * boxSize
        upPointX <- upPointRot[1] * cos(currAngle) +
                    upPointRot[2] * sin(currAngle) +
                    centerx + clusterStartPoint[1] * boxSize
        upPointY <- upPointRot[2] * cos(currAngle) -
                    upPointRot[1] * sin(currAngle) +
                    centery + clusterStartPoint[2] * boxSize

        upPoint <- c(upPointX, upPointY) %/% boxSize
        downPoint <- c(downPointX, downPointY) %/% boxSize

        ellipsePoints[(2 * pointNr - 1), ] <- upPoint
        ellipsePoints[(2 * pointNr), ] <- downPoint

        pointNr <- pointNr + 1
      }
    }
  }

  end.time <- Sys.time()
  evalTime <- round(end.time - start.time, 4)

  ellipsePointsClean <- ellipsePoints[rowSums(is.na(ellipsePoints)) != 2, ]
  resultPoints <- paste(toString(t(ellipsePointsClean), sep = ", "))

  if (representation == 'ellipse') {
    ellipseDirect <- conicfit::EllipseDirectFit(ellipsePointsClean)
    ellipseDirectG <- conicfit::AtoG(ellipseDirect)$ParG
    ellipseArea <- round(ellipseDirectG[3] * ellipseDirectG[4] * pi, 2)

    # center of the ellipse and two points at the axes cross
    C = c(ellipseDirectG[1], ellipseDirectG[2])
    A = c(C[1] - ellipseDirectG[3] * cos(ellipseDirectG[5]),
          C[2] - ellipseDirectG[3] * sin(ellipseDirectG[5]))
    B = c(C[1] + ellipseDirectG[4] * sin(ellipseDirectG[5]),
          C[2] - ellipseDirectG[4] * cos(ellipseDirectG[5]))

    percArea <- ellipseArea / (height * width ) * 100

    output <- paste("Levels of grey: ", levelsOfGray,
      ", Box size: ", boxSize,
      ", Line thickness: ", thickness,
      ", Error distribution: ", error,
      "<br><strong>Area</strong>: ", ellipseArea, " pixels (",
      round(percArea, 2), "% of the image area)<br>",
      "Ellipse center (", round(ellipseDirectG[1]), ", ",
      round(ellipseDirectG[2]), "), ",
      "semiaxes: (", round(ellipseDirectG[3], 2), ", ",
      round(ellipseDirectG[4], 2) , "), ",
      "angle: ", round(ellipseDirectG[5], 2),
      "<br>A = (", round(A[1]), ", ", round(A[2]), "), ",
      "B = (", round(B[1]), ", ", round(B[2]), ")",
      "<br>Parallel: ", parallel, " Slicing: ", slicing, " Time: ", evalTime,
      sep = "")
  } else if (representation == 'circle') {
    circle <- conicfit::CircleFitByPratt(ellipsePointsClean)
    circleArea <- circle[3]**2 * pi;

    percArea <- circleArea / (height * width ) * 100

    output <- paste("Levels of grey: ", levelsOfGray,
      ", Box size: ", boxSize,
      ", Line thickness: ", thickness,
      ", Error distribution: ", error,
      "<br><strong>Area</strong>: ", circleArea, " pixels (",
      round(percArea, 2), "% of the image area)<br>",
      "Circle center (", round(circle[1]), ", ", round(circle[2]), "), ",
      "radius: ", round(circle[3]),
      "<br>Parallel: ", parallel, " Slicing: ", slicing, " Time: ", evalTime,
      sep = "")

    # show circle as ellipse for drawing purposes
    ellipseArea <- circleArea
    ellipseDirectG <- c(circle[1], circle[2], circle[3], circle[3], 0)
  }

  # return(c(output, ellipseArea, ellipseDirectG, substring(resultPoints, 3)))
  return(c(output, ellipseArea, ellipseDirectG, resultPoints))
}

#' Performs area estimation of the numerically described object in plane.
#'
#' Use this function if you have a data set of uniformly distributed points on
#' an elliptical domain in the plane but captured with additive errors. The
#' estimation algorithm takes many horizontal and vertical, or star-shaped
#' slices of the object. Length estimation procedure is conducted on each slice
#' and in that way the set of edge points is obtained. An ellipse or a circle
#' is fitted to these edge points by function
#' \code{\link[conicfit]{EllipseDirectFit}} or
#' \code{\link[conicfit]{CircleFitByPratt}} from the package \code{conicfit}
#' and its semi-axes and area are returned as a result. Function optionally
#' plots input points, calculated edge points and the resulting ellipse or
#' circle.
#'
#' @param data Two-column data matrix containing the points that describe
#'   observed object. First column represents \emph{x} coordinate of the point,
#'   while second column represents \emph{y} coordinate.
#' @param nrSlices Number of slices applied for plain data cutting. Defaults to
#'   10.
#' @param error A character string specifying the error distribution. Must be
#'   one of "laplace", "gauss" or "student". Can be abbreviated.
#' @param var.est A character string specifying the method of error variance
#'   estimation. Must be given if \code{var} is not given. Can be "MM" (Method
#'   of Moments) or "ML" (Maximum Likelihood).
#' @param var Explicit error variance. Needs to be given if \code{var.est} is
#'   not given.
#' @param plot Logical parameter (TRUE or FALSE) that determines whether to
#'   plot given object, calculated edge points and the resulting ellipse.
#'   Defaults to FALSE.
#' @param representation A character string specifying the shape of an observed
#'   object. Can be "ellipse" or "circle". Can be abbreviated.
#' @param parallel Logical parameter (TRUE or FALSE) that determines whether to
#'   perform estimation procedure in a parallel manner. Can shorten
#'   estimation time if many border points need to be calculated. Defaults to
#'   FALSE.
#' @param slicing A character string specifying the method of slicing. Can be
#'  "hv" (horizontal and vertical slicing) or "star" (star-shaped slicing). Can
#'   be abbreviated.
#'
#' @return List containing:
#'   \itemize{
#'     \item area: Estimated area of the object,
#'     \item points: Set of calculated object's edge points,
#'     \item semiaxes: Resulting ellipse's semi-axes or circle radius.
#'   }
#'
#' @examples
#' # load a data set representing the ellipse with additive Gaussian error,
#' # run area estimation on it, and plot the results
#' inputfile <- system.file("extdata", "ellipse_3_4_0.1_gauss.txt", package = "LeArEst")
#' inputdata <- read.table(inputfile)
#' area <- areaest(inputdata, error = "gauss", var.est = "ML", plot = TRUE,
#'                 slicing = "hv", representation = "ellipse")
#'
#' # load a data set representing the ellipse with additive Laplacian error,
#' # run area estimation on it, and plot the results
#' inputfile <- system.file("extdata", "ellipse_3_4_0.1_laplace.txt", package = "LeArEst")
#' inputdata <- read.table(inputfile)
#' area <- areaest(inputdata, error = "laplace", var = 0.1, nrSlices = 5, plot = TRUE,
#'                 slicing = "star", representation = "ellipse")
#'
#' @importFrom utils read.csv
#' @importFrom graphics points
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#'
#' @export
areaest <- function(data,
  nrSlices = 10,
  error = c("laplace", "gauss", "student"),
  var.est = c("MM", "ML"),
  var = NULL,
  plot = FALSE,
  parallel = FALSE,
  slicing = c("hv", "star"),
  representation = c("ellipse", "circle")) {
  # convert matrix to dataframe
  mydata <- data.frame(data)
  colnames(mydata) <- c("x", "y")

  minx <- min(mydata$x)
  maxx <- max(mydata$x)
  miny <- min(mydata$y)
  maxy <- max(mydata$y)

  ellipsePoints <- matrix(NA, nrow = nrSlices * 4, ncol = 2)
  pointNr <- 1

  biggerDim <- max((maxx - minx), (maxy - miny))
  lineDistance <- biggerDim / nrSlices

  currentX <- minx + lineDistance / 2
  currentY <- miny + lineDistance / 2

  conf.level <- 0.95

  start.time <- Sys.time()

  # horizontal chopping
  if (slicing == "hv") {
    if (parallel) {
      no_cores <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(no_cores)
      doParallel::registerDoParallel(cl)
      ellipsePointsH <- foreach::foreach (currentX = seq(from = minx + lineDistance / 2,
                                          to = maxx, by = lineDistance)) %dopar% {
        linedata <- mydata[mydata$x > (currentX - lineDistance / 2) &
                           mydata$x < (currentX + lineDistance / 2), ]

        distancesVector <- sapply(linedata$y, '-', miny)

        # if (length(distancesVector) < 10) {
        #   c(NA, NA, NA, NA)
        # }

        dataCenter <- mean(distancesVector)
        data <- distancesVector - dataCenter

        if (!is.null(var)) {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = var.est,
            conf.level = conf.level)
        }

        leftPoint <- c(currentX, dataCenter + miny - unname(estResult$radius))
        rightPoint <- c(currentX, dataCenter + miny + unname(estResult$radius))

        # rm(linedata, distancesVector, dataCenter, data, estResult)

        c(leftPoint, rightPoint)
      }
      # razvrti tu listu u matricu:
      ellipsePoints <- matrix(unlist(ellipsePointsH), ncol = 2, byrow = TRUE)
    } else {
      while (currentX < maxx) {
        linedata <- mydata[mydata$x > (currentX - lineDistance / 2) &
                           mydata$x < (currentX + lineDistance / 2), ]

        distancesVector <- sapply(linedata$y, '-', miny)

        # if (length(distancesVector) < 10) {
        #   currentX <- currentX + lineDistance
        #   next
        # }

        dataCenter <- mean(distancesVector)
        data <- distancesVector - dataCenter

        if (!is.null(var)) {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = var.est,
            conf.level = conf.level)
        }

        leftPoint <- c(currentX, dataCenter + miny - estResult$radius)
        rightPoint <- c(currentX, dataCenter + miny + estResult$radius)

        ellipsePoints[(2 * pointNr - 1), ] <- leftPoint
        ellipsePoints[(2 * pointNr), ] <- rightPoint

        currentX <- currentX + lineDistance
        pointNr <- pointNr + 1
      }
    }

    # vertical chopping
    if (parallel) {
      ellipsePointsV <- foreach::foreach(currentY = seq(from = miny + lineDistance / 2,
                                         to = maxy, by = lineDistance)) %dopar% {
        linedata <- mydata[mydata$y > (currentY - lineDistance / 2) &
                           mydata$y < (currentY + lineDistance / 2), ]

        distancesVector <- sapply(linedata$x, '-', minx)

        # if (length(distancesVector) < 10) {
        #   c(NA, NA, NA, NA)
        # }

        dataCenter <- mean(distancesVector)
        data <- distancesVector - dataCenter

        if (!is.null(var)) {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = var.est,
            conf.level = conf.level)
        }

        leftPoint <- c(dataCenter + minx - unname(estResult$radius), currentY)
        rightPoint <- c(dataCenter + minx + unname(estResult$radius), currentY)

        # rm(linedata, distancesVector, dataCenter, data, estResult)

        c(leftPoint, rightPoint)
      }
      ellipsePointsV <- matrix(unlist(ellipsePointsV), ncol = 2, byrow = TRUE)
      ellipsePoints <- as.matrix(merge(ellipsePoints, ellipsePointsV,
                                 all = TRUE))
    } else {
      while (currentY < maxy) {
        linedata <- mydata[mydata$y > (currentY - lineDistance / 2) &
                           mydata$y < (currentY + lineDistance / 2), ]

        distancesVector <- sapply(linedata$x, '-', minx)

        # if (length(distancesVector) < 10) {
        #   currentY <- currentY + lineDistance
        #   next
        # }

        dataCenter <- mean(distancesVector)
        data <- distancesVector - dataCenter

        if (!is.null(var)) {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = var.est,
            conf.level = conf.level)
        }

        leftPoint <- c(dataCenter + minx - estResult$radius, currentY)
        rightPoint <- c(dataCenter + minx + estResult$radius, currentY)

        ellipsePoints[(2 * pointNr - 1), ] <- leftPoint
        ellipsePoints[(2 * pointNr), ] <- rightPoint

        currentY <- currentY + lineDistance
        pointNr <- pointNr + 1
      }
    }
  } else if (slicing == "star") {
    # common part:

    # translate input points so the origin is data mean
    centerx <- mean (mydata$x)
    centery <- mean (mydata$y)

    mydata$x <- mydata$x - centerx
    mydata$y <- mydata$y - centery

    lineDistance <- biggerDim / 10

    if (parallel) {
      no_cores <- parallel::detectCores() - 1
      cl <- parallel::makeCluster(no_cores)
      doParallel::registerDoParallel(cl)

      ellipsePoints <- foreach::foreach (currAngle = seq(
        from = 0, to = pi - 0.001, by = pi / nrSlices)) %dopar% {
          # rotate the points by currAngle
          mydataRot <- mydata
          if (currAngle != 0) {
            mydataRot$x <- mydata$x * cos(currAngle) - mydata$y * sin(currAngle)
            mydataRot$y <- mydata$y * cos(currAngle) + mydata$x * sin(currAngle)
          }

          minRoty <- min(mydataRot$y)
          linedata <- mydataRot[mydataRot$x > (-lineDistance / 2) &
              mydataRot$x < (lineDistance / 2), ]
          distancesVector <- sapply(linedata$y, '-', minRoty)
          dataCenter <- mean(distancesVector)
          data <- distancesVector - dataCenter

          if (!is.null(var)) {
            estResult <- lengthest(data, error = error, var = var,
              conf.level = conf.level)
          } else {
            estResult <- lengthest(data, error = error, var.est = var.est,
              conf.level = conf.level)
          }

          downPointRot <- c(0, dataCenter + minRoty - unname(estResult$radius))
          upPointRot <- c(0, dataCenter + minRoty + unname(estResult$radius))

          # rotate calculated points back and undo translate:
          downPointX <- downPointRot[1] * cos(currAngle) +
                        downPointRot[2] * sin(currAngle) + centerx
          downPointY <- downPointRot[2] * cos(currAngle) -
                        downPointRot[1] * sin(currAngle) + centery
          upPointX <- upPointRot[1] * cos(currAngle) +
                      upPointRot[2] * sin(currAngle) + centerx
          upPointY <- upPointRot[2] * cos(currAngle) -
                      upPointRot[1] * sin(currAngle) + centery

          upPoint <- c(upPointX, upPointY)
          downPoint <- c(downPointX, downPointY)

          c(upPoint, downPoint)
        }
      # razvrti tu listu u matricu:
      ellipsePoints <- matrix(unlist(ellipsePoints), ncol = 2, byrow = TRUE)
    } else {
      for (currAngle in seq(0, pi - 0.001, by = pi / nrSlices)) {
        # rotate the points by currAngle
        mydataRot <- mydata
        if (currAngle != 0) {
          mydataRot$x <- mydata$x * cos(currAngle) - mydata$y * sin(currAngle)
          mydataRot$y <- mydata$y * cos(currAngle) + mydata$x * sin(currAngle)
        }

        minRoty <- min(mydataRot$y)
        linedata <- mydataRot[mydataRot$x > (-lineDistance / 2) &
            mydataRot$x < (lineDistance / 2), ]
        distancesVector <- sapply(linedata$y, '-', minRoty)
        dataCenter <- mean(distancesVector)
        data <- distancesVector - dataCenter

        if (!is.null(var)) {
          estResult <- lengthest(data, error = error, var = var,
            conf.level = conf.level)
        } else {
          estResult <- lengthest(data, error = error, var.est = var.est,
            conf.level = conf.level)
        }

        downPointRot <- c(0, dataCenter + minRoty - unname(estResult$radius))
        upPointRot <- c(0, dataCenter + minRoty + unname(estResult$radius))

        # rotate calculated points back and undo translate:
        downPointX <- downPointRot[1] * cos(currAngle) +
                      downPointRot[2] * sin(currAngle) + centerx
        downPointY <- downPointRot[2] * cos(currAngle) -
                      downPointRot[1] * sin(currAngle) + centery
        upPointX <- upPointRot[1] * cos(currAngle) +
                    upPointRot[2] * sin(currAngle) + centerx
        upPointY <- upPointRot[2] * cos(currAngle) -
                    upPointRot[1] * sin(currAngle) + centery

        ellipsePoints[(2 * pointNr - 1), ] <- c(downPointX, downPointY)
        ellipsePoints[(2 * pointNr), ] <- c(upPointX, upPointY)

        pointNr <- pointNr + 1
      }
    }
    # undo translate for original points
    mydata$x <- mydata$x + centerx
    mydata$y <- mydata$y + centery
  }

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat ("Time for chopping:", time.taken)

  ellipsePointsClean <- ellipsePoints[rowSums(is.na(ellipsePoints)) != 2, ]

  if (representation == "ellipse") {
    ellipseDirect <- conicfit::EllipseDirectFit(ellipsePointsClean)
    ellipseDirectG <- conicfit::AtoG(ellipseDirect)$ParG

    ellipseArea <- ellipseDirectG[3] * ellipseDirectG[4] * pi
  } else if (representation == "circle") {
    circle <- conicfit::CircleFitByPratt(ellipsePointsClean)
    circleArea <- circle[3]**2 * pi;

    # show circle as ellipse for drawing purposes
    ellipseArea <- circleArea
    ellipseDirectG <- c(circle[1], circle[2], circle[3], circle[3], 0)
  }

  if (plot) {
    xLim <- range(min(minx, ellipsePointsClean[, 1]),
      max(maxx, ellipsePointsClean[, 1]))
    yLim <- range(min(miny, ellipsePointsClean[, 2]),
      max(maxy, ellipsePointsClean[, 2]))

    xyDirect<-conicfit::calculateEllipse(ellipseDirectG[1], ellipseDirectG[2],
      ellipseDirectG[3], ellipseDirectG[4], 180 / pi * ellipseDirectG[5])

    plot(xyDirect[, 1], xyDirect[, 2], type = 'l', col = 'cyan',
      xlab = "x", ylab = "y", lwd = 3, xlim = xLim, ylim = yLim)
    points(mydata)
    points(ellipsePointsClean, col = "red", lwd = 7)
  }

  return (list(
    "area" = ellipseArea,
    "points" = ellipsePointsClean,
    "semiaxes" = c(ellipseDirectG[3], ellipseDirectG[4])
  ))
}
