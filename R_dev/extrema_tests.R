

# # Update R version.
# install.packages('installr')
# library(installr)
#
# updateR()


# # Testing functions to locate local optima.
#
# install.packages('ggpmisc')
# library(ggpmisc)
#
# test_1 <- c(1,1,2,3,3,2,1,1,5,6,7,7,4,2)
#
# peak_locs <- ggpmisc:::find_peaks(x = test_1,
#                                   ignore_threshold = 0,
#                                   span = 3,
#                                   # strict = TRUE,
#                                   strict = FALSE
#                                   )
#
# seq(1, length(test_1))[peak_locs]
# test_1[peak_locs]
#
#
#
# # Try a 2D example
#
# test_mat <- as.array(t(t(test_1)) %*% t(test_1))


# peak_locs <- ggpmisc:::find_peaks(x = test_mat,
#                                   ignore_threshold = 0,
#                                   span = 3,
#                                   # strict = TRUE,
#                                   strict = FALSE)
#
# seq(1, length(test_mat))[peak_locs]
# test_mat[peak_locs]




# Another approach.
# https://stackoverflow.com/questions/11059104/given-a-2d-numeric-height-map-matrix-in-r-how-can-i-find-all-local-maxima

install.packages('raster')
library(raster)

## Construct an example matrix
set.seed(444)
msize <- 10
x <- matrix(sample(seq_len(msize), msize^2, replace=TRUE), ncol=msize)

## Convert it to a raster object
r <- raster(x)
extent(r) <- extent(c(0, msize, 0, msize) + 0.5)

## Find the maximum value within the 9-cell neighborhood of each cell
f <- function(X) max(X, na.rm=TRUE)
ww <- matrix(1, nrow=3, ncol=3) ## Weight matrix for cells in moving window
localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)

## Does each cell have the maximum value in its neighborhood?
r2 <- r==localmax

## Get x-y coordinates of those cells that are local maxima
maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
head(maxXY)


find_local_max <- function(x) {

  ## Convert it to a raster object
  r <- raster(x)
  # extent(r) <- extent(c(0, nrow(x), 0, ncol(x)) + 0.5)
  extent(r) <- extent(c(0, ncol(x), 0, nrow(x)) + 0.5)

  ## Find the maximum value within the 9-cell neighborhood of each cell
  # f <- function(X) max(X, na.rm = TRUE)
  # ww <- matrix(1, nrow = 3, ncol = 3) ## Weight matrix for cells in moving window
  localmax <- focal(r,
                    fun = function(X) max(X, na.rm = TRUE),
                    w = matrix(1, nrow = 3, ncol = 3),
                    pad = TRUE, padValue = NA)

  ## Does each cell have the maximum value in its neighborhood?
  r2 <- r == localmax

  # Which(r2 == 1)

  ## Get x-y coordinates of those cells that are local maxima
  maxXY <- xyFromCell(r2, Which(r2 == 1, cells = TRUE))


  # Obtain the proper y-coordinates.
  # Notice the switched coordinates.
  loc_max_out <- list(
    # y_ind = ncol(x) - maxXY[, 'y'] + 1,
    row = nrow(x) - maxXY[, 'y'] + 1,
    col = maxXY[, 'x'],
    value = rep(NA, nrow(maxXY))
  )

  # Append value of objective at local maxima.
  for (i in 1:nrow(maxXY)) {
    loc_max_out$value[i] <- x[loc_max_out$row[i], loc_max_out$col[i]]
    # loc_max_out$z_ind[i] <- x[loc_max_out$y_ind[i], loc_max_out$x_ind[i]]
  }


  return(loc_max_out)
}


test_1 <- c(1,1,2,3,3,2,1,1,5,6,7,7,4,2)
test_mat <- as.array(t(t(test_1)) %*% t(test_1))
test_mat <- test_mat + matrix(0.2*runif(length(test_1)*length(test_1)),
                              nrow = length(test_1), ncol = length(test_1))


loc_max <- find_local_max(x = test_mat)

# test_mat[loc_max]

# Case with global max.



test_2 <- c(1,1,2,2,3,3,4,3,2,2)
test_mat <- as.array(t(t(test_1)) %*% t(test_2))
test_mat <- test_mat + matrix(0.2*runif(length(test_1)*length(test_2)),
                              nrow = length(test_1), ncol = length(test_2))


loc_max <- find_local_max(x = test_mat)

# Test with missing values.
test_mat[7, 7] <- NA # Not a local max - no change.
test_mat[4, 5] <- NA # Local max - replaced with neighbor.

loc_max <- find_local_max(x = test_mat)



# Test for 1-dimensional problem.
test_1_array <- as.array(cbind(test_1 - 1, test_1 + 0.2*runif(length(test_1)), test_1 - 1))

loc_max <- find_local_max(x = as.array(test_1_array))
test_1_array[, 2]


test_2_array <- as.array(cbind(test_2 - 1, test_2 + 0.2*runif(length(test_2)), test_2 - 1))

loc_max <- find_local_max(x = as.array(test_2_array))
test_2_array[, 2]



#--------------------------------------------------
# Make a simple version that will work.
#--------------------------------------------------

set.seed(42)

test_1 <- c(1,1,2,3,3,2,1,1,5,6,7,7,4,2) + 0.2*runif(length(test_1))
test_2 <- c(1,1,2,2,3,3,4,3,2,2) + 0.2*runif(length(test_2))
test_mat <- as.array(t(t(test_1)) %*% t(test_2))
test_mat <- test_mat + matrix(0.2*runif(length(test_1)*length(test_2)),
                              nrow = length(test_1), ncol = length(test_2))

x <- test_mat

find_local_max_simple(test_mat)

find_local_max_simple <- function(x) {

  # Determine size of input.
  nrows <- dim(x)[1]
  ncols <- dim(x)[2]

  # Remove corner cases.
  if (nrows == 1 & ncols == 1) {
    return(list(row = 1, col = 1, value = x))
  } else if (nrows == 1) {
    x <- t(x)
    nrows <- ncols
    ncols <- 1
  }

  # Initialize matrix of local max indicators.
  loc_max_ind <- matrix(TRUE, nrow = nrows, ncol = ncols)

  # Proceed by ruling out locally dominated points.

  # Check above.
  # x_check <- rbind(matrix( - Inf, nrow = 1, ncol = ncols),
  #                  x[ - nrows, ])
  x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
  x_check[ - 1, ] <- x[ - nrows, ]
  loc_max_ind <- loc_max_ind & (x > x_check)

  # Check below.
  # x_check <- rbind(x[ - 1, ],
  #                  matrix( - Inf, nrow = 1, ncol = ncols))
  x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
  x_check[ - nrows, ] <- x[ - 1, ]
  loc_max_ind <- loc_max_ind & (x > x_check)


  # Check sides and diagonals only if 2-dimensional matrix.
  if (ncols > 1) {

    # Check left.
    # x_check <- cbind(matrix( - Inf, nrow = nrows, ncol = 1),
    #                  x[ , - ncols])
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ , - 1] <- x[ , - ncols]
    loc_max_ind <- loc_max_ind & (x > x_check)

    # Check right.
    # x_check <- cbind(x[ , - 1],
    #                  matrix( - Inf, nrow = nrows, ncol = 1))
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ , - ncols] <- x[ , - 1]
    loc_max_ind <- loc_max_ind & (x > x_check)


    # Check top left.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ - 1, - 1] <- x[ - nrows, - ncols]
    loc_max_ind <- loc_max_ind & (x > x_check)

    # Check top right.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ - 1, - ncols] <- x[ - nrows, - 1]
    loc_max_ind <- loc_max_ind & (x > x_check)


    # Check bottom left.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ - nrows, - 1] <- x[ - 1, - ncols]
    loc_max_ind <- loc_max_ind & (x > x_check)

    # Check bottom right.
    x_check <- matrix( - Inf, nrow = nrows, ncol = ncols)
    x_check[ - nrows, - ncols] <- x[ - 1, - 1]
    loc_max_ind <- loc_max_ind & (x > x_check)

  }


  # Assemble output into a list.
  loc_max_out <- list(
    row = matrix(rep(seq(nrows), ncols),
                 nrow = nrows,
                 ncol = ncols)[loc_max_ind],
    col = matrix(rep(seq(ncols), nrows),
                 nrow = nrows,
                 ncol = ncols, byrow = TRUE)[loc_max_ind],
    value = x[loc_max_ind])

  return(loc_max_out)
}


find_local_max_simple(test_mat)
find_local_max_simple(as.matrix(test_1))
find_local_max_simple(as.matrix(test_2))


# First version with existing packages:

# find_local_max <- function(x) {
#
#   ## Convert it to a raster object
#   print('Using raster function next')
#   r <- raster::raster(x)
#   print('Used raster function')
#
#   print('Using extent function next')
#   # extent(r) <- extent(c(0, nrow(x), 0, ncol(x)) + 0.5)
#   raster::extent(r) <- raster::extent(c(0, ncol(x), 0, nrow(x)) + 0.5)
#   print('Used extent function')
#
#   ## Find the maximum value within the 9-cell neighborhood of each cell
#   # f <- function(X) max(X, na.rm = TRUE)
#   # ww <- matrix(1, nrow = 3, ncol = 3) ## Weight matrix for cells in moving window
#   print('Using focal function next')
#   localmax <- raster::focal(r,
#                             fun = function(X) max(X, na.rm = TRUE),
#                             w = matrix(1, nrow = 3, ncol = 3),
#                             pad = TRUE, padValue = NA)
#   print('Used focal function')
#
#   ## Does each cell have the maximum value in its neighborhood?
#   r2 <- r == localmax
#
#   # Which(r2 == 1)
#
#   ## Get x-y coordinates of those cells that are local maxima
#   print('Using xyFromCell function next')
#   maxXY <- raster::xyFromCell(r2, raster::Which(r2 == 1, cells = TRUE))
#   print('Used xyFromCell function')
#
#
#   # Obtain the proper y-coordinates.
#   # Notice the switched coordinates.
#   loc_max_out <- list(
#     # y_ind = ncol(x) - maxXY[, 'y'] + 1,
#     row = nrow(x) - maxXY[, 'y'] + 1,
#     col = maxXY[, 'x'],
#     value = rep(NA, nrow(maxXY))
#   )
#
#   # Append value of objective at local maxima.
#   for (i in 1:nrow(maxXY)) {
#     loc_max_out$value[i] <- x[loc_max_out$row[i], loc_max_out$col[i]]
#     # loc_max_out$z_ind[i] <- x[loc_max_out$y_ind[i], loc_max_out$x_ind[i]]
#   }
#
#
#   return(loc_max_out)
# }


