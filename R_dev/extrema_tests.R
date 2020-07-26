

# # Update R version.
# install.packages('installr')
# library(installr)
#
# updateR()


# Testing functions to locate local optima.

install.packages('ggpmisc')
library(ggpmisc)

test_1 <- c(1,1,2,3,3,2,1,1,5,6,7,7,4,2)

peak_locs <- ggpmisc:::find_peaks(x = test_1,
                                  ignore_threshold = 0,
                                  span = 3,
                                  # strict = TRUE,
                                  strict = FALSE
                                  )

seq(1, length(test_1))[peak_locs]
test_1[peak_locs]



# Try a 2D example

test_mat <- as.array(t(t(test_1)) %*% t(test_1))


peak_locs <- ggpmisc:::find_peaks(x = test_mat,
                                  ignore_threshold = 0,
                                  span = 3,
                                  # strict = TRUE,
                                  strict = FALSE)

seq(1, length(test_mat))[peak_locs]
test_mat[peak_locs]




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
  extent(r) <- extent(c(0, nrow(x), 0, ncol(x)) + 0.5)

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

  loc_max_out <- list(
    x_ind = maxXY[, 'x'],
    y_ind = ncol(x) - maxXY[, 'y'] + 1,
    z_ind = rep(NA, nrow(maxXY))
  )

  # Append value of objective at local maxima.
  for (i in 1:nrow(maxXY)) {
    loc_max_out$z_ind[i] <- x[loc_max_out$x_ind[i], loc_max_out$y_ind[i]]
  }


  return(loc_max_out)
}


test_1 <- c(1,1,2,3,3,2,1,1,5,6,7,7,4,2)
test_mat <- as.array(t(t(test_1)) %*% t(test_1))

loc_max <- find_local_max(x = test_mat)

# test_mat[loc_max]


