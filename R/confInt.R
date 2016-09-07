confInt <-
function(object, boots = 100){

  # Initial values from MLE
  theta <- NULL
  theta$t <- object$out$t
  theta$v <- object$out$v

  pb <- txtProgressBar(min = 0, max = boots, style = 3)

  y <- rep(NA, boots)
  ymax <- w <- matrix(NA, nrow = boots, ncol = object$grid)

  for (i in 1:boots){
    theta$t <- bupdate(y = object$data$t, n = object$data$N, theta$t)
    theta$v <- bupdate(y = object$data$v, n = object$data$t, theta$v)
    f <- fraud(theta, object$data, resamples = object$resamples, grid = object$grid, bw = object$bw)
    y[i] <- f$fraud
    ymax[i, ] <- f$ymax
    w[i, ] <- f$w
    setTxtProgressBar(pb, i)
  }

  ymax <- apply(ymax, 2, mean)
  w    <- apply(w, 2, mean)
  out <- list(fraud = y, ymax = ymax, x = object$x, bins = object$bins, bw = object$bw, grid = object$grid, data = object$data, w = w, out = list(t = object$out$t, v = object$out$v))
  class(out) <- "out"
  return(out)

}
