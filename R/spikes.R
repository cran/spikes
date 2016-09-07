spikes <-
function(data, resamples = 1000, bw = 0.0001, grid = 1001, out = NULL){

  if (is.null(out)){
  cat("Estimating latent density of turnout", "\n")
  out$t <- est.density(y = data$t, n = data$N)
  dev.off()
  cat("\n")
  cat("Estimating latent density of support", "\n")
  out$v <- est.density(y = data$v, n = data$t)
  }

  cat("\n")
  cat("Resampling", "\n")

  output  <- fraud(out, data = data, resamples = resamples, grid = grid, bw = bw)

  out <- list(fraud = output$fraud, x = output$x, data = data, ymax = output$ymax, w = output$w, bw = bw, grid = grid, bins = output$bins, resamples = resamples, out = out, call = match.call())
  class(out) <- "out"
  cat(paste("Estimated percentage of fraudulent precincts: ", round(out$fraud, 2), sep = " "), "\n")
  out

}
