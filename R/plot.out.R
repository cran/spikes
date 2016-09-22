plot.out <-
function(x, main = NULL, ...){
  d <- function(z) density(z, from = 0, to = 1, cut = TRUE, n = x$grid, bw = x$bw, kernel = "gaussian")
  ydens <- d(x$data$v/x$data$t)
  plot(NULL, ylim = c(0, max(c(x$ymax, ydens$y)) + 0.5), xlim = c(0, 1), xaxt = "n", ylab = "Density", xlab = "Vote-share", cex.lab = 1.2, main = main)
  w <- which(ydens$y >= x$ymax)
  w <- w[x$w[w]!=0]
  abline(v = x$x[w], lwd = 6*x$w[w]/max(x$w[w]), col = "grey90")
  lines(x$x, x$ymax, type = "l", lwd = 2, lty = 1, col = "black")
  lines(x$x, ydens$y, col = "grey30", lwd = 3.5, lty = 1.5)
  legend("topleft", paste("F", roundr(mean(x$fraud), 2), sep = " = "), cex = 1.1, bty = "n")
  axis(1, at = seq(0, 1, by = 0.1))
}
