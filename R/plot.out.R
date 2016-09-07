plot.out <-
function(x, ...){
  d <- function(z) density(z, from = 0, to = 1, cut = TRUE, n = x$grid, bw = x$bw, kernel = "gaussian")
  ydens <- d(x$data$v/x$data$t)
  plot(NULL, ylim = c(0, max(c(x$ymax, ydens$y)) + 0.5), xlim = c(0, 1), xaxt = "n", ylab = "Density", xlab = "Vote-share", cex.lab = 1.2)
  w <- which(ydens$y >= x$ymax)
  w <- w[x$w[w]!=0]
  abline(v = x$x[w], lwd = 5*x$w[w]/max(x$w[w]), col = "grey90")
  lines(x$x, x$ymax, type = "l", lwd = 2, lty = 1, col = "black")
  lines(x$x, ydens$y, col = "dark red", lwd = 1.5, lty = 1)
  legend(-.05, max(x$ymax) + 0.25, paste("F", roundr(mean(x$fraud), 2), sep = " = "), cex = 1.5, bty = "n")
  axis(1, at = seq(0, 1, by = 0.1))
  b <- coefficients(lm(ydens$y ~ dmbeta(x$x,x$out$v)))
  lines(x$x, b[1] + b[2]*dmbeta(x$x, x$out$v), col = "white", lwd = 2)
}
