summary.out <-
function(object, ...){
  est <- roundr(c(mean(object$fraud), quantile(object$fraud, c(0.025, 0.975))), 2)
  c(est[1], paste("(", est[2], ", ", est[3], ")", sep = ""))
}
