a <-
function(x){(1+exp(x))/exp(x)}

init <-
  function(L){
    list(alpha = rgamma(L, runif(L, .1, 5), 1), beta = rgamma(L, runif(L, .1, 5), 1), w = rep(1/L, L))
  }


range01 <-
  function(x){
    if(length(x)==1) return(x)
    (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
  }


samp.mat <-
  function(d) {
    x <- runif(nrow(d))
    cumul.w <- d %*% upper.tri(diag(ncol(d)), diag = TRUE) /
      rowSums(d)
    i <- rowSums(x > cumul.w) + 1L
    (1:ncol(d))[i]
  }

y.v <-
  function(y){
    y=y[!is.na(y)]
    if (length(y) == 0) return(1)
    A=-a(y)*diag(length(y))
    A=cbind(0,A)
    A=rbind(A,1)
    diag(A)=1
    v=solve(A)[,length(y)+1]
    v
  }

roundr <-
  function(x, n = 2) {
    f <- function(x, n) sprintf(paste("%.", n, "f", sep=""), round(x, n))
    if (is.null(dim(x))) out <- f(x, n)
    if (!is.null(dim(x))) out <- apply(x, 2, f, n = n)
    out
  }

dmbetabinom <-
  function(y, n, theta){
    sapply(1:length(y), function(s) sapply(1:length(theta$w), function(k) dbetabinom(y[s], size = n[s], shape1 = theta$alpha[k], shape2 = theta$beta[k]))%*%theta$w)
  }

dmbeta <-
  function(x, theta){
    f <- function(x){
      if(x==0 | x==1) return(0)
      sapply(1:length(theta$w), function(k) dbeta(x, shape1 = theta$alpha[k], shape2 = theta$beta[k]))%*%theta$w
    }
    sapply(x, f)
  }

rmnorm <- function (n = 1, mean = rep(0, d), varcov, sqrt = NULL)
{
  sqrt.varcov <- if (is.null(sqrt))
    chol(varcov)
  else sqrt
  d <- if (is.matrix(sqrt.varcov))
    ncol(sqrt.varcov)
  else 1
  mean <- outer(rep(1, n), as.vector(matrix(mean, d)))
  drop(mean + t(matrix(rnorm(n * d), d, n)) %*% sqrt.varcov)
}

resample <- function(out, data){

    N <- nrow(data)

    # Turnout
    K <- length(out$t$alpha)
    W <- sapply(1:K, function(m) dbetabinom(data$t, size = data$N, shape1 = out$t$alpha[m], shape2 = out$t$beta[m])*out$t$w[m])
    w <- samp.mat(W)
    pt <- rbeta(N, data$t + out$t$alpha[w], data$N - data$t + out$t$beta[w])
    t  <- rbinom(N, data$N, pt)

    # Voting
    K <- length(out$v$alpha)
    W <- sapply(1:K, function(m) dbetabinom(data$v, size = data$t, shape1 = out$v$alpha[m], shape2 = out$v$beta[m])*out$v$w[m])
    w <- samp.mat(W)
    pv <- rbeta(N, data$v + out$v$alpha[w], data$t - data$v + out$v$beta[w])
    v  <- rbinom(N, t, pv)
    y  <- v/t
    y[t==0] <- 0
    y
  }

em.update <-
  function(y, n, theta){

    L <- length(theta$w)

    # Clustering weights
    W <- sapply(1:L, function(k) theta$w[k]*dbetabinom(y, size = n, shape1 = theta$alpha[k], shape2 = theta$beta[k]))
    W <- W/apply(W, 1, sum)
    theta$w <- apply(W, 2, sum)/length(y)

    # Weighted likelihood function
    llik <- function(pars){
      alpha <- pars[1:L]
      beta  <- pars[(L+1):(2*L)]
      sum(W*sapply(1:L, function(k) dbetabinom(x = y, size = n, shape1 = exp(alpha[k]), shape2 = exp(beta[k]), log = TRUE)))
    }

    out <- optim(log(c(theta$alpha, theta$beta)), llik, method = "L-BFGS-B", control=list(fnscale=-1), hessian = TRUE, lower = rep(-3, L), upper = rep(6, L))

    theta$alpha   <- exp(out$par[1:L])
    theta$beta    <- exp(out$par[(L+1):(2*L)])
    theta$hessian <- out$hessian
    theta
  }


fraud <- function(out, data, resamples, grid, bw){

    d <- function(x) density(x, from = 0, to = 1, cut = TRUE, n = grid, bw = bw, kernel = "gaussian")
    ydens <- d(data$v/data$t)

    x <- ydens$x # points at which density is estimated
    bins <- seq(-x[2]/2, 1 + x[2]/2, by = x[2] - x[1]) # bins

    W <- function(x) c(prop.table(table(cut(x, bins, levels = TRUE, include.lowest = TRUE))))

    rsamp <- function(t){
      y <- resample(out, data = data)
      list(ydens = d(y)$y, w = W(y))
    }

    outsamp <- lapply(1:resamples, rsamp)
    w <- apply(do.call("rbind", lapply(outsamp, function(x) x$w)), 2, mean)
    w <- W(data$v/data$t) - w
    w <- ifelse(w > 0, w, 0)
    ymax <- apply(do.call("rbind", lapply(outsamp, function(x) x$ydens)), 2, max)
    fraud <- 100*sum(w[ymax < ydens$y])

    list(fraud = fraud, ymax = ymax, x = x, bins = bins, resamples = resample, grid = grid, bw = bw, w = w)
  }

bupdate <-
  function(y, n, theta){

    rdirichlet <- function (n, alpha){
      l <- length(alpha)
      x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
      sm <- x %*% rep(1, l)
      return(x/as.vector(sm))
    }

    lpost <- function(theta, k) {
      sum(dbetabinom(y, size = n, shape1 = theta$alpha[k], shape2 = theta$beta[k], log = TRUE))
    }

    update.k <- function(theta){
      L <- length(theta$w)
      d <- sapply(1:L, function(m) dbetabinom(y, size = n, shape1 = theta$alpha[m], shape2 = theta$beta[m])*theta$w[m])
      d <- d/apply(d, 1, sum)
      k <- L - 1
      while(length(unique(k))!=L) k <- samp.mat(d)
      k
    }

    V <- solve(-theta$hessian)

    update.theta <- function(theta, k){
      new <- rmnorm(1, log(c(theta$alpha, theta$beta)), varcov = V)
      new <- matrix(new, ncol = 2)
      cand <- list(alpha = exp(new[,1]), beta = exp(new[,2]))
      alpha <- min(0, lpost(cand, k) - lpost(theta, k))
      if(runif(1) < exp(alpha)) {
        theta$alpha <- exp(new[,1])
        theta$beta  <- exp(new[,2])
      }
      theta
    }

    update.w <- function(k){
      theta$w <- c(rdirichlet(1, c(table(factor(k)) + 1)))
      theta
    }

    k <- update.k(theta)
    theta <- update.theta(theta, k)
    theta <- update.w(k)

    theta

  }


est.density <-
  function(y, n, Lmin = 1, Lmax = 5){

    lambda <- 0 # Correlation between consequtive estimates within cycle L
    Lambda <- 0 # Correlation between cycles L, L+1, ...
    L <- Lmin

    fn <- function(pars, L){
      theta$alpha <- exp(pars[1:L])
      theta$beta  <- exp(pars[(L+1):(2*L)])
      if(L == 1) theta$w <- 1
      if(L > 1) theta$w <- y.v(pars[(2*L+1):(3*L-1)])
      mean((dmbeta(x, theta) - ydens$y)^2)
    }

    optfn <- function(L){
      Li <- ifelse(L==1, 2, 3*L - 1)
      out <- optim(rnorm(Li), fn, method = "L-BFGS-B", L = L, lower = rep(-3, L), upper = rep(6, L))
      theta$alpha <- exp(out$par[1:L])
      theta$beta  <- exp(out$par[(L+1):(2*L)])
      if(L==1)  theta$w <- 1
      if(L > 1) theta$w <- y.v(out$par[(2*L+1):(3*L-1)])
      theta
    }


    while(L <= Lmax & Lambda < 0.99){

      theta <- init(L)
      lambda <- 0
      ydens <- density(y/n, bw = 0.001, n = 1001, from = 0.01, to = .99, cut = TRUE)
      x <- ydens$x
      plot(ydens, main = paste("Fitting with L = ", L, sep = ""))
      theta <- optfn(L) #Initialize parameters
      lines(x, dmbeta(x, theta), col = "blue", lty = 3)

      # Continue until the correlation between the consequtive iterates is below delta
      delta <- 0
      while (delta < 0.99) {
        old <- theta
        theta <- em.update(y = y, n = n, theta)
        delta <- cor(dmbeta(x, theta), dmbeta(x, old))
        lines(x, dmbeta(x, theta), col = "green")
      }

      if(L > 1) Lambda <- cor(dmbeta(x, theta), dmbeta(x, Old))
      Old <- theta
      L   <- L + 1
    }

    return(theta)
  }


# E-M algorithm:
em.update <- function(y, n, theta){

  L <- length(theta$w)

  # Clustering weights
  W <- sapply(1:L, function(k) theta$w[k]*dbetabinom(y, size = n, shape1 = theta$alpha[k], shape2 = theta$beta[k]))
  W <- W/apply(W, 1, sum)
  theta$w <- apply(W, 2, sum)/length(y)

  # Weighted likelihood function
  llik <- function(pars){
    alpha <- pars[1:L]
    beta  <- pars[(L+1):(2*L)]
    sum(W*sapply(1:L, function(k) dbetabinom(x = y, size = n, shape1 = exp(alpha[k]), shape2 = exp(beta[k]), log = TRUE)))
  }

  out <- optim(log(c(theta$alpha, theta$beta)), llik, method = "L-BFGS-B", control=list(fnscale=-1), hessian = TRUE, lower = rep(-3, L), upper = rep(6, L))

  theta$alpha   <- exp(out$par[1:L])
  theta$beta    <- exp(out$par[(L+1):(2*L)])
  theta$hessian <- out$hessian
  theta
}
