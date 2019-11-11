# confidence intervals for GEV and MEV distributions
# 
# hs, 10.2019
conf.int <- function(x, alpha = 0.05, return.periods = c(2,10,30,50,100), dist = c("GEV","MEV"), method = "boot", R = 502, verbose = FALSE){
     
     if(!inherits(x, "data.frame"))
          stop("conf.int: x must be data.frame")
     
     dist <- match.arg(dist)
     
     if(verbose & method == "boot"){
          cat("\n", "Using Bootstrap Method.\n")
     }
     
     switch(dist,
          "MEV" = {
               # start from fit with w and C
               # x must be data.frame with columns W,C,n
               if(!(all(colnames(x) %in% c("W","C","n"))))
                    stop("conf.int: x must be data.frame with columns 'W','C','n'")
               
               w=x$w #shape
               C=x$C #scale
               n=x$n # number of wet days
               l=length(as.numeric(na.omit(as.vector(x$data)))) # length of data series
               theta.hat <- c(shape=w,scale=C,n=n)
               
               
               # draw R samples of length l with w and C
               if (verbose) cat("\n", "Simulating data from fitted model.  Size = ", R, "\n")
               Z <- rweibull(n = l * R, shape = w, scale = C) 
               Z <- matrix(Z, l, R)
               if (verbose) cat("\n", "Simulated data found.\n")
               
               
               # fit MEVD to simulated samples
               if (verbose) cat("\n", "Fitting model to simulated data sets (this may take a while!).")
               bfun1 <- function(z, n){
                    fit <- fmev(z, n, threshold = 0, type ="simple", method="pwm") #
                    return(c(shape=fit$w,scale=fit$C,n=fit$n))
               }
               pars <- apply(Z, 2, bfun1, n = n)
               shape <- pars["shape",]
               scale <- pars["scale",]
               n <- pars["n",]
               
               
               # compute return levels from R w and C parameters
               th <- rbind(shape, scale, n)
               th.est <- theta.hat
               rlfun <- function(theta, q) rlmev(q = q, w = theta[1], c = theta[2], n = theta[3])
               sam <- apply(th, 2, rlfun, q = return.periods)
               rownames(sam) <- paste0(return.periods, "-year")
               theta.hat <- rlmev(q = return.periods, w = th.est[1], c = th.est[2], n = th.est[3])
               
          },
          "GEV" = {
               
               # start from fit with loc, scale and shape
               # x must be data.frame with columns loc,scale,shape
               if(!(all(colnames(x) %in% c("loc","scale","shape"))))
                    stop("conf.int: x must be data.frame with columns 'loc','scale','shape'")
               
               loc=x$loc
               scale=x$scale
               shape=x$shape
               l=length(return.periods) 
               theta.hat <- c(loc=loc,scale=scale,shape=shape)
               
               
               # draw R samples of length l with loc, scale and shape
               if (verbose) cat("\n", "Simulating data from fitted model.  Size = ", R, "\n")
               Z <- rgev(n = l * R, loc = loc, scale = scale, shape = shape) 
               Z <- matrix(Z, l, R)
               if (verbose) cat("\n", "Simulated data found.\n")
               
               
               # fit GEV to simulated samples
               if (verbose) cat("\n", "Fitting model to simulated data sets (this may take a while!).")
               bfun <- function(z){
                    return( gevmle(z) )
               }
               pars <- apply(Z, 2, bfun)
               loc <- pars["loc",]
               scale <- pars["scale",]
               shape <- pars["shape",]
               
               
               # compute return levels from R loc, scale and shape parameters
               th <- rbind(loc, scale, shape)
               th.est <- theta.hat
               rlfun <- function(theta, q) qgev(p = 1-1/q, loc = theta[1], scale = theta[2], shape = theta[3])
               sam <- apply(th, 2, rlfun, q = return.periods)
               rownames(sam) <- paste0(return.periods, "-year")
               theta.hat <- qgev(p= 1-1/return.periods, loc = th.est[1], scale = th.est[2], shape = th.est[3])

          }
     )
     
     # compute quantiles of simulated return levels
     out <- apply(sam, 1, quantile, probs = c(alpha/2,1 - alpha/2))
     out.names <- rownames(out)
     out <- rbind(out[1, ], theta.hat, out[2, ])
     rownames(out) <- c(out.names[1], "Estimate", out.names[2])
     colnames(out) <- rownames(sam)
     out <- t(out)
     if (verbose) cat("\n", "Finished fitting model to simulated data.\n")
     
     return(out)
     
}

