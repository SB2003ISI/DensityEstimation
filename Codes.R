rangexp=function(u,n=10){
  n*exp(-u)*((1-exp(-u))^(n-1))
}

f <- function(x) {
  rangexp(x, n = n_val)
}
#---Sample generator---
rangesamp=function(m,n=10){
  mat=matrix(rexp(m*(n+1),rate=1),nr=m)
  return(apply(mat,MARGIN = 1,FUN = function(x){return(max(x)-min(x))}))
}

x0=seq(from=0.5,to=9,by=0.5)  
den.exact=rangexp(x0)

#---Kernel Density Estimation---
x=rangesamp(m=500)
bw_ls=bw.ucv(x)
dens=density(x,bw=bw_ls,kernel = "gaussian")
plot(dens)
curve(rangexp,from=0,to=10,col="red",add=T)
density_at_new <- function(x_new, data, h, kernel=dnorm) {
  sapply(x_new, function(x0) {
    mean(kernel((x0 - data) / h)) / h
  })
}
den.kern=density_at_new(x_new=x0,data=x,h=bw_ls)
den.exact-den.kern

#----Edgeworth Expansion---
library(Deriv)
library(orthopolynom)
library(ggplot2)

# Define the Hermite function
hermite_eval <- function(order, x, type = c("probabilists", "physicists"), normalized = FALSE) {
  type <- match.arg(type)
  if (order < 0 || order != floor(order)) {
    stop("Order must be a nonnegative integer")
  }
  
  # Generate the polynomials up to the specified order.
  # The list will contain polynomials from order 0 to 'order'.
  if (type == "probabilists") {
    polys <- hermite.he.polynomials(n = order, normalized = normalized)
  } else {
    polys <- hermite.h.polynomials(n = order, normalized = normalized)
  }
  
  # Evaluate each polynomial at the point(s) x.
  # The polynomial of the desired order is at index order + 1.
  eval_list <- polynomial.values(polynomials = polys, x = x)
  result <- eval_list[[order + 1]]
  return(result)
}

#-A function for calculating derivatives--
cgf=function(t,n=10){
  return(sum(sapply(1:n,FUN = function(j){-log(1-(t/j))})))
}

cgfk=function(t,k,n=10){
  return(factorial(k-1)*mean(sapply(1:n,FUN = function(j){1/((j-t)^k)})))
}


den.edge <- function(y, n = 10) {
  # Compute location and scale adjustments
  mun <- sum(sapply(1:n, function(j) { 1/j }))
  sigman <- sqrt(sum(sapply(1:n, function(j) { 1/(j^2) })))
  
  # Standardize input y
  z <- (y - mun) / sigman
  
  # Cumulant generating function helper
  #K
  cgfk <- function(t, k) {
    factorial(k - 1) * mean(sapply(1:n, function(j) { (j - t)^(-k) }))
  }
  
  # Standardized cumulant estimator rho
  rho <- function(r) {
    cgfk(t = 0, k = r) / (cgfk(t = 0, k = 2)^(r / 2))
  }
  
  # Evaluate probabilists' Hermite polynomials at z
  H3 <- hermite_eval(order = 3, x = z, type = "probabilists")
  H4 <- hermite_eval(order = 4, x = z, type = "probabilists")
  H6 <- hermite_eval(order = 6, x = z, type = "probabilists")
  
  # Compute the density with Edgeworth expansion correction
  result <- dnorm(z) * (1 + (rho(3) * H3 / (6 * sqrt(n))) +
                          (rho(4) * H4 / (24 * n)) +
                          ((rho(3)^2) * H6 / (72 * n))) / sigman
  
  # Replace negative values with 0
  result <- pmax(result, 0)
  
  return(result)
}

den.edge(x0)
den.exact-den.edge(x0)

#---Visualization---
n_val=10
gridp=seq(from=0,to=9,length.out=10000)
plot(y=den.edge(gridp,n=n_val),x=gridp,type="l")
f <- function(x) {
  rangexp(x, n = n_val)
}
curve(f, from = 0, to = 10, col = "red", add = TRUE)


#---Approximation by samples---
#--Approximate the cumulants---
sample_cumulants <- function(x) {
  # Compute the mean (first cumulant)
  mu <- mean(x)
  
  # Compute central moments
  mu2 <- mean((x - mu)^2)         # Second central moment
  mu3 <- mean((x - mu)^3)         # Third central moment
  mu4 <- mean((x - mu)^4)         # Fourth central moment
  
  # Compute cumulants
  kappa1 <- mu
  kappa2 <- mu2
  kappa3 <- mu3
  kappa4 <- mu4 - 3 * mu2^2
  
  cumulants <- c(kappa1, kappa2, kappa3, kappa4)
  names(cumulants) <- c("kappa1", "kappa2", "kappa3", "kappa4")
  
  return(cumulants)
}

den.edge2=function(y,sam,n=10){
  samcum=sample_cumulants(sam)
  mu=samcum[1]
  sigma=sqrt(samcum[2])
  z=(y-mu)/sigma
  rho=function(r){
    samcum[r]/(sigma^r)
  }
  # Evaluate probabilists' Hermite polynomials at z
  H3 <- hermite_eval(order = 3, x = z, type = "probabilists")
  H4 <- hermite_eval(order = 4, x = z, type = "probabilists")
  H6 <- hermite_eval(order = 6, x = z, type = "probabilists")
  
  # Compute the density with Edgeworth expansion correction
  result <- dnorm(z) * (1 + (rho(3) * H3 / (6 * sqrt(n))) +
                          (rho(4) * H4 / (24 * n)) +
                          ((rho(3)^2) * H6 / (72 * n))) / sigma
  
  # Replace negative values with 0
  result <- pmax(result, 0)
  }

#--Visulaization--
n_val=10
sam=rangesamp(m=500,n=n_val)
gridp=seq(from=0,to=9,length.out=10000)
plot(y=den.edge2(gridp,sam = sam,n=n_val),x=gridp,type="l")
f <- function(x) {
  rangexp(x, n = n_val)
}
curve(f, from = 0, to = 10, col = "red", add = TRUE)


#------Saddle Point Approximation----

den.saddle1 <- function(x0, n = 10) {
  # Cumulant generating function
  cgf <- function(t) {
    sum(-log(1 - t / (1:n)))
  }
  
  # Second derivative of CGF
  cgfk2 <- function(t) {
    sum((1:n - t)^(-2))
  }
  
  # Compute density for each x0
  result <- sapply(x0, function(x) {
    # Saddlepoint equation: K'(t) - x = 0
    f <- function(t) sum(1 / (1:n - t)) - x
    
    # Find root in (-Inf, 1)
    t_0 <- tryCatch(
      uniroot(f, c(-1e6, 1 - 1e-6))$root,
      error = function(e) NA
    )
    
    if (is.na(t_0)) return(NA)
    
    # Saddlepoint approximation formula
    exp(cgf(t_0) - t_0 * x) / sqrt(2 * pi * cgfk2(t_0))
  })
  
  return(result)
}

#--Visulaization--
n_val=10
gridp=seq(from=0,to=9,length.out=10000)
plot(y=den.saddle1(gridp,n=n_val),x=gridp,type="l")
f <- function(x) {
  rangexp(x, n = n_val)
}
curve(f, from = 0, to = 10, col = "red", add = TRUE)

#--------------
#---Saddle Point approximation 2---
den.saddle2 <- function(x0, n = 10) {
  # Cumulant generating function
  #nK
  cgf <- function(t) {
    sum(-log(1 - t / (1:n)))
  }
  
  #K
  cgfk <- function(t, k) {
    factorial(k - 1) * mean(sapply(1:n, function(j) { (j - t)^(-k) }))
  }
  
  
  rho <- function(t,r) {
    n*cgfk(t = t, k = r) / ((cgfk(t = t, k = 2)*n)^(r / 2))
  }
  
  # Compute density for each x0
  result <- sapply(x0, function(x) {
    # Saddlepoint equation: K'(t) - x = 0
    f <- function(t) sum(1 / (1:n - t)) - x
    
    # Find root in (-Inf, 1)
    t_0 <- tryCatch(
      uniroot(f, c(-1e6, 1 - 1e-6))$root,
      error = function(e) NA
    )
    
    if (is.na(t_0)) return(NA)
    
    # Saddlepoint approximation formula
    (exp(cgf(t_0) - t_0 * x) / sqrt(2 * pi *n* cgfk(t=t_0,k=2)))*(1+((3*rho(t=t_0,r=4)-5*rho(t=t_0,r=3)^2)/(24*n)))
  })
  
  return(result)
}

#---Visualization---
n_val=100
gridp=seq(from=0,to=9,length.out=10000)
plot(y=den.saddle2(gridp,n=n_val),x=gridp,type="l")
f <- function(x) {
  rangexp(x, n = n_val)
}
curve(f, from = 0, to = 10, col = "red", add = TRUE)


#----Aggregation----

library(ggplot2)

density_estimation_gg <- function(x, method = "exact", data = NULL, n = 10, plotit = TRUE) {
  # x      : vector of points at which to evaluate the density.
  # method : one of "exact", "kernel", "edgeworth", "edgeworth2", "saddle1", or "saddle2".
  # data   : sample data required for kernel or sample-based edgeworth methods.
  # n      : parameter used in the functions (e.g., sample size for rangexp).
  # plotit : logical. If TRUE, the function returns a ggplot object of the estimated density.
  
  # Check that required functions are available:
  required_funcs <- c("rangexp", "density_at_new", "den.edge", "den.edge2", "den.saddle1", "den.saddle2")
  missing_funcs <- required_funcs[!sapply(required_funcs, exists)]
  if (length(missing_funcs) > 0) {
    stop("The following required functions are not available: ", paste(missing_funcs, collapse = ", "))
  }
  
  # Compute the density based on the method selected
  if (method == "exact") {
    dens_method <- rangexp(x, n = n)
  } else if (method == "kernel") {
    if (is.null(data)) {
      stop("For the kernel method, please supply sample data in the 'data' argument.")
    }
    bw_ls <- bw.ucv(data)
    dens_method <- density_at_new(x_new = x, data = data, h = bw_ls)
  } else if (method == "edgeworth") {
    dens_method <- den.edge(x, n = n)
  } else if (method == "edgeworth2") {
    if (is.null(data)) {
      stop("For the edgeworth2 method, please supply sample data in the 'data' argument.")
    }
    dens_method <- den.edge2(x, sam = data, n = n)
  } else if (method == "saddle1") {
    dens_method <- den.saddle1(x, n = n)
  } else if (method == "saddle2") {
    dens_method <- den.saddle2(x, n = n)
  } else {
    stop("Invalid method. Choose one of: exact, kernel, edgeworth, edgeworth2, saddle1, saddle2.")
  }
  
  # Create a data frame for the method's density
  df_method <- data.frame(x = x, density = dens_method, Method = method)
  
  # If method is not "exact", also compute the exact density for reference
  if (method != "exact") {
    df_exact <- data.frame(x = x, density = rangexp(x, n = n), Method = "exact")
    df_all <- rbind(df_method, df_exact)
  } else {
    df_all <- df_method
  }
  
  # If plotit is FALSE, return the data frame
  if (!plotit) {
    return(df_all)
  }
  
  # Define fill colors for the methods
  fill_colors <- c("exact" = "#D55E00", "kernel" = "#0072B2",
                   "edgeworth" = "#009E73", "edgeworth2" = "#E69F00",
                   "saddle1" = "#56B4E9", "saddle2" = "#CC79A7")
  
  # Use only the colors needed for the plot
  fill_colors_used <- fill_colors[names(fill_colors) %in% unique(df_all$Method)]
  
  # Build the ggplot using geom_area for filled density curves
  p <- ggplot(df_all, aes(x = x, y = density, fill = Method)) +
    geom_area(alpha = 0.5, color = "black") +
    labs(title = paste("Density Estimation using", method, "method"),
         x = "x", y = "Density", fill = "Method") +
    theme_minimal(base_size = 14) +
    scale_fill_manual(values = fill_colors_used)
  
  return(p)
}


gridp <- seq(0, 10, length.out = 1000)
density_estimation_gg(x = gridp, method = "exact", n = 10, plotit = TRUE)

# Assuming 'rangesamp' is defined and returns sample data:
sample_data <- rangesamp(m = 500, n = 10)
gridp <- seq(0, 10, length.out = 1000)
density_estimation_gg(x = gridp, method = "kernel", data = sample_data, n = 10, plotit = TRUE)
density_estimation_gg(x = gridp, method = c("edgeworth"), data = NULL, n = 10, plotit = TRUE)
density_estimation_gg(x = gridp, method = c("saddle1"), data = NULL, n = 10, plotit = TRUE)
density_estimation_gg(x = gridp, method = c("saddle2"), data = NULL, n = 10, plotit = TRUE)







