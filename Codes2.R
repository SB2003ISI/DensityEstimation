f <- function(x) {
  sum(1 / (1:10 - x)) - 0.5
}

result <- uniroot(f, c(-20, -10))
root <- result$root
print(root)


uniroot(f, c(1.001, 1.999))

#------
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
round(den.saddle1(x0=seq(from=0.5,to=9,by=0.5)),5)











