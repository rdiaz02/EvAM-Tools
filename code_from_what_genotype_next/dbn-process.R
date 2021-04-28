library(OncoBN)


do_DBN <- function(data){
  
  out <- fitCPN(data)
  out 
}

# 
# N <- 50
# na <- N
# nb <- N + round( 10 * runif(1))
# nab <- .6 * N + round( 10 * runif(1))
# n00 <- round( 10 * runif(1))
# dB <- matrix(
#   c(
#     rep(c(1, 0), na) 
#     , rep(c(0, 1), nb)
#     , rep(c(1, 1), nab)
#     , rep(c(0, 0), n00)
#   ), ncol = 2, byrow = TRUE
# )
# colnames(dB) <- LETTERS[1:2]
# 
# z <- fitCPN(dB, model="DBN", algorithm="DP")

