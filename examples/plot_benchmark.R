
x <- readRDS(sprintf("b1/times.rds"))

plt_data <- sapply(x, function(x){
    tmp_data <- matrix(unlist(x), ncol = 3, byrow = TRUE)
    return(c(
        apply(tmp_data, 2, mean)
        , apply(tmp_data, 2, function(x)sd(x))
    ))
})

titles <- c("Theta 2 trm_SM", "SIM 10 000", "SIM 100 000")
png("b1/benchmark.png")
par(mfrow=c(3,1), mai = rep(0.5, 4))
for(i in 1:3){
    plot(c(3, 5, 7, 9, 11, 13)
        , plt_data[i,]
        , xlab = "#Genes"
        , ylab = "Time (s)"
        , main = titles[i])

    arrows(c(3, 5, 7, 9, 11, 13)
        , y0 = plt_data[i,] - plt_data[i+3,]
        , y1 = plt_data[i,] + plt_data[i+3,]
        , code=3, length=0.02, angle = 90
        )
}

dev.off()