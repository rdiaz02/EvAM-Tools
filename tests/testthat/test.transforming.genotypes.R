source("../../code_from_what_genotype_next/utils.R")

int_state1 <- 0
int_state2 <- 15
int_state3 <- 9
int_state4 <- -4
int_state5 <- "abcd"

binary_state1 <- c(0)
binary_state2 <- c(1, 1, 1, 1)
binary_state3 <- c(1, 0, 0, 1)
binary_state4 <- c(1, 1, 1, -3)
binary_state5 <- c()

str_state1 <- "WT"
str_state2 <- "A, B, C, D"
str_state3 <- "A, D"
str_state4 <- "B, A, D, C"
str_state5 <- 42

test_that("Transforming 2binary works",{
    expect_equal(int2binary(int_state1), binary_state1)
    expect_equal(int2binary(int_state2), binary_state2)
    expect_equal(int2binary(int_state3), binary_state3)
    
    expect_equal(int2binary(int_state3, n = 5), c(0, 1, 0, 0, 1))
    expect_equal(int2binary(int_state3, n = -5), binary_state3)
    expect_error(int2binary(int_state4, n = -4))
    expect_error(int2binary("abc", n = -4))

    expect_equal(str2binary(str_state1), binary_state1)
    expect_equal(str2binary(str_state2), binary_state2)
    expect_equal(str2binary(str_state3), binary_state3)
    expect_equal(str2binary(str_state4), binary_state2)
    expect_equal(str2binary("ABCD", sep = ""), binary_state2)
    expect_equal(str2binary("A-B-C-D", sep = "-"), binary_state2)
    expect_equal(str2binary("ROOT", sep = "-", wt = "ROOT"), binary_state1)
    expect_equal(str2binary("A-B-C-D", sep = "-", n = 5), c(0, binary_state2))
    expect_error(int2binary(int_state5))
    expect_error(str2binary(str_state5))
})

test_that("Transforming to string works", {
    expect_equal(int2str(int_state1), str_state1)
    expect_equal(int2str(int_state2), str_state2)
    expect_equal(int2str(int_state3), str_state3)

    expect_equal(binary2str(binary_state1), str_state1)
    expect_equal(binary2str(binary_state2), str_state2)
    expect_equal(binary2str(binary_state3), str_state3)
    expect_error(binary2str(binary_state4))
    expect_error(binary2str(binary_state5))
    expect_error(int2str(int_state5))

})

test_that("Transforming to integer works", {
    expect_equal(str2int(str_state1), int_state1)
    expect_equal(str2int(str_state2), int_state2)
    expect_equal(str2int(str_state3), int_state3)

    expect_equal(binary2int(binary_state1), int_state1)
    expect_equal(binary2int(binary_state2), int_state2)
    expect_equal(binary2int(binary_state3), int_state3)
    expect_error(binary2int(binary_state4))
    expect_error(binary2int(binary_state5))
    expect_error(str2int(str_state5))
})
