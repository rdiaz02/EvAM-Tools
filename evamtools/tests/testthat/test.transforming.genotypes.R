# source("../../R/utils.R")

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
    expect_equal(int2binary(0, n = 5), c(0, 0, 0, 0, 0))
    
    expect_equal(int2binary(int_state3, n = 5), c(1, 0, 0, 1, 0))
    expect_equal(int2binary(int_state3, n = -5), binary_state3)
    expect_error(int2binary(int_state4, n = -4, "Only integers >= 0 are valid"))
    expect_error(int2binary("abc", n = -4), "Invalid input")

    expect_equal(str2binary(str_state1), binary_state1)
    expect_equal(str2binary(str_state2), binary_state2)
    expect_equal(str2binary(str_state3), binary_state3)
    expect_equal(str2binary(str_state4), binary_state2)
    expect_equal(str2binary("ABCD", sep = ""), binary_state2)
    expect_equal(str2binary("A-B-C-D", sep = "-"), binary_state2)
    expect_equal(str2binary("ROOT", sep = "-", wt = "ROOT"), binary_state1)
    expect_equal(str2binary("A-B-C-D", sep = "-", n = 5), c( binary_state2, 0))
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

test_that("Transforming back and forth works", {
    expect_equal(int2str(str2int("A")), "A")
    expect_equal(int2str(str2int("WT")), "WT")
    expect_equal(binary2str(str2binary("A")), "A")
    expect_equal(binary2str(str2binary("AB"
        , n = 5, sep = ""), sep = ""), "AB")
    expect_equal(binary2str(str2binary("WT"
        , n = 5, sep = "")), "WT")

    expect_equal(binary2int(int2binary(16)), 16)
    expect_equal(binary2int(int2binary(0)), 0)
    expect_equal(str2int(int2str(13), n = 12), 13)
    expect_equal(str2int(int2str(0), n = 12), 0)

    expect_equal(int2binary(binary2int(c(0,0,0)), n = 3), c(0,0,0))
    expect_equal(int2binary(binary2int(c(0,1, 0)), n = 3), c(0,1, 0))
    expect_equal(str2binary(binary2str(c(0,0,0)), n = 7), c(0, 0, 0, 0, 0, 0, 0))
    expect_equal(str2binary(binary2str(c(0, 1, 0, 0)), n = 7), c(0,1, 0, 0, 0, 0, 0))
})


# test_that("Generate the sorted genoytpes correctly",{
#     expect_equal(generate_sorted_genotypes(0), c("WT"))
#     expect_equal(generate_sorted_genotypes(1), c("WT", "A"))
#     expect_equal(generate_sorted_genotypes(2), c("WT", "A", "B", "A, B"))
#     expect_equal(generate_sorted_genotypes(2, sep = ""), c("WT", "A", "B", "AB"))
#     expect_equal(generate_sorted_genotypes(3, sep = ""), c("WT", "A", "B", "C", "AB", "AC", "BC", "ABC"))
#     expect_error(generate_sorted_genotypes(-1, sep = ""), "Number of genes should be >= 0")
#     expect_equal(generate_sorted_genotypes(3, sep = "", index.return = TRUE), list(x = c("WT", "A", "B", "C", "AB", "AC", "BC", "ABC")
#         , ix = c(0, 1, 2, 4, 3, 5, 6, 7))
#     )

# })
