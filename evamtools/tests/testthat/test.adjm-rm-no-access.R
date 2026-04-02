t1 <- Sys.time()

## Direct unit tests for adjm_rm_no_access().
##
## The function takes a weighted fitness graph (TRM as a named matrix,
## rows = from, cols = to) and iteratively removes any genotype whose
## column-sum is zero — meaning nothing transitions into it — except
## WT, which is always the starting state and has no incoming edges
## by definition.
##
## The while loop is needed because removing one genotype can orphan
## others that were only reachable through it (cascading removal).

## As of 2026-04-02 adj_rm_no_access is only used for OncoBN
## models that seem to be the only ones with the need for this function
## (none of the other models show the strange behavior that requires
## the use of this function). See function cpm2tm.



test_that("adjm_rm_no_access: no removal needed", {
    ## All non-WT genotypes have at least one incoming edge.
    ## Matrix should be returned unchanged.
    m <- matrix(
        c(0,.1,.8, 0,
          0, 0, 0,.7,
          0, 0, 0,.4,
          0, 0, 0, 0),
        nrow = 4, ncol = 4, byrow = TRUE,
        dimnames = list(c("WT", "A", "B", "A, B"), c("WT", "A", "B", "A, B"))
    )
    result <- evamtools:::adjm_rm_no_access(m)
    expect_equal(sort(rownames(result)), c("A", "A, B", "B", "WT"))
    expect_equal(result, m)
})


test_that("adjm_rm_no_access: WT is never removed", {
    ## WT always has colSum = 0 (nothing transitions into it).
    ## It must be preserved regardless.
    m <- matrix(
        c(0, .91,
          0, 0),
        nrow = 2, ncol = 2, byrow = TRUE,
        dimnames = list(c("WT", "A"), c("WT", "A"))
    )
    result <- evamtools:::adjm_rm_no_access(m)
    expect_true("WT" %in% rownames(result))
    expect_equal(sort(rownames(result)), c("A", "WT"))
})


test_that("adjm_rm_no_access: single direct orphan removed", {
    ## B has no incoming edges (colSum = 0); A does.
    ## Only B should be removed.
    m <- matrix(
        c(0, .8, 0,
          0, 0, 0,
          0, 0, 0),
        nrow = 3, ncol = 3, byrow = TRUE,
        dimnames = list(c("WT", "A", "B"), c("WT", "A", "B"))
    )
    result <- evamtools:::adjm_rm_no_access(m)
    expect_equal(sort(rownames(result)), c("A", "WT"))
    expect_equal(result["WT", "A"], 0.8)
})


test_that("adjm_rm_no_access: cascading removal exercises the while loop", {
    ## I am not sure this can really happen with OncoBN.
    ## B has no incoming edges.
    ## A, B is only reachable via B, so after B is removed A, B becomes
    ## orphaned and must be removed in the next loop iteration.
    ## Only WT and A should remain.
    ##
    ##   WT -> A (rate 1)
    ##   B  -> A, B (rate 2)   [B itself is orphaned]
    m <- matrix(
        c(0, .5, 0,  0,
          0, 0,  0,  0,
          0, 0,  0,  .8,
          0, 0,  0,  0),
        nrow = 4, ncol = 4, byrow = TRUE,
        dimnames = list(c("WT", "A", "B", "A, B"),
                        c("WT", "A", "B", "A, B"))
    )
    result <- evamtools:::adjm_rm_no_access(m)
    expect_equal(sort(rownames(result)), c("A", "WT"))
    expect_equal(result["WT", "A"], 0.5)
})


test_that("adjm_rm_no_access: multiple independent orphans removed in one pass", {
    ## Both B and C have no incoming edges and are independent.
    ## Both should be removed; WT and A remain.
    m <- matrix(
        c(0, .4, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0),
        nrow = 4, ncol = 4, byrow = TRUE,
        dimnames = list(c("WT", "A", "B", "C"),
                        c("WT", "A", "B", "C"))
    )
    result <- evamtools:::adjm_rm_no_access(m)
    expect_equal(sort(rownames(result)), c("A", "WT"))
})


test_that("adjm_rm_no_access: longer cascade (3 levels)", {
    ## Actually let's set up:
    ##   WT -> A (p = 0.7)
    ##   D  -> BD (p = 0.6)   [D is orphaned]
    ##   BD  -> BCD (p = 0.8)   [after D removed, BD orphaned, then BCD orphaned]
    ##
    ## Expected survivors: WT and A.
    m <- matrix(
        c(0, .7, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 0, .6, 0,
          0, 0, 0, 0, 0.8,
          0, 0, 0, 0, 0),
        nrow = 5, ncol = 5, byrow = TRUE,
        dimnames = list(c("WT", "A", "D", "B, D ", "B, C, D"),
                        c("WT", "A", "D", "B, D ", "B, C, D"))
    )
    result <- evamtools:::adjm_rm_no_access(m)
    expect_equal(sort(rownames(result)), c("A", "WT"))
    expect_equal(result["WT", "A"], 0.7)
})


test_that("adjm_rm_no_access: all non-WT genotypes removed", {
    ## No transitions at all: only WT should survive.
    m <- matrix(
        c(0, 0, 0,
          0, 0, 0,
          0, 0, 0),
        nrow = 3, ncol = 3, byrow = TRUE,
        dimnames = list(c("WT", "A", "B"), c("WT", "A", "B"))
    )
    expect_message(evamtools:::adjm_rm_no_access(m),
                   "No transitions from WT to anything in this model")

})


cat("\n Done test.adjm-rm-no-access.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
