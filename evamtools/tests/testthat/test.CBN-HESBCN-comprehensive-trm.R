t1 <- Sys.time()

## TODO
## - Fix the multiple statements per line
## - Fix the no spaces around "="
## - brute_force_trm: Use sparse matrix
## - brute_force_trm: Return in the "canonicalized order"

## - Myself, by hand:
##     - check list of accessible genotypes
##     - check the entries


## Comprehensive test of CBN and HESBCN transition rate matrices using a
## single complex 9-gene DAG with parent-set sizes 4, 3, 2, and 1.
##
## DAG structure (lambdas assigned alphabetically, 2 to 10):
##
##   Root -> C(2), K(7), N(8), Z(10)       [Single, Root children]
##   C, K, N, Z -> F(4)                    [4-parent node]
##   C, K, N    -> G(5)                    [3-parent node]
##   F -> I(6)                             [Single]
##   G -> T(9)                             [Single]
##   I, T -> D(3)                          [2-parent node]
##
## C, K, N are parents of BOTH F and G — deliberately chosen so
## XOR constraints interact across the two branches.
##
## CBN: one case (all AND).
## HESBCN: 9 cases — 3 pure (all AND / all OR / all XOR) plus all 6
## permutations of AND, OR, XOR assigned to (F, G, D).
##
## The reference TRM is produced by brute-force over all 2^9 subsets
## using a simple, direct implementation of the accessibility rules and
## rate assignment.  It is completely independent of cpm2tm.

test_that("CBN and HESBCN: comprehensive TRM, 9-gene DAG
(4+3+2-parent nodes, all AND/OR/XOR combinations)", {

    genes   <- c("C", "D", "F", "G", "I", "K", "N", "T", "Z")
    lambdas <- c(C=2, D=3, F=4, G=5, I=6, K=7, N=8, T=9, Z=10)

    ## Parent sets of each gene (empty = Root child, always acquirable)
    parents_list <- list(
        C = character(0), K = character(0), N = character(0), Z = character(0),
        F = c("C", "K", "N", "Z"),
        G = c("C", "K", "N"),
        I = "F",
        T = "G",
        D = c("I", "T")
    )

    ## ------------------------------------------------------------------
    ## Reference: brute-force over all 2^9 = 512 subsets
    ## ------------------------------------------------------------------

    genot_name <- function(S) {
        if (length(S) == 0L) "WT" else paste(sort(S), collapse = ", ")
    }

    ## rels: named list, gene -> "Single" | "AND" | "OR" | "XOR"
    ## A genotype S is accessible iff for every gene g in S, g's parent
    ## constraint (given rels[[g]]) is satisfied by the rest of S.
    is_accessible <- function(S, rels) {
        for (g in S) {
            pars <- parents_list[[g]]
            if (length(pars) == 0L) next   ## Root child: always OK
            n_p <- sum(pars %in% S)
            rel <- rels[[g]]
            ok <- if (rel %in% c("Single", "AND")) {
                      n_p == length(pars)
                  } else if (rel == "OR") {
                      n_p >= 1L
                  } else if (rel == "XOR") {
                      n_p == 1L
                  } else {
                      stop("Unknown relation: ", rel)
                  }
            if (!ok) return(FALSE)
        }
        TRUE
    }

    ## Enumerate accessible genotypes, then build TRM.
    ## TRM[from, to] = lambda of added gene X, when:
    ##   (a) X's own parent constraint is satisfied in "from", AND
    ##   (b) the resulting genotype "to" is itself accessible
    ##       (needed for XOR: adding a second XOR-parent blocks a child)
    brute_force_trm <- function(rels) {
        n <- length(genes)
        acc_genots <- list()
        for (k in 0L:(2L^n - 1L)) {
            S <- genes[as.logical(intToBits(as.integer(k))[seq_len(n)])]
            if (is_accessible(S, rels))
                acc_genots <- c(acc_genots, list(S))
        }
        gnames <- vapply(acc_genots, genot_name, "")
        m <- length(acc_genots)

        trm <- matrix(0.0, m, m, dimnames = list(gnames, gnames))
        for (i in seq_len(m)) {
            S <- acc_genots[[i]]
            for (X in genes) {
                if (X %in% S) next
                pars <- parents_list[[X]]
                n_p  <- sum(pars %in% S)
                rel  <- rels[[X]]
                can_add <- if (rel %in% c("Single", "AND")) {
                               n_p == length(pars)
                           } else if (rel == "OR") {
                               n_p >= 1L
                           } else if (rel == "XOR") {
                               n_p == 1L
                           } else {
                               FALSE
                           }
                if (!can_add) next
                dest <- genot_name(sort(c(S, X)))
                ## XOR on another gene may make dest inaccessible
                if (dest %in% gnames)
                    trm[gnames[[i]], dest] <- lambdas[[X]]
            }
        }
        trm
    }

    ## ------------------------------------------------------------------
    ## Comparison: accessible genotypes + full TRM matrix
    ## ------------------------------------------------------------------
    compare_with_ref <- function(cpm2tm_out, ref_trm, label) {
        computed <- as.matrix(cpm2tm_out$weighted_fgraph)
        ## (1) Accessible genotype sets must match exactly
        expect_equal(
            sort(rownames(computed)), sort(rownames(ref_trm)),
            label = paste0(label, ": accessible genotypes")
        )
        ## (2) Full TRM entry-by-entry, aligned by sorted names
        g <- sort(rownames(computed))
        expect_equal(
            computed[g, g], ref_trm[g, g],
            tolerance = 1e-10,
            label = paste0(label, ": TRM entries")
        )
    }

    ## ------------------------------------------------------------------
    ## Shared edge skeleton (From/To/lambda values are the same for
    ## CBN and all HESBCN variants; only the column name and Relation differ)
    ## ------------------------------------------------------------------
    edge_from <- c(
        "Root", "Root", "Root", "Root",  ## -> C, K, N, Z
        "C", "K", "N", "Z",              ## -> F (4 parents)
        "C", "K", "N",                   ## -> G (3 parents)
        "F", "G",                        ## -> I, T (single)
        "I", "T"                         ## -> D (2 parents)
    )
    edge_to <- c(
        "C", "K", "N", "Z",
        "F", "F", "F", "F",
        "G", "G", "G",
        "I", "T",
        "D", "D"
    )
    edge_lam <- c(
        2, 7, 8, 10,         ## C, K, N, Z
        4, 4, 4, 4,          ## F
        5, 5, 5,             ## G
        6, 9,                ## I, T
        3, 3                 ## D
    )

    ## ------------------------------------------------------------------
    ## CBN (rerun_lambda format; semantics identical to HESBCN all-AND)
    ## ------------------------------------------------------------------
    cbn_input <- list(edges = data.frame(
        From = edge_from, To = edge_to,
        rerun_lambda = edge_lam,
        stringsAsFactors = FALSE
    ))
    cbn_rels <- list(
        C="Single", D="AND", F="AND", G="AND",
        I="Single", K="Single", N="Single", T="Single", Z="Single"
    )
    cbn_ref <- brute_force_trm(cbn_rels)
    compare_with_ref(cpm2tm(cbn_input), cbn_ref, "CBN")

    ## ------------------------------------------------------------------
    ## HESBCN helpers
    ## ------------------------------------------------------------------
    make_hesbcn <- function(rel_F, rel_G, rel_D) {
        rel_vec <- c(
            rep("Single", 4),            ## Root -> C, K, N, Z
            rep(rel_F, 4),               ## -> F
            rep(rel_G, 3),               ## -> G
            "Single", "Single",          ## -> I, T
            rep(rel_D, 2)                ## -> D
        )
        edges <- data.frame(
            From = edge_from, To = edge_to,
            Lambdas = edge_lam, Relation = rel_vec,
            stringsAsFactors = FALSE
        )
        parent_set <- c(
            C="Single", D=rel_D, F=rel_F, G=rel_G,
            I="Single", K="Single", N="Single", T="Single", Z="Single"
        )
        list(edges = edges, parent_set = parent_set)
    }

    make_rels <- function(rel_F, rel_G, rel_D) {
        list(
            C="Single", D=rel_D, F=rel_F, G=rel_G,
            I="Single", K="Single", N="Single", T="Single", Z="Single"
        )
    }

    ## ------------------------------------------------------------------
    ## HESBCN cases: 3 pure + all 6 permutations of AND/OR/XOR
    ## Column order: (rel for F's 4 parents, rel for G's 3 parents,
    ##                rel for D's 2 parents)
    ## ------------------------------------------------------------------
    hesbcn_cases <- list(
        list("AND", "AND", "AND", "HESBCN all AND"),
        list("OR",  "OR",  "OR",  "HESBCN all OR"),
        list("XOR", "XOR", "XOR", "HESBCN all XOR"),
        list("AND", "OR",  "XOR", "HESBCN AND/OR/XOR"),
        list("AND", "XOR", "OR",  "HESBCN AND/XOR/OR"),
        list("OR",  "AND", "XOR", "HESBCN OR/AND/XOR"),
        list("OR",  "XOR", "AND", "HESBCN OR/XOR/AND"),
        list("XOR", "AND", "OR",  "HESBCN XOR/AND/OR"),
        list("XOR", "OR",  "AND", "HESBCN XOR/OR/AND")
    )

    for (cc in hesbcn_cases) {
        rF <- cc[[1]]; rG <- cc[[2]]; rD <- cc[[3]]; lbl <- cc[[4]]
        ref <- brute_force_trm(make_rels(rF, rG, rD))
        compare_with_ref(cpm2tm(make_hesbcn(rF, rG, rD)), ref, lbl)
    }

    ## HESBCN all-AND must give the same weighted_fgraph as CBN
    hesbcn_and <- make_hesbcn("AND", "AND", "AND")
    expect_equal(
        as.matrix(cpm2tm(hesbcn_and)$weighted_fgraph),
        as.matrix(cpm2tm(cbn_input)$weighted_fgraph),
        label = "HESBCN all-AND equals CBN"
    )
})

cat("\n Done test.CBN-HESBCN-comprehensive-trm.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
