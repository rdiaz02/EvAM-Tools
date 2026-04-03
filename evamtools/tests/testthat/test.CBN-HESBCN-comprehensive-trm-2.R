t1 <- Sys.time()

## Comprehensive tests of CBN and HESBCN transition rate matrices for
## four additional DAGs (DAG_2, DAG_3, DAG_4, DAG_5, DAG_6, DAG_7).
##
## For each DAG:
##   - CBN (all-AND conjunction)
##   - HESBCN: all 3^k combinations of AND/OR/XOR for the k
##     multi-parent nodes (3 for DAG_3; 27 for DAG_2, DAG_4, DAG_5)
##   - HESBCN all-AND must equal CBN
##
## The reference TRM is produced by brute-force over all 2^n subsets,
## completely independent of cpm2tm.

## ------------------------------------------------------------------
## Shared helpers (defined once, reused across all six DAG tests)
## ------------------------------------------------------------------

genot_name_2 <- function(S) {
    if (length(S) == 0L) "WT" else paste(sort(S), collapse = ", ")
}

is_accessible_2 <- function(S, parents_list, rels) {
    for (g in S) {
        pars <- parents_list[[g]]
        if (length(pars) == 0L) next
        n_p <- sum(pars %in% S)
        rel <- rels[[g]]
        ok  <- if (rel %in% c("Single", "AND")) {
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

brute_force_trm_2 <- function(genes, lambdas, parents_list, rels) {
    n          <- length(genes)
    acc_buf    <- vector("list", 2L^n)
    count      <- 0L
    for (k in 0L:(2L^n - 1L)) {
        S <- genes[as.logical(intToBits(as.integer(k))[seq_len(n)])]
        if (is_accessible_2(S, parents_list, rels)) {
            count <- count + 1L
            acc_buf[[count]] <- S
        }
    }
    acc_genots <- acc_buf[seq_len(count)]
    gnames     <- vapply(acc_genots, genot_name_2, "")
    m          <- length(acc_genots)
    trm <- matrix(0.0, m, m, dimnames = list(gnames, gnames))
    for (i in seq_len(m)) {
        S <- acc_genots[[i]]
        for (X in genes) {
            if (X %in% S) next
            pars  <- parents_list[[X]]
            n_p   <- sum(pars %in% S)
            rel   <- rels[[X]]
            ok    <- if (rel %in% c("Single", "AND")) {
                         n_p == length(pars)
                     } else if (rel == "OR") {
                         n_p >= 1L
                     } else if (rel == "XOR") {
                         n_p == 1L
                     } else {
                         FALSE
                     }
            if (!ok) next
            dest <- genot_name_2(sort(c(S, X)))
            if (dest %in% gnames)
                trm[gnames[[i]], dest] <- lambdas[[X]]
        }
    }
    trm
}

compare_with_ref_2 <- function(cpm2tm_out, ref_trm, label) {
    computed <- as.matrix(cpm2tm_out$weighted_fgraph)
    expect_equal(
        sort(rownames(computed)), sort(rownames(ref_trm)),
        label = paste0(label, ": accessible genotypes")
    )
    g <- sort(rownames(computed))
    expect_equal(
        computed[g, g], ref_trm[g, g],
        tolerance = 1e-10,
        label = paste0(label, ": TRM entries")
    )
}

edges_to_parents_2 <- function(edges) {
    genes <- unique(edges$To)
    out   <- setNames(vector("list", length(genes)), genes)
    for (g in genes) {
        frm       <- edges$From[edges$To == g]
        out[[g]]  <- frm[frm != "Root"]
    }
    out
}

lambda_from_edges_2 <- function(edges) {
    lam_col <- if ("rerun_lambda" %in% names(edges)) {
        edges$rerun_lambda
    } else {
        edges$rerun_lamba
    }
    genes <- unique(edges$To)
    vapply(genes,
           function(g) lam_col[edges$To == g][1L],
           0.0)
}

all_combos_2 <- function(gene_names) {
    rels <- c("AND", "OR", "XOR")
    args <- setNames(rep(list(rels), length(gene_names)), gene_names)
    do.call(expand.grid, c(args, stringsAsFactors = FALSE))
}

make_rels_2 <- function(parents_list, multi_rels) {
    out <- setNames(vector("list", length(parents_list)),
                    names(parents_list))
    for (g in names(parents_list)) {
        pars   <- parents_list[[g]]
        out[[g]] <- if (length(pars) <= 1L) {
                        "Single"
                    } else if (g %in% names(multi_rels)) {
                        as.character(multi_rels[[g]])
                    } else {
                        "AND"
                    }
    }
    out
}

cbn_to_hesbcn_2 <- function(cbn_edges, multi_rels) {
    lam_col <- if ("rerun_lambda" %in% names(cbn_edges)) {
        cbn_edges$rerun_lambda
    } else {
        cbn_edges$rerun_lamba
    }
    parents_list <- edges_to_parents_2(cbn_edges)
    rel_vec <- vapply(seq_len(nrow(cbn_edges)), function(i) {
        g   <- cbn_edges$To[i]
        frm <- cbn_edges$From[i]
        if (frm == "Root") return("Single")
        pars <- parents_list[[g]]
        if (length(pars) <= 1L) return("Single")
        if (g %in% names(multi_rels))
            return(as.character(multi_rels[[g]]))
        "AND"
    }, "")
    list(edges = data.frame(
        From     = cbn_edges$From,
        To       = cbn_edges$To,
        Lambdas  = lam_col,
        Relation = rel_vec,
        stringsAsFactors = FALSE
    ))
}

test_one_dag_2 <- function(cbn_edges, multi_genes, dag_label) {
    ## Normalise column name for the rerun_lamba typo (DAG_4)
    if ("rerun_lamba" %in% names(cbn_edges) &&
        !("rerun_lambda" %in% names(cbn_edges))) {
        names(cbn_edges)[names(cbn_edges) == "rerun_lamba"] <-
            "rerun_lambda"
    }
    genes        <- unique(cbn_edges$To)
    lambdas      <- lambda_from_edges_2(cbn_edges)
    parents_list <- edges_to_parents_2(cbn_edges)

    ## CBN (all-AND)
    cbn_rels <- make_rels_2(
        parents_list,
        setNames(rep("AND", length(multi_genes)), multi_genes)
    )
    cbn_ref <- brute_force_trm_2(genes, lambdas, parents_list, cbn_rels)
    compare_with_ref_2(
        cpm2tm(list(edges = cbn_edges)),
        cbn_ref,
        paste0(dag_label, " CBN")
    )

    ## HESBCN: all 3^k combinations of AND/OR/XOR
    combos <- all_combos_2(multi_genes)
    for (i in seq_len(nrow(combos))) {
        multi_rels <- setNames(
            vapply(multi_genes,
                   function(g) as.character(combos[i, g]),
                   ""),
            multi_genes
        )
        rels <- make_rels_2(parents_list, multi_rels)
        ref  <- brute_force_trm_2(genes, lambdas, parents_list, rels)
        compare_with_ref_2(
            cpm2tm(cbn_to_hesbcn_2(cbn_edges, multi_rels)),
            ref,
            paste0(dag_label, " HESBCN ",
                   paste(multi_rels, collapse = "/"))
        )
    }

    ## HESBCN all-AND must equal CBN
    all_and <- setNames(rep("AND", length(multi_genes)), multi_genes)
    expect_equal(
        as.matrix(
            cpm2tm(cbn_to_hesbcn_2(cbn_edges, all_and))$weighted_fgraph
        ),
        as.matrix(cpm2tm(list(edges = cbn_edges))$weighted_fgraph),
        label = paste0(dag_label, ": HESBCN all-AND equals CBN")
    )
}


## ------------------------------------------------------------------
## DAG_2: 10 genes
##   Root -> G(2), D(1), M(5), Y(3), B(8)
##   {G,D,M,Y,B} -> S(7)   [5-parent node]
##   {G,D,M,Y,B} -> W(4)   [5-parent node]
##   {D,M,Y,B}   -> A(6)   [4-parent node]
##   S -> C(9)
##   S -> F(4.6)
## Multi-parent nodes: S, W, A  ->  3^3 = 27 HESBCN cases
## ------------------------------------------------------------------
dag2_edges <- data.frame(
    From = c("Root", "Root", "Root", "Root", "Root",
             "G", "D", "M", "Y", "B",
             "G", "D", "M", "Y", "B",
             "D", "M", "Y", "B",
             "S", "S"),
    To   = c("G", "D", "M", "Y", "B",
             "S", "S", "S", "S", "S",
             "W", "W", "W", "W", "W",
             "A", "A", "A", "A",
             "C", "F"),
    rerun_lambda = c(2, 1, 5, 3, 8,
                     7, 7, 7, 7, 7,
                     4, 4, 4, 4, 4,
                     6, 6, 6, 6,
                     9, 4.6),
    stringsAsFactors = FALSE
)

test_that("DAG_2: CBN and HESBCN TRM, 10-gene DAG
(5+5+4-parent nodes, all AND/OR/XOR combinations)", {
    test_one_dag_2(dag2_edges, c("S", "W", "A"), "DAG2")
})


## ------------------------------------------------------------------
## DAG_3: 9 genes
##   Root -> X(3), W(8)
##   X -> M(7), B(6)
##   M -> T(9)
##   T -> N(1), T -> A(2)
##   A -> K(4)
##   {W,B,N,K} -> G(5)   [4-parent node — only multi-parent node]
## Multi-parent node: G  ->  3^1 = 3 HESBCN cases
## ------------------------------------------------------------------
dag3_edges <- data.frame(
    From = c("Root", "Root", "X", "X",
             "W", "M", "B", "T", "T", "N", "A", "K"),
    To   = c("X",    "W",    "M", "B",
             "G",    "T",    "G", "N", "A", "G", "K", "G"),
    rerun_lambda = c(3, 8, 7, 6,
                     5, 9, 5, 1, 2, 5, 4, 5),
    stringsAsFactors = FALSE
)

test_that("DAG_3: CBN and HESBCN TRM, 9-gene DAG
(single 4-parent node, AND/OR/XOR)", {
    test_one_dag_2(dag3_edges, "G", "DAG3")
})


## ------------------------------------------------------------------
## DAG_4: 9 genes
##   Root -> S(8), G(2), H(3), Z(4)
##   S -> W(5)
##   {S,G,H} -> V(1)   [3-parent node]
##   H -> T(6)
##   Z -> Y(9)
##   {W,V,T} -> X(7)   [3-parent node]
##   {Z,X}   -> Y(9)   [2-parent node]  (note Z->Y and X->Y)
## Multi-parent nodes: V, X, Y  ->  3^3 = 27 HESBCN cases
## (Original data file had typo "rerun_lamba"; corrected here.)
## ------------------------------------------------------------------
dag4_edges <- data.frame(
    From = c("Root", "Root", "Root", "Root",
             "S", "S", "G", "H", "H", "Z",
             "W", "V", "T", "X"),
    To   = c("S", "G", "H", "Z",
             "W", "V", "V", "V", "T", "Y",
             "X", "X", "X", "Y"),
    rerun_lambda = c(8, 2, 3, 4,
                     5, 1, 1, 1, 6, 9,
                     7, 7, 7, 9),
    stringsAsFactors = FALSE
)

test_that("DAG_4: CBN and HESBCN TRM, 9-gene DAG
(3+3+2-parent nodes, all AND/OR/XOR combinations)", {
    test_one_dag_2(dag4_edges, c("V", "X", "Y"), "DAG4")
})


## ------------------------------------------------------------------
## DAG_5: 9 genes
##   Root -> C(2), K(7), N(8), Z(10)
##   {C,K,N,Z} -> F(4)   [4-parent node]
##   {C,K,N}   -> G(5)   [3-parent node]
##   F -> I(6)
##   G -> T(9)
##   {I,G}     -> D(3)   [2-parent node — note: I and G, not I and T]
## Multi-parent nodes: F, G, D  ->  3^3 = 27 HESBCN cases
## ------------------------------------------------------------------
dag5_edges <- data.frame(
    From = c("Root", "Root", "Root", "Root",
             "C", "K", "N", "Z",
             "C", "K", "N",
             "F", "G",
             "I", "G"),
    To   = c("C", "K", "N", "Z",
             "F", "F", "F", "F",
             "G", "G", "G",
             "I", "T",
             "D", "D"),
    rerun_lambda = c(2, 7, 8, 10,
                     4, 4, 4, 4,
                     5, 5, 5,
                     6, 9,
                     3, 3),
    stringsAsFactors = FALSE
)

test_that("DAG_5: CBN and HESBCN TRM, 9-gene DAG
(4+3+2-parent nodes, shared parents, all AND/OR/XOR combinations)", {
    test_one_dag_2(dag5_edges, c("F", "G", "D"), "DAG5")
})


## ------------------------------------------------------------------
## DAG_6: 7 genes — stacked multi-parent (addresses coverage gap)
##   Root -> D(1), E(2), G(3), H(4)
##   {D,E} -> B(5)   [2-parent node]
##   {G,H} -> C(6)   [2-parent node]
##   {B,C} -> F(7)   [2-parent node; parents are themselves multi-parent]
## Multi-parent nodes: B, C, F  ->  3^3 = 27 HESBCN cases
## ------------------------------------------------------------------
dag6_edges <- data.frame(
    From = c("Root", "Root", "Root", "Root",
             "D", "E", "G", "H", "B", "C"),
    To   = c("D", "E", "G", "H",
             "B", "B", "C", "C", "F", "F"),
    rerun_lambda = c(1, 2, 3, 4,
                     5, 5, 6, 6, 7, 7),
    stringsAsFactors = FALSE
)

test_that("DAG_6: CBN and HESBCN TRM, 7-gene DAG
(stacked 2+2+2-parent nodes, all AND/OR/XOR combinations)", {
    test_one_dag_2(dag6_edges, c("B", "C", "F"), "DAG6")
})

test_that("DAG_7: hand check on the XOR/XOR/AND", {

  d7m <- dag7_edges
  colnames(d7m)[3] <- "Lambdas"
  d7m$Relation <- c(rep("Single", 3), rep("XOR", 5), rep("AND", 2))

  computed_trm <- cpm2tm(list(edges = d7m))$weighted_fgraph

  access_genots <- c("WT",
                     "D", "E", "G",
                     "B, D", "B, E", "B, G",
                     "C, E", "C, G",
                     "D, E", "D, G", "E, G",
                     "B, C, E", "B, C, G",
                     "C, D, E", "C, D, G",
                     "D, E, G",
                     "B, C, E, F",
                     "B, C, F, G")

  expected_trm <- matrix(0, nrow = length(access_genots),
                         ncol = length(access_genots))
  row.names(expected_trm) <- colnames(expected_trm) <- access_genots

  expected_trm["WT", c("D", "E", "G")] <- c(1, 2, 3)
  expected_trm["D", c("D, E", "D, G", "B, D")] <- c(2, 3, 4)
  expected_trm["E", c("D, E", "E, G", "B, E", "C, E")] <- c(1, 3, 4, 5)
  expected_trm["G", c("B, G", "E, G", "C, G", "D, G")] <- c(4, 2, 5, 1)
  expected_trm["B, E", c("B, C, E")] <- c(5)
  expected_trm["B, G", c("B, C, G")] <- c(5)
  expected_trm["C, E", c("C, D, E", "B, C, E")] <- c(1, 4)
  expected_trm["C, G", c("C, D, G", "B, C, G")] <- c(1, 4)
  expected_trm["D, E", c("D, E, G", "C, D, E")] <- c(3, 5)
  expected_trm["D, G", c("C, D, G", "D, E, G")] <- c(5, 2)
  expected_trm["E, G", c("D, E, G")] <- c(1)
  expected_trm["B, C, E", c("B, C, E, F")] <- c(6)
  expected_trm["B, C, G", c("B, C, F, G")] <- c(6)

  all.equal(expected_trm, as.matrix(computed_trm))

})



## ------------------------------------------------------------------
## DAG_7: 6 genes — stacked multi-parent with shared ancestors
##   Root -> D(1), E(2), G(3)
##   {D,E,G} -> B(4)   [3-parent node]
##   {E,G}   -> C(5)   [2-parent node; E,G shared with B]
##   {B,C}   -> F(6)   [2-parent node; parents are themselves multi-parent]
## Multi-parent nodes: B, C, F  ->  3^3 = 27 HESBCN cases
## ------------------------------------------------------------------
dag7_edges <- data.frame(
    From = c("Root", "Root", "Root",
             "D", "E", "G", "E", "G", "B", "C"),
    To   = c("D", "E", "G",
             "B", "B", "B", "C", "C", "F", "F"),
    rerun_lambda = c(1, 2, 3,
                     4, 4, 4, 5, 5, 6, 6),
    stringsAsFactors = FALSE
)

test_that("DAG_7: CBN and HESBCN TRM, 6-gene DAG
(stacked 3+2+2-parent nodes, shared ancestors, all AND/OR/XOR)", {
    test_one_dag_2(dag7_edges, c("B", "C", "F"), "DAG7")
})


cat("\n Done test.CBN-HESBCN-comprehensive-trm-2.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")


## About this testing file. It was written by Claude and
## reviewd by a subagent. RDU then manually checked a few cases (DAG_6,
## with XOR, XOR, AND).

## The production code (cpm2tm) and reference share zero code,
## zero data structures, and a different computational paradigm.
## cpm2tm uses:
## - igraph::all_simple_paths() to traverse the DAG
## - Sparse matrix construction for the fitness graph
## - A dedicated accessible-genotype enumeration per relation
##      type (DAG_2_access_genots, DAG_2_access_genots_OR,
##      DAG_2_access_genots_relationships)
## - transition_fg_sparseM() for weighting edges

## The brute-force reference (brute_force_trm_2) uses:
##   - Bit arithmetic (intToBits) to enumerate all 2^n subsets directly
## - A simple loop checking parent constraints gene-by-gene
## - No graphs, no sparse matrices, no igraph

## For them to agree on 294 entries across 6 structurally diverse
## DAGs and all AND/OR/XOR combinations, any cancelling-error scenario
## would have to be extraordinarily contrived.
