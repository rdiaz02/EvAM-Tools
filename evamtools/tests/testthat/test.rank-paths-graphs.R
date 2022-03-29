t1 <- Sys.time()

# #' Test the procedure of node labeling for the plot_genot_fg function
# #' It labels nodes as those involved in the most transited paths

# #' that star in the root of the tree and that lead to any leave 
test_that("Modifying matrix with top_paths",{
    base_adj_matrix <- matrix(0, nrow = 6, ncol = 6)
    rownames(base_adj_matrix) <- colnames(base_adj_matrix) <- c("WT", "A", "A, B",
                                                      "A, C", "C", "C, D")
    adj_matrix <- base_adj_matrix

    adj_matrix["WT", "A"] <- 1000
    adj_matrix["WT", "C"] <- 500
    adj_matrix["A", "A, B"] <- 500
    adj_matrix["A", "A, C"] <- 200
    adj_matrix["C", "A, C"] <- 100
    adj_matrix["C", "C, D"] <- 200

    g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "directed",
                                             weighted = TRUE)
    
    all_genotypes <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    all_paths <- rank_paths(g, log_weights = FALSE)$paths

    l1 <- c("WT", "A", "A, B", "", "", "")
    names(l1) <- all_genotypes
    adj1 <- base_adj_matrix
    adj1["WT", "A"] <- 1000
    # adj1["WT", "C"] <- 1
    adj1["A", "A, B"] <- 500
    # adj1["A", "A, C"] <- 1
    # adj1["C", "A, C"] <- 1
    # adj1["C", "C, D"] <- 1
    adj1 <- Matrix(adj1, sparse = TRUE)

    l2 <- c("WT", "A", "A, B", "A, C", "", "")
    names(l2) <- all_genotypes
    adj2 <- base_adj_matrix
    adj2["WT", "A"] <- 1000
    # adj2["WT", "C"] <- 1
    adj2["A", "A, B"] <- 500
    adj2["A", "A, C"] <- 200
    # adj2["C", "A, C"] <- 1
    # adj2["C", "C, D"] <- 1
    adj2 <- Matrix(adj2, sparse = TRUE)

    l3 <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    names(l3) <- all_genotypes
    adj3 <- base_adj_matrix
    adj3["WT", "A"] <- 1000
    adj3["WT", "C"] <- 500
    adj3["A", "A, B"] <- 500
    adj3["A", "A, C"] <- 200
    # adj3["C", "A, C"] <- 1
    adj3["C", "C, D"] <- 200
    adj3 <- Matrix(adj3, sparse = TRUE)

    l4 <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    names(l4) <- all_genotypes
    adj4 <- base_adj_matrix
    adj4["WT", "A"] <- 1000
    adj4["WT", "C"] <- 500
    adj4["A", "A, B"] <- 500
    adj4["A", "A, C"] <- 200
    adj4["C", "A, C"] <- 100
    adj4["C", "C, D"] <- 200
    adj4 <- Matrix(adj4, sparse = TRUE)

    new_matrix_0 <-    compute_matrix_from_top_paths(g, all_paths, top_paths = 0)
    new_matrix_1 <-    compute_matrix_from_top_paths(g, all_paths, top_paths = 1)
    new_matrix_2 <-    compute_matrix_from_top_paths(g, all_paths, top_paths = 2)
    new_matrix_3 <-    compute_matrix_from_top_paths(g, all_paths, top_paths = 3)
    new_matrix_4 <-    compute_matrix_from_top_paths(g, all_paths, top_paths = 4)
    new_matrix_40 <-   compute_matrix_from_top_paths(g, all_paths, top_paths = 40)
    new_matrix_null <- compute_matrix_from_top_paths(g, all_paths, top_paths = NULL)

    expect_equal(as.matrix(new_matrix_1), as.matrix(adj1), ignore_attr = TRUE)
    expect_equal(as.matrix(new_matrix_2), as.matrix(adj2), ignore_attr = TRUE)
    expect_equal(as.matrix(new_matrix_3), as.matrix(adj3), ignore_attr = TRUE)
    expect_equal(as.matrix(new_matrix_4), as.matrix(adj4), ignore_attr = TRUE)
    expect_equal(as.matrix(new_matrix_40), as.matrix(adj4), ignore_attr = TRUE)
    expect_equal(as.matrix(new_matrix_null), as.matrix(adj4), ignore_attr = TRUE)
    expect_equal(as.matrix(new_matrix_0), as.matrix(adj4), ignore_attr = TRUE)
})

test_that("Returns the correct number and type of vertex labels", {

    base_adj_matrix <- matrix(0, nrow = 6, ncol = 6)
    rownames(base_adj_matrix) <- colnames(base_adj_matrix) <- c("WT", "A", "A, B",
                                                      "A, C", "C", "C, D")
    adj_matrix <- base_adj_matrix

    adj_matrix["WT", "A"] <- 1000
    adj_matrix["WT", "C"] <- 500
    adj_matrix["A", "A, B"] <- 500
    adj_matrix["A", "A, C"] <- 200
    adj_matrix["C", "A, C"] <- 100
    adj_matrix["C", "C, D"] <- 200

    g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
    
    all_genotypes <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    all_paths <- rank_paths(g, log_weights = FALSE)$paths

    l1 <- c("WT", "A", "A, B", "", "", "")
    names(l1) <- all_genotypes
    l2 <- c("WT", "A", "A, B", "A, C", "", "")
    names(l2) <- all_genotypes
    l3 <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    names(l3) <- all_genotypes
    l4 <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    names(l4) <- all_genotypes

    adj_0 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 0)
    g0 <- igraph::graph_from_adjacency_matrix(adj_0
            , weighted = TRUE)

    adj_1 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 1)
    g1 <- igraph::graph_from_adjacency_matrix(adj_1
            , weighted = TRUE)
    
    adj_2 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 2)
    g2 <- igraph::graph_from_adjacency_matrix(adj_2
            , weighted = TRUE)
    
    adj_3 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 3)
    g3 <- igraph::graph_from_adjacency_matrix(adj_3
            , weighted = TRUE)

    adj_4 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 4)
    g4 <- igraph::graph_from_adjacency_matrix(adj_4
            , weighted = TRUE)

    adj_40 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 40)
    g40 <- igraph::graph_from_adjacency_matrix(adj_40
            , weighted = TRUE)
    
    adj_null <- compute_matrix_from_top_paths(g, all_paths, top_paths = NULL)
    g_null <- igraph::graph_from_adjacency_matrix(adj_null
            , weighted = TRUE)

    labels_0 <- compute_vertex_labels(g0, all_paths, top_paths = 0)
    labels_1 <- compute_vertex_labels(g1, all_paths, top_paths = 1)
    labels_2 <- compute_vertex_labels(g2, all_paths, top_paths = 2)
    labels_3 <- compute_vertex_labels(g3, all_paths, top_paths = 3)
    labels_4 <- compute_vertex_labels(g4, all_paths, top_paths = 4)
    labels_40 <- compute_vertex_labels(g40, all_paths, top_paths = 40)
    labels_null <- compute_vertex_labels(g_null, all_paths, top_paths = NULL)

    expect_equal(labels_1$vertex_labels, l1)
    expect_equal(labels_2$vertex_labels, l2)
    expect_equal(labels_3$vertex_labels, l3)
    expect_equal(labels_4$vertex_labels, l4)
    expect_equal(labels_4$vertex_labels, labels_40$vertex_labels)
    expect_equal(labels_4$vertex_labels, labels_null$vertex_labels)
    expect_equal(labels_4$vertex_labels, labels_0$vertex_labels)

    expect_equal(all(labels_1$edge_labels == ""), TRUE)
    expect_equal(all(labels_2$edge_labels == ""), TRUE)
    expect_equal(all(labels_3$edge_labels == ""), TRUE)
    expect_equal(all(labels_4$edge_labels == ""), TRUE)
    expect_equal(all(labels_40$edge_labels == ""), TRUE)
    expect_equal(all(labels_null$edge_labels == ""), TRUE)
})

test_that("Returns the correct number and type of edges labels", {

    base_adj_matrix <- matrix(0, nrow = 6, ncol = 6)
    rownames(base_adj_matrix) <- colnames(base_adj_matrix) <- c("WT", "A", "A, B",
                                                      "A, C", "C", "C, D")
    adj_matrix <- base_adj_matrix

    adj_matrix["WT", "A"] <- 1000
    adj_matrix["WT", "C"] <- 500
    adj_matrix["A", "A, B"] <- 500
    adj_matrix["A", "A, C"] <- 200
    adj_matrix["C", "A, C"] <- 100
    adj_matrix["C", "C, D"] <- 200

    g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
    
    all_genotypes <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    all_paths <- rank_paths(g, log_weights = FALSE)$paths

    l1 <- c("+A", "+B")
    l2 <- c("+A", "+B", "+C")
    l3 <- c("+A", "+B", "+C", "+C", "+D")
    l4 <- c("+A", "+B", "+C", "+A", "+C", "+D")

    adj_0 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 0)
    g0 <- igraph::graph_from_adjacency_matrix(adj_0
            , weighted = TRUE)

    adj_1 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 1)
    g1 <- igraph::graph_from_adjacency_matrix(adj_1
            , weighted = TRUE)
    
    adj_2 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 2)
    g2 <- igraph::graph_from_adjacency_matrix(adj_2
            , weighted = TRUE)
    
    adj_3 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 3)
    g3 <- igraph::graph_from_adjacency_matrix(adj_3
            , weighted = TRUE)

    adj_4 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 4)
    g4 <- igraph::graph_from_adjacency_matrix(adj_4
            , weighted = TRUE)

    adj_40 <- compute_matrix_from_top_paths(g, all_paths, top_paths = 40)
    g40 <- igraph::graph_from_adjacency_matrix(adj_40
            , weighted = TRUE)
    
    adj_null <- compute_matrix_from_top_paths(g, all_paths, top_paths = NULL)
    g_null <- igraph::graph_from_adjacency_matrix(adj_null
            , weighted = TRUE)

    labels_0 <- compute_vertex_labels(g0, all_paths, top_paths = 0
        , type="acquisition")
    labels_1 <- compute_vertex_labels(g1, all_paths, top_paths = 1
        , type="acquisition")
    labels_2 <- compute_vertex_labels(g2, all_paths, top_paths = 2
        , type="acquisition")
    labels_3 <- compute_vertex_labels(g3, all_paths, top_paths = 3
        , type="acquisition")
    labels_4 <- compute_vertex_labels(g4, all_paths, top_paths = 4
        , type="acquisition")
    labels_40 <- compute_vertex_labels(g40, all_paths, top_paths = 40
        , type="acquisition")
    labels_null <- compute_vertex_labels(g_null, all_paths, top_paths = NULL
        , type="acquisition")

    expect_equal(labels_1$edge_labels, l1)
    expect_equal(labels_2$edge_labels, l2)
    expect_equal(labels_3$edge_labels, l3)
    expect_equal(labels_4$edge_labels, l4)
    expect_equal(labels_4$edge_labels, labels_40$edge_labels)
    expect_equal(labels_4$edge_labels, labels_null$edge_labels)
    expect_equal(labels_4$edge_labels, labels_0$edge_labels)

    base_vertex <- rep("", length(all_genotypes))
    expect_equal(labels_1$vertex_labels, base_vertex)
    expect_equal(labels_2$vertex_labels, base_vertex)
    expect_equal(labels_3$vertex_labels, base_vertex)
    expect_equal(labels_4$vertex_labels, base_vertex)
    expect_equal(labels_4$vertex_labels, base_vertex)
    expect_equal(labels_4$vertex_labels, base_vertex)
    expect_equal(labels_4$vertex_labels, base_vertex)
})


test_that("Simple test that we recover correct rank of paths and their prob", {
    m1 <- matrix(0, nrow = 7, ncol = 7)
    rownames(m1) <- colnames(m1) <- c("WT", "A", "B", "A, B",
                                      "B, D", "B, E", "A, B, C")
    m1["WT", "A"] <- .1
    m1["A", "A, B"] <- 1
    m1["A, B", "A, B, C"] <- 1
    m1["WT", "B"] <- 0.9
    m1["B", "B, D"] <- 0.6
    m1["B", "B, E"] <- 0.4
    m1 <- Matrix(m1, sparse = TRUE)
    ## m1
    ## So three paths, with these probabilities,
    ## WT -> A -> AB -> ABC : 0.1
    ## WT -> B -> BD        : 0.9 * 0.6
    ## WT -> B -> BE        : 0.9 * 0.4

    
    m1g <- igraph::graph_from_adjacency_matrix(m1, weighted = TRUE,
                                               mode = "directed")
    rp <- rank_paths(m1g, log_weights = TRUE)

    expect_identical(lapply(rp$paths, igraph::as_ids),
                     list(c("WT", "B", "B, D"),
                          c("WT", "B", "B, E"),
                          c("WT", "A", "A, B", "A, B, C")))
    expect_equal(exp(rp$weights), c(0.9 * 0.6, 0.9 * 0.4, 0.1))

    expect_equivalent(paths_probs_2_df(trans_mat_2_paths_probs(m1),
                                       order = "path")[, "Prob"],
                 c(0.1, 0.54, 0.36))
                 

    m2 <- matrix(0, nrow = 9, ncol = 9)
    rownames(m2) <- colnames(m2) <- c("WT", "A", "B",
                                      "B, C",
                                      "B, D",
                                      "B, D, E", "B, D, F",
                                      ## zeros, but ensure working OK
                                      "B, C, D", "B, D, E, F"
                                      )
    m2["WT", "A"] <- 0.3
    m2["WT", "B"] <- 0.7
    m2["B", "B, C"] <- 0.2
    m2["B", "B, D"] <- 0.8
    m2["B, D", "B, D, E"] <- 0.4
    m2["B, D", "B, D, F"] <- 0.6
    m2 <- Matrix(m2, sparse = TRUE)

    expect_equivalent(paths_probs_2_df(trans_mat_2_paths_probs(m2),
                                       order = "path")[, "Prob"],
                 c(0.3, 0.7 * 0.2, 0.7 * 0.8 * 0.4, 0.7 * 0.8 * 0.6))
})


cat("\n Done test.rank-paths-graphs.R. Seconds = ",
    as.vector(difftime(Sys.time(), t1, units = "secs")), "\n")
