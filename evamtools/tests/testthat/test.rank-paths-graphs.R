#' Test the procedure of node labeling for the plot_genot_fg function
#' It labels nodes as those involved in the most transited paths
#' that star in the root of the tree and that lead to any leave 

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
    all_paths <- rank_paths(g)

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

    labels_0 <- compute_vertex_labels(g, all_paths, top_paths = 0)
    labels_1 <- compute_vertex_labels(g, all_paths, top_paths = 1)
    labels_2 <- compute_vertex_labels(g, all_paths, top_paths = 2)
    labels_3 <- compute_vertex_labels(g, all_paths, top_paths = 3)
    labels_4 <- compute_vertex_labels(g, all_paths, top_paths = 4)
    labels_40 <- compute_vertex_labels(g, all_paths, top_paths = 40)
    labels_null <- compute_vertex_labels(g, all_paths, top_paths = NULL)
    
    expect_equal(labels_1$vertex_labels, l1)
    expect_equal(labels_2$vertex_labels, l2)
    expect_equal(labels_3$vertex_labels, l3)
    expect_equal(labels_4$vertex_labels, l4)
    expect_equal(labels_4$vertex_labels, labels_40$vertex_labels)
    expect_equal(labels_4$vertex_labels, labels_null$vertex_labels)
    expect_equal(labels_4$vertex_labels, labels_0$vertex_labels)

    expect_equal(labels_1$adj_matrix, adj1, check.attributes=FALSE)
    expect_equal(labels_2$adj_matrix, adj2, check.attributes=FALSE)
    expect_equal(labels_3$adj_matrix, adj3, check.attributes=FALSE)
    expect_equal(labels_4$adj_matrix, adj4, check.attributes=FALSE)
    expect_equal(labels_40$adj_matrix, adj4, check.attributes=FALSE)
    expect_equal(labels_null$adj_matrix, adj4, check.attributes=FALSE)
    expect_equal(labels_0$adj_matrix, adj4, check.attributes=FALSE)
})

cat("\n Done test.rank-paths-graphs.R \n")
