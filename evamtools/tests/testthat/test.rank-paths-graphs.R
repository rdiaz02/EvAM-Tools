
test_that("Returns the correct number and type of vertex labels", {

    adj_matrix <- matrix(0, nrow = 6, ncol = 6)
    rownames(adj_matrix) <- colnames(adj_matrix) <- c("WT", "A", "A, B",
                                                      "A, C", "C", "C, D")
    adj_matrix["WT", "A"] <- 1000
    adj_matrix["WT", "C"] <- 500
    adj_matrix["A", "A, B"] <- 500
    adj_matrix["A", "A, C"] <- 200
    adj_matrix["C", "A, C"] <- 100
    adj_matrix["C", "C, D"] <- 200

    g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)
    
    all_genotypes <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    all_paths <- rank_paths(g)

    l1 <- c("WT", "A", "A, B", "", "", "")
    names(l1) <- all_genotypes

    l2 <- c("WT", "A", "A, B", "A, C", "", "")
    names(l2) <- all_genotypes

    l3 <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    names(l3) <- all_genotypes

    l4 <- c("WT", "A", "A, B", "A, C", "C", "C, D")
    names(l4) <- all_genotypes



    labels_1 <- compute_vertex_labels(g, all_paths, top_paths = 1)$vertex_labels
    labels_2 <- compute_vertex_labels(g, all_paths, top_paths = 2)$vertex_labels
    labels_3 <- compute_vertex_labels(g, all_paths, top_paths = 3)$vertex_labels
    labels_4 <- compute_vertex_labels(g, all_paths, top_paths = 4)$vertex_labels
    labels_40 <- compute_vertex_labels(g, all_paths, top_paths = 40)$vertex_labels
    labels_null <- compute_vertex_labels(g, all_paths, top_paths = NULL)$vertex_labels

    expect_equal(labels_1, l1)
    expect_equal(labels_2, l2)
    expect_equal(labels_3, l3)
    expect_equal(labels_4, l4)
    expect_equal(labels_4, labels_40)
    expect_equal(labels_4, labels_null)
})

cat("\n Done test.rank-paths-graphs.R \n")
