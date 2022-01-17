source("../External-code/HESBCN/hesbcn-external.R")
## TODO
# Create something similar to cbn-process for checking if the package is installed


#' @title Run HESBCN
#' 
#' @param data Cross secitonal data. Matrix of genes (columns)
#' and individuals (rows)
#' @param n_steps Number of steps to run. Default: 100000
#' @param tmp_folder Folder name where the oput is located. 
#' It will be place under /tmp/HESBCN/tmp_folder
#' @param seed Seed to run the experiment
#' @param clean_dir Whether to delete the folder upon completion
#' 
#' @return A list with the adjacency matrix, the lambdas, the parent set
#' and a data.frame with From-To edges and associated lambdas.
do_HESBCN <- function(data, n_steps=100000, 
    tmp_folder="", seed=NULL, clean_dir=FALSE){
    # Setting tmp folder

    date_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    random_letters <- paste(c("_",
        tmp_folder,
        "_",
        LETTERS[floor(runif(4, min = 1, max = 26))]), collapse = "")
    tmp_folder <- paste(c(date_time, random_letters), collapse = "")

    tmp_folder <- file.path("/", "tmp", "HESBCN", tmp_folder)
    dir.create(tmp_folder, recursive = TRUE)
    orig_folder <- getwd()

    setwd(tmp_folder)

    orig_gene_names <- colnames(data)
    colnames(data) <- LETTERS[1:ncol(data)]

    write.csv(data, "input.txt", row.names = FALSE, quote = FALSE)
 
    # Launching
    print("Running HESBCN")
    if (is.null(seed)) {
        command <- sprintf("h-esbcn -d input.txt -o output.txt -n %1.f", 
            n_steps)
    } else if (is.numeric(seed) & seed > 0) {
        command <- sprintf("h-esbcn -d input.txt -o output.txt -n %1.f -s %1.f", 
            n_steps, seed)
    }
    system(command, ignore.stdout = TRUE)

    # Reading output
    model_info <- import.hesbcn("output.txt", genes = orig_gene_names)

    # Updating gene names
    # gene_names <- orig_gene_names
    # names(model_info$parent_set) <- gene_names
    # rownames(model_info$lambdas_matrix) <-
    # colnames(model_info$lambdas_matrix) <-
    # rownames(model_info$adjacency_matrix) <-
    # colnames(model_info$adjacency_matrix) <- c("Root", gene_names)

    indexes_array <- data.frame(which(model_info$lambdas_matrix > 0, arr.ind = TRUE))
    indexes_list <- which(model_info$lambdas_matrix > 0, arr.ind = TRUE)
    lambdas <- model_info$lambdas_matrix[indexes_list]
    from <- rownames(model_info$lambdas_matrix)[indexes_array$row]
    to <- colnames(model_info$lambdas_matrix)[indexes_array$col]
    edges <- paste(from, to, sep = "->")
    adjacency_matrix_2 <- data.frame(From = from, To = to, Edge = edges, Lambdas = lambdas)
    model_info$edges <- adjacency_matrix_2


    # Transforming data from model
    # weighted_fgraph <- generate_trans_matrix(model_info$hesbcn_out, "Lambdas")

    # ##TODO: include thetas for out-of-the-path mutations?
    # trans_mat_genots <- rowScaleMatrix(weighted_fgraph)

    # Housekeeping
    setwd(orig_folder)
    if (clean_dir) {
        unlink(tmp_folder, recursive = TRUE)
    }

    return(model_info)
}


