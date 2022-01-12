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

# #' @title Read output of HESBCN
# #'  
# #' Copied from hhttps://github.com/BIMIB-DISCo/PMCE/blob/main/Utilities/R/utils.R
# #' and extended to match the output with that of others CPMs
# #'
# #' @param file Filename with the output
# #' @param genes Custon genes names
# #' 
# #' @return A list with the adjacency matrix, the lambdas, the parent set
# #' and a data.frame with From-To edges and associated lambdas.
# import.hesbcn.evamtools <- function( file, genes = NULL ) {

#     # read results from file
#     results <- read.table(file = file,
#         header = FALSE, sep = "\n", 
#         stringsAsFactors = FALSE, check.names = FALSE)

#     # build adjacency matrix
#     start <- 2
#     end <- grep("accepted",results[,1])-1
#     adjacency_matrix <- NULL
#     for(entry in results[start:end,1]) {
#         adjacency_matrix = rbind(adjacency_matrix,
#             as.numeric(strsplit(entry,split = ",")[[1]]))
#     }

#     tmp_genes <- LETTERS[1:nrow(adjacency_matrix)]

#     rownames(adjacency_matrix) <- paste0(tmp_genes)
#     colnames(adjacency_matrix) <- paste0(tmp_genes)

#     # get lambda values and parent set types
#     lambdas <- results[(grep("Best Lambdas", results[, 1]) + 1), 1]
#     lambdas <- strsplit(lambdas,split = " ")[[1]][1:(length(strsplit(lambdas, split = " ")[[1]])-1)]
#     lambdas <- as.numeric(lambdas)
#     parent_set <- strsplit(results[(grep("best theta types",results[,1])+1),1],split=",")[[1]]
#     for (i in 1:length(parent_set)) {
#         if(parent_set[i] == "-1") {
#             parent_set[i] <- "Single"
#         }
#         else if(parent_set[i] == "0") {
#             parent_set[i] = "AND"
#         }
#         else if(parent_set[i] == "1") {
#             parent_set[i] <- "OR"
#         }
#         else if(parent_set[i] == "2") {
#             parent_set[i] <- "XOR"
#         }
#     }
#     # names(parent_set) <- paste0("Gene_", seq_len(adjacency_matrix))   
#     names(parent_set) <- tmp_genes

#     # make final results and normalize lambda values
#     adjacency_matrix <- cbind(rep(0, nrow(adjacency_matrix)), adjacency_matrix)
#     adjacency_matrix <- rbind(rep(0, ncol(adjacency_matrix)), adjacency_matrix)
#     rownames(adjacency_matrix)[1] <- "Root"
#     colnames(adjacency_matrix)[1] <- "Root"
#     adjacency_matrix["Root", as.numeric(which(colSums(adjacency_matrix) == 0))[-1]] = 1
#     lambdas_matrix <- adjacency_matrix
#     for (i in 2:nrow(lambdas_matrix)) {
#         curr_in_lambda <- lambdas[(i - 1)]
#         # we assume equally distributed lambdas among components 
#         # of a single formula
#         lambdas_matrix[as.numeric(which(adjacency_matrix[, i] == 1)), i] <-
#         curr_in_lambda / length(which(adjacency_matrix[, i] == 1))
#         # browser()
#     }

#     # set genes names
#     if (!is.null(genes)) {
#         rownames(adjacency_matrix) <- c("Root",genes)
#         colnames(adjacency_matrix) <- c("Root",genes)
#         rownames(lambdas_matrix) <- c("Root",genes)
#         colnames(lambdas_matrix) <- c("Root",genes)
#         names(parent_set) <- genes
#     }

#     # adjacency_matrix_2
#     # indexes_array <- data.frame(which(lambdas_matrix > 0, arr.ind = TRUE))
#     # indexes_list <- which(lambdas_matrix > 0, arr.ind = TRUE)
#     # lambdas <- lambdas_matrix[indexes_list]
#     # from <- rownames(lambdas_matrix)[indexes_array$row]
#     # to <- colnames(lambdas_matrix)[indexes_array$col]
#     # edges <- paste(from, to, sep = "->")
#     # adjacency_matrix_2 <- data.frame(From = from, To = to, Edge = edges, Lambdas = lambdas)

#     # estimated epsilon
#     epsilon <- as.numeric(gsub("Best epsilon = ", "", 
#         results[nrow(results), 1]))
    
#     # return results
#     return(list(adjacency_matrix = adjacency_matrix,
#         lambdas_matrix = lambdas_matrix,
#         parent_set = parent_set,
#         epsilon = epsilon
#         # ,
#         # edges = adjacency_matrix_2
#         )
#         )
# }


# # compute predictability of the HESBCN
# "compute.predictability" <- function( hesbcn ) {

#     # compute needed statistics to calculate predictability from graph
#     results = list()
#     graph = graph_from_adjacency_matrix(hesbcn$adjacency_matrix)
#     paths_predictability = NULL
#     number_path = 0
#     for(i in 1:nrow(hesbcn$adjacency_matrix)) {
#         # consider all leaves
#         if(all(hesbcn$adjacency_matrix[i,]==0)) {
#             leaf = rownames(hesbcn$adjacency_matrix)[i] # consider each leaf
#             paths = all_simple_paths(graph=graph,from="Root",to=leaf,mode="out")
#             paths_nodes = list()
#             transitions = list()
#             lambdas = list()
#             for(j in 1:length(paths)) {
#                 curr_path = names(paths[[j]]) # consider each simple path
#                 number_path = number_path + 1
#                 curr_nodes = NULL
#                 curr_transitions = NULL
#                 curr_lambdas = NULL
#                 curr_path_predictability = 1
#                 for(k in 2:length(curr_path)) { # consider each node visited in this path
#                     curr_node = curr_path[k]
#                     curr_nodes = c(curr_nodes,curr_node)
#                     curr_in_lambda = sum(hesbcn$lambdas_matrix[,curr_node])
#                     curr_out_lambda = sum(hesbcn$lambdas_matrix[curr_node,])
#                     if(k==length(curr_path)) {
#                         curr_out_lambda = 1
#                     }
#                     curr_transitions = c(curr_transitions,(curr_in_lambda/curr_out_lambda))
#                     curr_lambdas = c(curr_lambdas,curr_in_lambda)
#                     curr_path_predictability = curr_path_predictability * (curr_in_lambda/curr_out_lambda)
#                 }
#                 # save current results
#                 paths_nodes[[j]] = curr_nodes
#                 transitions[[j]] = curr_transitions
#                 lambdas[[j]] = curr_lambdas
#                 paths_predictability = c(paths_predictability,curr_path_predictability)
#             }
#             results[[paste0("Evolution_",(length(results)+1))]] = list(leaf_node=leaf,simple_paths=paths_nodes,transitions_vector=transitions,lambdas=lambdas)
#         }
#     }
#     paths_predictability = paths_predictability / sum(paths_predictability) # normalize path predictabilities to sum to 1
#     normalized_paths_predictability = paths_predictability

#     # compute graph metrics
#     graph_predictability = 0.0
#     for(path in normalized_paths_predictability) {
#         graph_predictability = graph_predictability + (path*log(path))
#     }
#     graph_entropy = - graph_predictability
#     max_entropy = - log((1/number_path))
#     graph_predictability = 1 - (graph_entropy/max_entropy)

#     # return results
#     predictability = list(statistics=results,paths_probabilities=normalized_paths_predictability,graph_entropy=graph_entropy,number_path=number_path,graph_predictability=graph_predictability)
#     return(predictability)

# }

# # enumerate valid genotypes given the HESBCN
# "valid.genotypes" <- function( hesbcn, predictability ) {

#     # consider all valid genotypes
#     valid_genotypes = NULL
#     template_genotype = rep("*",(nrow(hesbcn$adjacency_matrix)-1))
#     names(template_genotype) = rownames(hesbcn$adjacency_matrix)[2:nrow(hesbcn$adjacency_matrix)]
#     for(i in 1:length(predictability$statistics)) { # consider each independent evolution
#         for(j in 1:length(predictability$statistics[[i]][["simple_paths"]])) { # consider each path within the current evolution
#             curr_genotype = template_genotype # start from empty/template genotype
#             path_time = 0.00
#             for(k in 1:length(predictability$statistics[[i]][["simple_paths"]][[j]])) {
#                 curr_gene = predictability$statistics[[i]][["simple_paths"]][[j]][[k]]
#                 curr_genotype[curr_gene] = 1
#                 if(hesbcn$parent_set[curr_gene]=="AND") {
#                     curr_parents_extra_nodes = names(which(hesbcn$adjacency_matrix[,curr_gene]==1)) # consider all parents of current node
#                     curr_genotype[curr_parents_extra_nodes] = 1 # AND --> all parents are 1
#                 }
#                 if(hesbcn$parent_set[curr_gene]=="XOR") {
#                     curr_parents_extra_nodes = names(which(hesbcn$adjacency_matrix[,curr_gene]==1)) # consider all parents of current node
#                     curr_parents_extra_nodes = curr_parents_extra_nodes[which(!curr_parents_extra_nodes==predictability$statistics[[i]][["simple_paths"]][[j]][(k-1)])]
#                     curr_genotype[curr_parents_extra_nodes] = 0 # XOR --> only one parent is 1
#                 }
#                 path_time = path_time + (1/predictability$statistics[[i]][["lambdas"]][[j]][[k]])
#                 valid_genotypes = rbind(valid_genotypes,c(curr_genotype,path_time))
#             }
#         }
#     }
#     rownames(valid_genotypes) = paste0("Genotype_",1:nrow(valid_genotypes))
#     colnames(valid_genotypes)[ncol(valid_genotypes)] = "Progression_Time"
    
#     # return results
#     return(valid_genotypes)

# }

# # estimate patients' corrected genotypes from the HESBCN
# "corrected.genotypes" <- function( genotypes, patients, epsilon ) {

#     # estimate patients' corrected genotypes
#     corrected_genotypes = NULL
#     for(i in 1:nrow(patients)) {
#         # compute likelihood of patients' attachment to corrected genotypes
#         patients_attachments_likelihoods = NULL
#         curr_patient = as.character(patients[i,])
#         for(j in 1:nrow(genotypes)) {
#             curr_genotype = as.character(genotypes[j,1:(ncol(genotypes)-1)])
#             consistent = length(which(curr_patient[which(curr_genotype!="*")]==curr_genotype[which(curr_genotype!="*")]))
#             inconsistent = length(which(curr_patient[which(curr_genotype!="*")]!=curr_genotype[which(curr_genotype!="*")]))
#             wild_card = which(curr_genotype=="*")
#             wild_card_p = 1
#             if(length(wild_card)>0) {
#                 for(k in wild_card) {
#                     if(as.numeric(curr_patient[k])==0) {
#                         curr_prob = length(which(patients[,k]==0))/nrow(patients)
#                     }
#                     if(as.numeric(curr_patient[k])==1) {
#                         curr_prob = length(which(patients[,k]==1))/nrow(patients)
#                     }
#                     wild_card_p = wild_card_p * ((curr_prob*(1-epsilon))+((1-curr_prob)*epsilon))
#                 }
#             }
#             curr_attachment_likelihood = ((1-epsilon)^consistent) * (epsilon^inconsistent) * wild_card_p
#             patients_attachments_likelihoods = c(patients_attachments_likelihoods,curr_attachment_likelihood)
#         }
#         max_likelihood_attachment = which(patients_attachments_likelihoods==max(patients_attachments_likelihoods))[1]
#         curr_corrected_genotype = as.character(genotypes[max_likelihood_attachment,1:(ncol(genotypes)-1)])
#         corrected_genotypes = rbind(corrected_genotypes,curr_corrected_genotype)
#     }
#     rownames(corrected_genotypes) = rownames(patients)
#     colnames(corrected_genotypes) = colnames(genotypes)[1:(ncol(genotypes)-1)]
    
#     # return results
#     return(corrected_genotypes)

# }

# generate_trans_matrix <- function(data, parameter_column_name){

#   tmp <- try(df_2_access_genots_and_graph_OR(data[, c("From", "To")]))
  
#   if(inherits(tmp, "try-error")) {
#     stop("how is this happening? there was edges component!")
#   } else {
#     accessible_genots <- tmp$accessible_genots
#     fgraph <- unrestricted_fitness_graph_sparseM(accessible_genots)
#   }
  
#   if (is.null(data[,parameter_column_name])){
#     stop("No such column")
#   }
  
#   weights <- unique(data[, c("To", parameter_column_name)])
#   rownames(weights) <- weights[, "To"]
#   print(weights)
#   print(typeof(fgraph))
#   print(typeof(weights))
#   weighted_fgraph <- transition_fg_sparseM(fgraph, weights)

#   return(weighted_fgraph)
# }