## Copyright 2022 Pablo Herrera Nieto, Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU Affero General Public License (AGPLv3.0) as published by
## the Free Software Foundation, either version 3 of the License, or (at your
## option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## You should have received a copy of the GNU Affero General Public License along
## with this program.  If not, see <http://www.gnu.org/licenses/>.

## Given the Genotype, Count data.frame, return it in standard order.
## If single row, or empty, do nothing.
reorder_to_standard_order_count_df <- function(x) {
    x <- x[x$Counts > 0, ]
    if (nrow(x) <= 1) return(x)
    counts_tmp <- x$Counts
    names(counts_tmp) <- x$Genotype
    counts_tmp <- reorder_to_standard_order(counts_tmp)
    ret_tmp <- na.omit(data.frame(Genotype = names(counts_tmp),
                                  Counts = counts_tmp))
    attributes(ret_tmp)$na.action <- NULL
    stopifnot(nrow(ret_tmp) == nrow(x))
    ## rownames(ret_tmp) <- 1:nrow(ret_tmp)
    return(ret_tmp)
}

## Given a data.frame that contains column Genotype, return it in standard order.
## This could be done much more efficiently (we could just put the genotypes
## in standard order, instead of creating all, etc). 
reorder_to_standard_order_arbitrary_df <- function(x) {
    if (nrow(x) <= 1) return(x)
    counts_tmp <- x$Genotype
    names(counts_tmp) <- x$Genotype
    counts_tmp <- na.omit(reorder_to_standard_order(counts_tmp))
    ret_tmp <- data.frame(Genotype = names(counts_tmp), x)
    stopifnot(isTRUE(all(ret_tmp$Genotype == ret_tmp$Genotype.1)))
    rm_cols <- which(colnames(ret_tmp) == "Genotype.1")
    ret_tmp <- ret_tmp[, -rm_cols]
    ## There might, or might not, be an Index column. We create it/overwrite it
    ret_tmp$Index <- seq_len(nrow(ret_tmp))
    g_index <- which(colnames(ret_tmp) == "Genotype")
    i_index <- which(colnames(ret_tmp) == "Index")
    col_order <- c(i_index, g_index, (1:ncol(ret_tmp))[-c(i_index, g_index)])
    ret_tmp <- ret_tmp[, col_order]
    stopifnot(nrow(ret_tmp) == nrow(x))
    return(ret_tmp)
}




## From one-row-per-subject, as matrix of 0/1, to grouped counts of genotypes
get_csd <- function(complete_csd,
                    default_counts = .ev_SHINY_dflt$template_data$csd_counts){
    if (is.null(complete_csd)) return(default_counts)
    csd <- data_to_counts(complete_csd, out = "data.frame", omit_0 = TRUE)
    rownames(csd) <- csd$Genotype
    return(csd)
}

modify_dag <- function(dag, from_node, to_node, operation, parent_set,
                       dag_model="HESBCN",
                       default_dag = .ev_SHINY_dflt$template_data$dag,
                       default_DAG_parent_set = .ev_SHINY_dflt$template_data$DAG_parent_set){
    value_from_operation <- c(1, 0)
    names(value_from_operation) <- c("add", "remove")

    if (operation == "clear") {
        return(list(dag = default_dag,
                    DAG_parent_set = default_DAG_parent_set)) 
    }

    if (is.null(from_node) | is.null(to_node) | is.null(dag)){
        stop("From and To options and DAG have to be defined.")
    } 


    if(!all(names(parent_set) == colnames(dag)[-1])){
        stop("Dag and parent set should have information of the same genes.")
    }
    
    ## from_node <- ifelse(from_node == "Root", "WT", from_node)
    if (!(from_node %in% colnames(dag)) | !(to_node %in% colnames(dag))) {
        stop("Both From and To options have to be valid gene names.")
    } else if(!all(unique(as.vector(dag)) %in% c(0, 1))){
        stop("The adjacency matrix should only contain 0 and 1.")
    } else if (from_node == to_node){
        stop("Both From and To options must be different.")
    } else if(dag[to_node, from_node] == 1 && operation == "add") {
        stop("Relationships cannot be bidirectional.")
    } else if (dag[from_node, to_node] == 1 && operation == "add") {
        stop("That edge is already present.")
    } else if (dag["Root", to_node] == 1 && operation == "add") {
        stop("A direct children of Root cannot have multiple parents.")
    } 

    tmp_dag <- dag
    tmp_dag[from_node, to_node] = value_from_operation[operation]

    ##Taking subcomponent starting from WT
    g1 <- igraph::graph_from_adjacency_matrix(tmp_dag, mode = "directed")

    for(tmp_g in igraph::decompose(g1)){
        ## if ("WT" %in% V(tmp_g)$name) {
        if ("Root" %in% V(tmp_g)$name) {          
            tmp_vertices <- igraph::as_data_frame(tmp_g)
        }
    }
    
    tmp_dag2 <- default_dag*0
    colnames(tmp_dag2) <- rownames(tmp_dag2) <- colnames(dag) 
    if(nrow(tmp_vertices>0)){
        for(i in 1:nrow(tmp_vertices)){
            tmp_row <- tmp_vertices[i, ]
            tmp_dag2[tmp_row$from, tmp_row$to] <- 1
        }
    }

    ## Restructure the DAG
    number_of_parents <- colSums(tmp_dag2)
    number_of_children <- rowSums(tmp_dag2)
    for(i in colnames(tmp_dag2)[-1]){
        if(number_of_children[i] > 0 & number_of_parents[i] == 0){
            ## We add link to Root
            tmp_dag2["Root", i] <- 1
        }
    }

    if(all(dag == tmp_dag2)) {
        stop("This operation had no effect. ",
             "Some possible causes and solutions: ",
             "\n\n                          ",
             "A) You tried to draw a disconnected DAG, ",
             "where the From node is not already ",
             "part of the DAG (i.e., has no ancestor): ",
             "all nodes except Root must have an ancestor.  ",
             "\n\n                            ",
             "B) You tried to remove the edge from Y to X, ",
             "where X has descendants that have other ancestors: ",
             "This operation is sometimes not possible,  ",
             "especially if Y is Root. ",
             "Remove first the edges from X to its descendants.  ",
             "\n\n                             ",           
             "C) You tried to remove a non-existent edge. ",
             "\n\n",
             "D) Something else. Start removing from the leave ",
             "nodes."
             )
    }
    g2 <- igraph::graph_from_adjacency_matrix(tmp_dag2, mode = "directed")
    if (dag_model %in% c("HESBCN", "OncoBN") && !igraph::is_dag(g2)) {
        stop("This relationship breaks the DAG. Revise it. ",
             "For example, did you try to introduce a cycle?")
    } else if (dag_model %in% c("OT") && any(number_of_parents > 1)) {
        stop("This operation does not give a tree. ",
             "With OT nodes cannot have multiple parents.")
    }


    ## Recompute parent set
    number_of_parents <- colSums(tmp_dag2)[-1]
    tmp_parent_set <- parent_set
    tmp_parent_set[number_of_parents <= 1] <- "Single" ##Default is Single
    tmp_parent_set[number_of_parents > 1 & !(parent_set %in% c("AND", "OR", "XOR"))] <- "AND" ##Default is AND

    if (dag_model == "OncoBN")  {
        if (length(unique(tmp_parent_set)) > 2) {
            stop("This operation does not give a valid OncoBN model. ",
                 "The OncoBN model can only include ",
                 "one type of relationship ",
                 "(conjunctive ---AND--- or disjunctive ---OR---) ",
                 "but not both.")
        } else if (any(tmp_parent_set == "XOR" )) {
            ## I think this is impossible from here, but just in case
            stop("This operation does not give a valid OncoBN model. ",
                 "The OncoBN model cannot include ",
                 "XOR relationships.")
        }
    }
    
    
    return(list(dag = tmp_dag2, parent_set = tmp_parent_set))
}

modify_lambdas_and_parent_set_from_table <- function(dag_data, info,
                                                     lambdas, dag, parent_set,
                                                     dag_model) {
    if (!all(names(lambdas) == names(parent_set))) {
        stop("Lambdas and parent set should have information about the same genes")
    }
    if (!all(names(lambdas) == colnames(dag)[-1])) {
        stop("DAG should have information about the same genes as parent set")
    }

    ## Different lambdas
    if (dag_model %in% c("OncoBN", "HESBCN")) {
        new_lambdas <- suppressWarnings(as.numeric(info[info["col"] == 3, "value"]))
    } else if (dag_model == "OT") {
        new_lambdas <- suppressWarnings(as.numeric(info[info["col"] == 2, "value"]))
    }

    if (any(is.na(new_lambdas))) {
        stop("There are missing lambdas")
    }
    if (dag_model == "HESBCN") {
        old_lambdas <- dag_data$Lambdas
    } else if (dag_model == "OT") {
        old_lambdas <- dag_data$Weight
    } else if (dag_model == "OncoBN") {
        old_lambdas <- dag_data$theta
    }
    
    all_genes <- dag_data$To

    if (length(new_lambdas) &&
        (length(new_lambdas) != length(old_lambdas)))
        stop("Inconsistent state of DAG structure after editing. ",
             "\n\n This could be the result of addition/removals. ",
             "Try adding/removing genes one by one, using ",
             "'Add edge' and 'Remove edge' instead of passing ",
             "0s to lambdas/weights. "
             )

    if (length(new_lambdas) == 0) {
        changed_genes <- character()
    } else {
        ## yes, the if does the same, but an if with a logical(0) sucks
        changed_genes <- all_genes[(new_lambdas != old_lambdas)
                                   & (new_lambdas > 0)]
    }
    info_from <- info[info["col"] == 0,"value"]
    info_to <- info[info["col"] == 1,"value"]
    if (!(all(c(info_from, info_to, dag_data$To, dag_data$From) %in%
              c("Root", names(lambdas))))) {
        stop("There are unkown genes")
    }

    if (length(new_lambdas) == 0) {
        changed_lambdas <- numeric()
    } else {  
        changed_lambdas <- new_lambdas[(new_lambdas != old_lambdas)
                                       & (new_lambdas > 0)]
    }
    
    tmp_lambdas <- lambdas
    tmp_lambdas[changed_genes] <- changed_lambdas

    if (dag_model %in% c("OT", "OncoBN")
        & (any(tmp_lambdas <= 0) | any(tmp_lambdas >= 0.99995))){
        stop("thetas/probabilities should be between 0 and 1.")
    }

    ##Relationships
    number_of_parents <- colSums(dag)[-1]
    if (dag_model %in% c("HESBCN", "OncoBN")) {
        tmp_parent_set <- parent_set
        number_of_parents <- number_of_parents[number_of_parents > 0]
        new_relationships <- info[info["col"] == 2, "value"]
        names(new_relationships) <- info[info["col"] == 1, "value"]

        genes_in_relationships <- names(new_relationships)
        freq_genes_in_relationships <- table(genes_in_relationships)

        old_relationships <- stats::setNames(dag_data$Relation, dag_data$To)
        if (dag_model %in% c("HESBCN"))  
            new_relationships[!(new_relationships %in%
                                c("Single", "AND", "OR", "XOR"))] <- "AND"
        if (dag_model %in% c("OncoBN"))  
            new_relationships[!(new_relationships %in%
                                c("Single", "AND", "OR"))] <- "OR"
        
        new_relationships[names(number_of_parents[number_of_parents < 2])] <- "Single"
        new_relationships[setdiff(unique(genes_in_relationships), 
                                  names(freq_genes_in_relationships[freq_genes_in_relationships > 1]))] <- "Single"
        
        multiple_parents <- names(number_of_parents[number_of_parents > 1])

        for (i in multiple_parents){
            changed_relationships <-
                setdiff(new_relationships[genes_in_relationships == i], c(parent_set[i]))
            if (length(changed_relationships) > 0){
                tmp_parent_set[i] <- changed_relationships[[1]]
            }
        }
        ## consistency check
        if (dag_model == "OncoBN")  {
            if (length(unique(tmp_parent_set)) > 2) {
                stop("This operation does not give a valid OncoBN model. ",
                     "The OncoBN model can only include ",
                     "one type of relationship ",
                     "(conjunctive ---AND--- or disjunctive ---OR---) ",
                     "but not both.")
            } else if (any(tmp_parent_set == "XOR" )) {
                ## I think this is impossible from here, but just in case
                stop("This operation does not give a valid OncoBN model. ",
                     "The OncoBN model cannot include ",
                     "XOR relationships.")
            }
        }        
    } else if (dag_model == "OT") {
        tmp_parent_set <- parent_set
        tmp_parent_set[1:length(tmp_parent_set)] <- "Single"
    } 
    

    return(list(lambdas = tmp_lambdas, parent_set = tmp_parent_set))
}

get_mhn_data <- function(thetas, noise = 0, N = 10000){
    if (any(is.null(thetas))) stop("thetas should be defined")
    n_genes <- ncol(thetas)
    gene_names <- colnames(thetas)

    mhn_trm <- MHN_from_thetas(thetas)$MHN_trans_rate_mat
    mhn_probs <- probs_from_trm(mhn_trm)
    tmp_samples_as_vector <- genot_probs_2_pD_ordered_sample(x = mhn_probs,
                                                             ngenes = n_genes,
                                                             gene_names = gene_names,
                                                             N = N,
                                                             out = "vector"
                                                             )
    data_with_noise <- genotypeCounts_to_data(tmp_samples_as_vector,
                                              e = noise)
    csd_counts <- data_to_counts(data_with_noise, out="data.frame")
    return(list(csd_counts = csd_counts,
                data = data_with_noise))
}

generate_sample_from_dag <- function(data, parent_set, noise = 0, N = 10000,
                                     dag_model = "HESBCN", epos = 0.01) {

    if (nrow(data) == 0) stop("The DAG does not contain any edge.")
    if (any(is.null(data))) stop("Data should be defined")

    gene_names <- names(parent_set)
    n_genes <- length(parent_set)
    
    if (dag_model == "HESBCN") {
        dag_trm <- HESBCN_model_2_output(data, parent_set)$HESBCN_trans_rate_mat
        dag_probs <- probs_from_trm(dag_trm)
    } else if (dag_model == "OT") {
        data$OT_edgeWeight <- data$Weight
        dag_probs <- OT_model_2_output(data, epos)$OT_predicted_genotype_freqs
    } else if (dag_model == "OncoBN") {
        epsilon <- epos
        dag_probs <- OncoBN_model_2_output(data, epsilon)$OncoBN_predicted_genotype_freqs
    }

    tmp_samples_as_vector <- genot_probs_2_pD_ordered_sample(x = dag_probs,
                                                             ngenes = n_genes,
                                                             gene_names = gene_names,
                                                             N = N,
                                                             out = "vector"
                                                             )
    data_with_noise <- genotypeCounts_to_data(tmp_samples_as_vector,
                                              e = noise)
    csd_counts <- data_to_counts(data_with_noise, out="data.frame")
    return(list(csd_counts = csd_counts,
                data = data_with_noise))
}

create_tabular_data <- function(data) {
  available_methods <- c("OT", "OncoBN", "CBN", "MHN", "HESBCN", "MCCBN",
                         "HyperTraPS"
                         ##, "BML" ## Nope, BML is not here
                         )
    attr_to_make_tabular <- c("trans_mat", "trans_rate_mat",
                              "predicted_genotype_freqs"
                              )
    for(i in names(data)){
        if(grepl("sampled_genotype_counts", i)
           && !(is.null(data[i]))){
            attr_to_make_tabular <- c(attr_to_make_tabular
                                    , "sampled_genotype_counts"
                                      ## , "obs_genotype_transitions"
                                      )
        }
    }

    tabular_data <- list()
    for (attr in attr_to_make_tabular){
        if (attr %in% c("sampled_genotype_counts",
                        "predicted_genotype_freqs")) {
            all_counts <- data.frame(Genotype = c(NULL))

            for (method in available_methods){
                tmp_data <- data[[paste0(method, "_", attr)]]
                if (!any(is.null(tmp_data)) && !any(is.na(tmp_data))) {
                    if (any(is.null(all_counts$Genotype))){
                        all_counts <- data.frame(Genotype = names(tmp_data)) #They include all genotypes
                        rownames(all_counts) <- all_counts$Genotype
                    }
                    ## To avoid relying in indexes
                    all_counts[, method] <- rep(0, nrow(all_counts))
                    tmp_data <- round(tmp_data, 3)
                    for (genotype in names(tmp_data)){
                        all_counts[genotype, method] <- tmp_data[genotype]
                    }
                }
            }
            all_counts <-
                data.frame(Index = standard_rank_genots_1(all_counts$Genotype),
                           all_counts)
            tabular_data[[attr]] <- all_counts[order(all_counts$Index), ]

        } else if (attr %in% c("trans_rate_mat",
                               ## "obs_genotype_transitions",
                               "trans_mat")) {
            df <- data.frame(From = character(),
                             To = character(), 
                             OT = numeric(), 
                             OncoBN = numeric(), 
                             CBN = numeric(), 
                             HESBCN = numeric(), 
                             MCCBN = numeric(), 
                             MHN = numeric(),
                             HyperTraPS = numeric(),
                       ## No, BML is never here
                       ## BML = numeric(),
                             stringsAsFactors = FALSE) 

            for (method in available_methods){
                                        #1 Matrix to vector
                tmp_data <- data[[paste0(method, "_", attr)]]
                if (!is.null(tmp_data)) {
                    indexes <- which(as.matrix(tmp_data) > 0, arr.ind = TRUE)
                    genotypes <- colnames(tmp_data)
                    transitions <- mapply(function(x, y) paste0(x, " -> ", y)
                                        , genotypes[indexes[, "row"]], genotypes[indexes[, "col"]])
                                        # counts <- tmp_data[tmp_data > 0]
                    counts <- mapply(function(x, y) tmp_data[x, y]
                                   , genotypes[indexes[, "row"]], genotypes[indexes[, "col"]])
                    names(counts) <- transitions

                                        #2 Adding empty row for each transitions
                    for(transition in transitions){
                        ## New row
                        if(all(is.na(df[transition, ]))){
                            genes <- strsplit(transition, " -> ")[[1]]
                            df[transition, ] <- c(genes[[1]], genes[[2]],
                                                  rep(0, length(available_methods)))
                            df[transition, method] <- counts[[transition]]
                        } else {
                            df[transition, method] <- counts[[transition]]
                        }
                    }
                }
            }
            ## Dropping empty columns
            for (method in available_methods){
                tmp_col <- as.numeric(df[[method]])
                if (sum(abs(tmp_col)) == 0) df[[method]] <- NULL
                else df[[method]] <- round(tmp_col, 3)
            }

            df <- data.frame(Index = standard_rank_genots_2(df$From, df$To),
                             df)
            ## Sorting 
            ## order_by_counts <- order(rowSums(df[,3:ncol(df)]),
            ##                          decreasing = TRUE)

            ## tabular_data[[attr]] <- df[order_by_counts, ]
            tabular_data[[attr]] <- df[order(df$Index), ]
        }
    }
    return(tabular_data)
}

## FIXME: try to simplify and clarify this function.
## Fills empty attributes of a csd data object

## Yes, could have passed all default_* as part of a single list. Whatever.
to_stnd_csd_dataset <- function(data,
                                default_max_genes = .ev_SHINY_dflt$max_genes,
                                default_template_data = .ev_SHINY_dflt$template_data
                                ) {
    ## Be completely explicit upfront!!!
    default_data <- default_template_data$data
    default_lambdas <- default_template_data$lambdas
    default_DAG_parent_set <- default_template_data$DAG_parent_set
    default_dag <- default_template_data$dag
    default_thetas <- default_template_data$thetas
    default_this_d_dag_model <- "HESBCN"
    
    if(is.null(data)) return(default_template_data)

    if(typeof(data) != "list"){
        stop("Input data should be a list")
    }
    
    new_data <- list()
    
    if(!is.null(data$gene_names)){
        tmp_names <- data$gene_names
    } else if(!is.null(colnames(data$data))){
        tmp_names <- colnames(data$data)
    } else if(!is.null(colnames(data$dag))){
        tmp_names <- colnames(data$dag)[-1]
    } else if(!is.null(colnames(data$thetas))){
        tmp_names <- colnames(data$thetas)
    } else {
        tmp_names <- c()
    }
    new_data$gene_names <-
        c(tmp_names,
          LETTERS[(length(tmp_names) + 1) : default_max_genes ])[1:default_max_genes]
    
    new_data$name <- data$name  

    if(is.null(data$lambdas)) {
        new_data$lambdas <- default_lambdas
    }else{
        if(!is.numeric(data$lambdas)){
            stop("Lambdas should only contain numbers")
        }
        new_lambdas <- default_lambdas
        new_lambdas[1:length(data$lambdas)] <- data$lambdas
        new_data$lambdas <- new_lambdas
    }
    names(new_data$lambdas) <- new_data$gene_names

    if(is.null(data$DAG_parent_set)) {
        new_data$DAG_parent_set <- default_DAG_parent_set
    }else {
        if(!all(data$DAG_parent_set %in% c("Single", "AND", "OR", "XOR")))
            stop("Parent set must include only 'Single', 'AND', 'OR' or 'XOR'")
        new_parent_set <- default_DAG_parent_set
        new_parent_set[1:length(data$DAG_parent_set)] <- data$DAG_parent_set
        new_data$DAG_parent_set <- new_parent_set
    }
    names(new_data$DAG_parent_set ) <- new_data$gene_names

    if(is.null(data$dag)) {
        new_data$dag <- default_dag
    } else {
        if(!all(unique(data$dag) %in% c(0, 1))){
            stop("Adjacency matrix should be binary: only 0 and 1")
        }
        to_keep <- which(colSums(data$dag)>0 | rowSums(data$dag)>0)
        tmp_dag <- data$dag[to_keep, to_keep]
        n_genes <- ncol(data$dag)
        new_dag <- default_dag 
        new_dag[to_keep, to_keep] <- tmp_dag
        new_data$dag <- new_dag

        if(!is_dag(graph_from_adjacency_matrix(tmp_dag, mode="directed"))){
            stop("The graph is not a DAG")
        }

        ## Revising parent set
        ## Only genes with multiple parent can have something different from "Single" relationship
        new_data$DAG_parent_set[colSums(new_dag)[-1] <= 1] <- "Single"
    }

    rownames(new_data$dag) <- colnames(new_data$dag) <- c("Root", new_data$gene_names)

    if(is.null(data$thetas)) {
        new_data$thetas <- default_thetas
    } else {
        if(!is.numeric(data$thetas)){
            stop("Theta matrix should only contain numbers")
        }
        n_genes <- ncol(data$thetas)
        new_thetas <- default_thetas 
        new_thetas[1:n_genes, 1:n_genes] <- data$thetas
        new_data$thetas <- new_thetas
    }
    rownames(new_data$thetas) <- colnames(new_data$thetas) <- new_data$gene_names

    if (is.null(data$data)) {
        new_data$data <- default_data
    } else {
        if(!all(unique(unlist(data$data)) %in% c(0, 1))){
            stop("Data should be binary: only 0 and 1")
        }
        new_data$data <- data$data
        colnames(new_data$data) <- new_data$gene_names[1:ncol(new_data$data)]
    }
    new_data$csd_counts <- get_csd(new_data$data)
    if (is.null(data$this_d_dag_model)) {
        new_data$this_d_dag_model <- default_this_d_dag_model
    } else {
        new_data$this_d_dag_model <- data$this_d_dag_model
    }
    return(new_data)
}

to_stnd_csd_all_datasets <- function(datasets){
    all_new_data <- list()
    for(i in c("upload", "csd", "dag", "matrix")){
        tmp_data <- datasets[[i]] 
        for(j in names(tmp_data)) all_new_data[[i]][[j]] <-     
                                      to_stnd_csd_dataset(tmp_data[[j]])
    }
    return(all_new_data)
}



## ## No longer used. See A21_gnn_numfix. So in the future this might be useful
## ## if we want to be more permissive. Thus DO NOT DELETE!
## get_display_freqs <- function(freqs, n_genes, gene_names, input2build) {
##     if (input2build %in% c("upload", "dag", "matrix"))
##         stop("get_display_freqs should not be called when ",
##              "input is either dag or upload or matrix (mhn); ",
##              "DAG has a different procedure for dealing with number of genes ",
##              "and upload does not use that setting at all.")
##     if (is.null(freqs)) return(.ev_SHINY_dflt$template_data$csd_counts)
##     if (nrow(freqs) == 0) return(.ev_SHINY_dflt$template_data$csd_counts)


##     ## Assumes genes in the order given
##     valid_gene_names <- c("WT", gene_names[1:n_genes])
##     ## Verify assumption. This is an assumption about input.
##     ##  gene_names must have been set correctly
##     ##  This is check A1_gnn

##     gene_names_in_freqs <- setdiff(unlist(strsplit(freqs$Genotype, ", ")), "WT")
##     if (length(gene_names_in_freqs)) {
##         if (n_genes >= (length(gene_names_in_freqs))) {
##             ## if (!all(gene_names_in_freqs %in% valid_gene_names)) {
##             ##     browser()
##             ## }
##             stopifnot(all(gene_names_in_freqs %in% valid_gene_names))
##         } else {
##             stopifnot(evamtools:::evam_string_sort(gene_names_in_freqs)[1:n_genes] %in% valid_gene_names)
##         }
##     }

##     ## Why would this be necessary?  To make sure size of data reduced when
##     ##    changing number of gens and we use genotype freqs.  But this works poorly
##     ##    with DAGs; it leads to incorrect behavior like having in dag only A and C,
##     ##    and asking for two genes. Only A would be shown.
##     ## So don't use this function with DAGs :-)

##     selected_rows <- vapply(freqs$Genotype,
##                             function(x) {
##                                 genes <- strsplit(x, ", ")[[1]]
##                                 return(all(genes %in% valid_gene_names))
##                             },
##                             logical(1))

##     freqs <- freqs[selected_rows, , drop = FALSE]
##     ## Remove 0 count rows
##     return(freqs[freqs$Counts > 0, , drop = FALSE])

##     }
