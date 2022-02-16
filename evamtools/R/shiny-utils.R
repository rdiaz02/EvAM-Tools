

freqs2csd <- function(freqs, gene_names){
    csd <- NULL
    if(nrow(freqs) > 0){
        csd2 <- apply(
            freqs, 1
            , function(x){
                mut <- x[[1]]
                freq <- as.numeric(x[[2]])
                genot <- rep(0, length(gene_names))
                if(mut != "WT"){
                    mut <-  strsplit(mut, ", ")[[1]]
                    genot[which(gene_names %in% mut)] <- 1
                }
                csd <- matrix(rep(genot, freq), ncol= length(gene_names), byrow = TRUE)
                return(list(csd))
            })
        
        csd <- csd2[[1]][[1]]
        if(nrow(freqs) > 1){ 
            for (i in 2:nrow(freqs)){
                csd <- rbind(csd, csd2[[i]][[1]])
            }
        } 
        colnames(csd) <- gene_names
    }

    return(csd)
}

get_display_freqs <- function(freqs, n_genes, gene_names){
    if(is.null(freqs)) return(SHINY_DEFAULTS$template_data$csd_freqs)
    if(nrow(freqs) == 0) return(SHINY_DEFAULTS$template_data$csd_freqs)
    valid_gene_names <- c("WT", gene_names[1:n_genes])

    selected_rows <- sapply(freqs$Genotype, function(x){
        genes <- strsplit(x, ", ")[[1]]
        return(all(genes %in% valid_gene_names))
    })

    return(freqs[selected_rows, ])
}

get_csd <- function(complete_csd){
    if(is.null(complete_csd)) return(NULL)
    csd <- data.frame(OncoSimulR::sampledGenotypes(complete_csd))
    rownames(csd) <- csd$Genotype
    return(csd)
}


get_mhn_data <- function(n_genes, n_samples, gene_names, thetas = NULL){
    if(is.null(thetas)) thetas <- evamtools:::Random.Theta(n=n_genes)
    rownames(thetas) <- colnames(thetas) <- gene_names
    samples <- floor(evamtools:::Finite.Sample(evamtools:::Generate.pTh(thetas), n_samples)*n_samples)
    trm <- evamtools:::theta_to_trans_rate_3_SM(thetas,
                                    inner_transition = evamtools:::inner_transitionRate_3_1)
    state_names <- vapply(1:(ncol(trm)), function(x){
        x <- x - 1
        if(x == 0) state_name <- "WT"
        else state_name <- paste(gene_names[which(evamtools:::int2binary(x, n_genes) == 1)], collapse = ", ")
        return(state_name)
    }, character(1))
    rownames(trm) <- colnames(trm) <- state_names
    samples <- data.frame("Genotype" = state_names, "Freq" = samples)
    rownames(samples) <- samples$Genotype
    return(list(thetas = thetas, trm = trm, samples = samples))
}

create_tabular_data <- function(data, type){
    available_methods <- c("Source", "OT", "CBN", "MHN", "HESBCN")
    # , "DBN", "MCCBN")
    if(type %in% c("freqs")){
        all_counts <- data.frame(Genotype = data[["MHN_genotype_freqs"]]$Genotype)
        for(name in names(data)){
            if(grepl("_genotype_freqs", name)
            #  & !grepl("^OT", name) ##Now we have genotypes frequencies for OT
             ){
                method_name <- strsplit(name, "_")[[1]][[1]]
                all_counts[[method_name]] <- data[[name]]$Counts
            } 
        }
        
        order_by_counts <- sort(rowSums(all_counts[-1]), 
        decreasing = TRUE, index.return = TRUE)$ix
    
        all_counts[order_by_counts, ]
        return(all_counts[order_by_counts, ])

    } else if(type %in% c("trans_rate_mat", "genotype_transitions", "trans_mat", "td_trans_mat")){
        var2var <- c("trans_rate_mat", "genotype_transitions", "trans_mat", "td_trans_mat")
        names(var2var) <- c("trans_rate_mat", "genotype_transitions", "trans_mat", "td_trans_mat")
        
        var2use <- var2var[type]
        ## 1 Methods to compute
        methods2compute <- vapply(available_methods, function(x){
            if(!is.null(data[[sprintf("%s_%s", x, var2use)]])){
                return(!any(is.na(data[[sprintf("%s_%s", x, var2use)]])))
            }
            return(FALSE)
        }, logical(1))

        methods2compute <- names(methods2compute)[methods2compute]
        if(type == "genotype_transitions"){
            methods2compute <- setdiff(methods2compute, "OT")
        }
            
        ## 2 Fill the data frame
        all_genotypes <- matrix(0, nrow = 0, ncol = 2)
        all_methods <- matrix(0, nrow = 0, ncol = length(methods2compute))
        n_methods2compute <- length(methods2compute)
        base_vector <- rep(0, n_methods2compute)
        names(base_vector) <- methods2compute
        ## I use MHN because it has all the genotypes
        for(i in rownames(data[[sprintf("MHN_%s", var2use)]])){
            for(j in rownames(data[[sprintf("MHN_%s", var2use)]])){
                row_data <- sapply(methods2compute, function(x){
                    tryCatch({
                        return(data[[sprintf("%s_%s", x, var2use)]][i, j])
                    }, error = function(e){
                        return(0)
                    })
                })

                all_genotypes <- rbind(all_genotypes, c(i, j))
                all_methods <- rbind(all_methods, unname(row_data))
            }
        }
        
        selected_rows <- rowSums(abs(all_methods))>0
        all_methods <- round(all_methods[selected_rows, ], 2)
        all_genotypes <- all_genotypes[selected_rows, ]

        all_the_data <- data.frame(From = all_genotypes[, 1]
            , To = all_genotypes[, 2])
        for(i in 1:n_methods2compute){
            all_the_data[[methods2compute[i]]] <- all_methods[, i]
        }

        colnames(all_the_data) <- c("From", "To", methods2compute)

        order_by_counts <- sort(rowSums(all_methods), 
        decreasing = TRUE, index.return = TRUE)$ix
    
        return(all_the_data[order_by_counts, ])

    } else if(type %in% c("lambdas")){
        lambda_field <- c("Lambdas", "OT_edgeWeight", "rerun_lambda", "Lambdas", "lambda", "Thetas")
        names(lambda_field) <- c("Source", "OT", "CBN", "HESBCN", "MCCBN")

        gene_names <- sort(unique(data$OT_model$To))
        all_counts <- data.frame(Gene = gene_names)
        for(name in names(data)){
            if(grepl("_model", name)){
                method_name <- strsplit(name, "_")[[1]][[1]]
                if(!is.null(data[[name]]) & !is.na(data[[name]])){
                    tmp_data <- data[[name]][[lambda_field[method_name]]]
                    names(tmp_data) <- data[[name]]$To
                    all_counts[[method_name]] <- round(tmp_data[all_counts$Gene], 2)
                }
            }
        }

        return(all_counts)
    }

    order_by_counts <- sort(rowSums(all_counts[-1]), 
        decreasing = TRUE, index.return = TRUE)$ix
    
    all_counts[order_by_counts, ]
    return(to_return)
}


## Fills empty attributes of a csd data object
standarize_dataset <- function(data){
  if(is.null(data)) return(SHINY_DEFAULTS$template_data)

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

  new_data$gene_names <- c(tmp_names 
    , LETTERS[(length(tmp_names) + 1) : SHINY_DEFAULTS$max_genes ])[1:SHINY_DEFAULTS$max_genes]
  
  new_data$name <- data$name  

  if(is.null(data$lambdas)) {
    new_data$lambdas <- SHINY_DEFAULTS$template_data$lambdas
  }else{
    if(!is.numeric(data$lambdas)){
      stop("Lambdas should only contain numbers")
    }
    new_lambdas <- SHINY_DEFAULTS$template_data$lambdas
    new_lambdas[1:length(data$lambdas)] <- data$lambdas
    new_data$lambdas <- new_lambdas
  }
  names(new_data$lambdas) <- new_data$gene_names

  if(is.null(data$dag_parent_set)) {
    new_data$dag_parent_set <- SHINY_DEFAULTS$template_data$dag_parent_set
  }else {
    if(!all(data$dag_parent_set %in% c("Single", "AND", "OR", "XOR")))
      stop("Parent set must include only 'Single', 'AND', 'OR' or 'XOR'")
    new_parent_set <- SHINY_DEFAULTS$template_data$dag_parent_set
    new_parent_set[1:length(data$dag_parent_set)] <- data$dag_parent_set
    new_data$dag_parent_set <- new_parent_set
  }
  names(new_data$dag_parent_set ) <- new_data$gene_names

  if(is.null(data$dag)) {
    new_data$dag <- SHINY_DEFAULTS$template_data$dag
  } else {
    if(!all(unique(data$dag) %in% c(0, 1))){
      stop("Adjacency matrix should be binary: only 0 and 1")
    }
    to_keep <- which(colSums(data$dag)>0 | rowSums(data$dag)>0)
    tmp_dag <- data$dag[to_keep, to_keep]
    n_genes <- ncol(data$dag)
    new_dag <- SHINY_DEFAULTS$template_data$dag 
    new_dag[to_keep, to_keep] <- tmp_dag
    new_data$dag <- new_dag

    if(!is_dag(graph_from_adjacency_matrix(tmp_dag, mode="directed"))){
      stop("The graph is not a DAG")
    }

    ## Revising parent set
    ## Only genes with multiple parent can have something different from "Single" relationship
    new_data$dag_parent_set[colSums(new_dag)[-1] <= 1] <- "Single"
  }
  rownames(new_data$dag) <- colnames(new_data$dag) <- c("WT", new_data$gene_names)

  if(is.null(data$thetas)) {
    new_data$thetas <- SHINY_DEFAULTS$template_data$thetas
  } else {
    if(!is.numeric(data$thetas)){
      stop("Theta matrix should only contain numbers")
    }
    n_genes <- ncol(data$thetas)
    new_thetas <- SHINY_DEFAULTS$template_data$thetas 
    new_thetas[1:n_genes, 1:n_genes] <- data$thetas
    new_data$thetas <- new_thetas
  }
  rownames(new_data$thetas) <- colnames(new_data$thetas) <- new_data$gene_names

  if(is.null(data$data)) {
    new_data$data <- SHINY_DEFAULTS$template_data$data
  } else {
    if(!all(unique(as.vector(data$data)) %in% c(0, 1))){
      stop("Data should be binary: only 0 and 1")
    }
    new_data$data <- data$data
    colnames(new_data$data) <- new_data$gene_names[1:ncol(new_data$data)]
  }
  
  new_data$csd_freqs <- get_csd(new_data$data)
  return(new_data)
}

standarize_all_datasets <- function(datasets){
  all_new_data <- list()
  for(i in c("csd", "dag", "matrix")){
    tmp_data <- datasets[[i]] 
    for(j in names(tmp_data)) all_new_data[[i]][[j]] <-     
      standarize_dataset(tmp_data[[j]])
  }
  return(all_new_data)
}


# plot_model <- function(cpm_output, mod){
#     ## DAG relationships colors 
#     standard_relationship <- "gray73"
#     colors_relationships <- c(standard_relationship, standard_relationship, "cornflowerblue", "coral2")
#     names(colors_relationships) <- c("Single", "AND", "OR", "XOR")
    
#     model_data2plot <- evamtools:::process_data(cpm_output, mod)
#     ## Plotting data
#     if(!is.null(model_data2plot$dag_tree)) {

#         if(!is.null(model_data2plot$parent_set)){
#             for(i in names(model_data2plot$parent_set)){
#                 igraph::E(model_data2plot$dag_tree)[.to(i)]$color <- colors_relationships[model_data2plot$parent_set[[i]]]
#             }
#         } else igraph::E(model_data2plot$dag_tree)$color <- standard_relationship
#         plot(model_data2plot$dag_tree
#             , layout = igraph::layout.reingold.tilford
#             , vertex.size = 55 
#             , vertex.label.color = "black"
#             , vertex.label.family = "Helvetica"
#             , vertex.label.cex = 1.5
#             , font.best = 2
#             , vertex.frame.width = 0.5
#             , vertex.color = "white"
#             , vertex.frame.color = "black" 
#             , vertex.label.cex = 1
#             , edge.arrow.size = 0
#             , edge.width = 5
#             )
#         # par(mar=c(0, 0, 0, 0))
#         if(!is.null(model_data2plot$parent_set)){
#             legend("topleft", legend = names(colors_relationships),
#                 col = colors_relationships, lty = 1, lwd = 5, bty = "n")
#         }
#     }else if(!is.null(model_data2plot$theta)) {
#         op <- par(mar=c(3, 3, 5, 3), las = 1)
#         plot(model_data2plot$theta, cex = 1.4, digits = 2, key = NULL
#             , axis.col = list(side = 3)
#             , xlab = "Effect of this (effector)"
#             , ylab = " on this (affected)"
#             , main = ""
#             , mgp = c(2, 1, 0))
#         par(op)
#     }
#     title(mod, cex.main = 1.8)
# }




# plot_dag <- function(dag, parent_set){
#     if (is.null(dag)) return()
#     standard_relationship <- "gray73"
#     colors_relationships <- c(standard_relationship, standard_relationship, "cornflowerblue", "coral2")
#     names(colors_relationships) <- c("Single", "AND", "OR", "XOR")
#     dag <- igraph::graph_from_adjacency_matrix(dag, mode = "directed")
#     dag <- igraph::decompose(dag)[[1]]
#     ## Plotting data
#     if(!is.null(parent_set)){
#         for(i in igraph::E(dag)){
#             igraph::E(dag)[i]$color <- colors_relationships[parent_set[[igraph::head_of(dag, igraph::E(dag)[i])$name]]]
#         }
#     } else igraph::E(dag)$color <- standard_relationship
        
#     plot(dag
#         , layout = igraph::layout.reingold.tilford
#         , vertex.size = 30 
#         , vertex.label.color = "black"
#         , vertex.label.family = "Helvetica"
#         , vertex.label.cex = 1.5
#         , font.best = 2
#         , vertex.frame.width = 0.5
#         , vertex.color = "white"
#         , vertex.frame.color = "black" 
#         , vertex.label.cex = 1
#         , edge.arrow.size = 0
#         , edge.width = 5
#         )
#     if(!is.null(parent_set)){
#         legend("topleft", legend = names(colors_relationships),
#             col = colors_relationships, lty = 1, lwd = 5, bty = "n")
#     }
#     title("Direct acyclic graph", cex.main = 1.8)
# }


# check_if_csd <- function(data){
#     tmp_names <- c("data", "dag", "dag_parent_set", "gene_names", "name", "type", "thetas")
#     types <- c("csd", "dag", "matrix")
#     if(is.null(data)) return(FALSE)
#     if(all(tmp_names %in% names(data))){
#         if((data$type %in% types 
#         & all(unique(c(data$data)) %in% c(0, 1)))){
#             return(TRUE)
#         } else {return(FALSE)}
#     } else { ## We assume a data.frame
#         tmp_data <- unlist(data)
#         names(tmp_data) <- NULL
#         return(all(unique(tmp_data) %in% c(0, 1)))
#     }
# }


# available_cpms <- function(data){
#     data$csd_data <- NULL
#     cpm_names <- unique(sapply(names(data), function(x) str_split(x, "_")[[1]][[1]]))
#     return(cpm_names)
# }