## Reasons for not using this function:

## That it takes a gene_names is questionable:
## gene_names should be extracted from freqs.
## Doing it this way suggests not self-contained and it is unclear
## why not self-contained.

## Could errors be introduced if genes not ordered or in numbers that differ from
## expected? I do not know and there are no tests. For example: what would happen
## if freqs had columns A, B, but gene_names was A, D, E? Or if gene_names only
## had A?  One could add checks for that, but since the role of this is to
## produce tabular data, it seems sensible to assume that all the info is in
## freqs. (I tried adding a check but I think it is just simpler to used
## genotypeCounts_to_data).

## There is no checking for freqs having specifically named
## columns: it is assumed the first is genotype and the second is
## Counts/Frequency.

## Finally: using genotypeCounts_to_data we only output columns
## that have at least one mutation, never columns filled with only 0s.
freqs2csd <- function(freqs, gene_names) {
    genes_in_freqs <- setdiff(unique(unlist(strsplit(freqs[, 1], ", "))),
                              "WT")
    stopifnot(isTRUE(all(genes_in_freqs %in% gene_names)))
    
    csd <- NULL
    if (nrow(freqs) > 0) {
        csd2 <- apply(
            freqs, 1,
            function(x) {
                mut <- x[[1]]
                freq <- as.numeric(x[[2]])
                genot <- rep(0, length(gene_names))
                if (mut != "WT") {
                    mut <-  strsplit(mut, ", ")[[1]]
                    genot[which(gene_names %in% mut)] <- 1
                }
                csd <- matrix(rep(genot, freq), ncol = length(gene_names),
                              byrow = TRUE)
                return(list(csd))
            })
        
        csd <- csd2[[1]][[1]]
        if (nrow(freqs) > 1) {
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
    if(is.null(complete_csd)) return(SHINY_DEFAULTS$template_data$csd_freqs)
    csd <- data.frame(evamtools:::data_to_counts(complete_csd, out = "data.frame", omit_0 = TRUE))
    rownames(csd) <- csd$Genotype
    return(csd)
}

modify_dag <- function(dag, from_node, to_node, operation, parent_set){
  value_from_operation <- c(1, 0)
  names(value_from_operation) <- c("add", "remove")

  if(operation == "clear") {
    return(list(dag=SHINY_DEFAULTS$template_data$dag
      , dag_parent_set=SHINY_DEFAULTS$template_data$dag_parent_set))
  }

  if(is.null(from_node) | is.null(to_node) | is.null(dag)){
    stop("From and To options and DAG have to be defined")
  } 


  if(!all(names(parent_set) == colnames(dag)[-1])){
    stop("Dag and parent set should have information of the same genes")
  }

  from_node <- ifelse(from_node == "Root", "WT", from_node)
  if(!(from_node %in% colnames(dag)) | !(to_node %in% colnames(dag))){
    stop("Both From and To options have to be valid gene names")
  } else if(!all(unique(as.vector(dag)) %in% c(0, 1))){
    stop("The adjacency matrix should only contain 0 and 1")
  } else if (from_node == to_node){
    stop("Both From and To options must be different")
  } else if(dag[to_node, from_node] == 1 && operation == "add"){
    stop("Relationships cannot be bidirectional")
  } else if (dag[from_node, to_node] == 1 && operation == "add"){
    stop("That edge is already present")
  } else if(dag["WT", to_node] == 1 && operation == "add"){
    stop("A direct children of Root cannot have multiple parents")
  } 
  
  tmp_dag <- dag
  tmp_dag[from_node, to_node] = value_from_operation[operation]

   ##Taking subcomponent starting from WT
  g1 <- igraph::graph_from_adjacency_matrix(tmp_dag, mode = "directed")
  # browser()
  for(tmp_g in igraph::decompose(g1)){
    if("WT" %in% V(tmp_g)$name){
      tmp_vertices <- igraph::as_data_frame(tmp_g)
    }
  }
  
  tmp_dag2 <- SHINY_DEFAULTS$template_data$dag
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
      ## We add link to WT
      tmp_dag2["WT", i] <- 1
    }
  }

  if(all(dag == tmp_dag2)){
    stop("This operation had no effect.")
  }
  g2 <- igraph::graph_from_adjacency_matrix(tmp_dag2, mode = "directed")
  if(!igraph::is_dag(g2)){
    stop("This relationship breaks the DAG. Revise it.")
  }
  ## Recompute parent set
  number_of_parents <- colSums(tmp_dag2)[-1]
  tmp_parent_set <- parent_set
  tmp_parent_set[number_of_parents > 1] <- "AND" ##Default is AND
  tmp_parent_set[number_of_parents <= 1] <- "Single" ##Default is AND
  tmp_parent_set[number_of_parents > 1 & !(parent_set %in% c("AND", "OR", "XOR"))] <- "AND" ##Recapture previous relationships
  return(list(dag = tmp_dag2, parent_set = tmp_parent_set))
}

modify_lambdas_and_parent_set_from_table <- function(dag_data, info, lambdas, dag, parent_set){
  if(!all(names(lambdas) == names(parent_set))){
    stop("Lambdas and parent set should have information about the same genes")
  }
  if(!all(names(lambdas) == colnames(dag)[-1])){
    stop("DAG should have information about the same genes as parent set")
  }

  ## Different lambdas
  new_lambdas <- as.numeric(info[info["col"] == 3,"value"])
  if(any(is.na(new_lambdas))){
    stop("There are missing lambdas")
  }
  old_lambdas <- dag_data$Lambdas
  all_genes <- dag_data$To
  changed_genes <- all_genes[new_lambdas != old_lambdas
    & new_lambdas > 0]
  info_from <- info[info["col"] == 0,"value"]
  info_to <- info[info["col"] == 1,"value"]
  if(!(all(c(info_from, info_to, dag_data$To, dag_data$From) %in% c("Root", names(lambdas))))){
    stop("There are unkown genes")
  }
  changed_lambdas <- new_lambdas[new_lambdas != old_lambdas
    & new_lambdas > 0]

  tmp_lambdas <- lambdas
  tmp_lambdas[changed_genes] <- changed_lambdas

  ##Relationships
  number_of_parents <- colSums(dag)[-1]
  number_of_parents <- number_of_parents[number_of_parents > 0]
  tmp_parent_set <- parent_set
  new_relationships <- info[info["col"] == 2,"value"]
  names(new_relationships) <- info[info["col"] == 1,"value"]
  
  #Cleaning unwanted value
  genes_in_relationships <- names(new_relationships)
  freq_genes_in_relationships <- table(genes_in_relationships)

  old_relationships <- stats::setNames(dag_data$Relation, dag_data$To)
  new_relationships[!(new_relationships %in% c("Single", "AND", "OR", "XOR"))] <- "AND"
  new_relationships[names(number_of_parents[number_of_parents<2])] <- "Single"
  new_relationships[setdiff(unique(genes_in_relationships), 
    names(freq_genes_in_relationships[freq_genes_in_relationships > 1]))] <- "Single"
  
  multiple_parents <- names(number_of_parents[number_of_parents > 1])

  for(i in multiple_parents){
    changed_relationships <- setdiff(new_relationships[genes_in_relationships == i], c(parent_set[i]))
    if (length(changed_relationships) > 0){
      tmp_parent_set[i] <- changed_relationships[[1]]
    }
  }
  # for(idx in c(1:length(new_relationships))){
  #   tmp_parent_set[names(new_relationships[idx])] <- new_relationships[idx]
  # }
  # if(is.null(new_relationships)  && !is.null(data$dag_parent_set)){
  #   tmp_parent_set <- data$dag_parent_set
  # } else {
  #   tmp_parent_set[!(tmp_parent_set %in% c("Single", "AND", "OR", "XOR"))] <- "AND"
  #   tmp_parent_set[number_of_parents <= 1] <- "Single"
  # }
  return(list(lambdas = tmp_lambdas, parent_set = tmp_parent_set))
}

get_mhn_data <- function(thetas, N = 10000){
  if (any(is.null(thetas))) stop("Thetas should be defined")
  n_genes <- ncol(thetas)
  gene_names <- colnames(thetas)

  mhn_trm <- MHN_from_thetas(thetas)$MHN_trans_rate_mat
  mhn_probs <- probs_from_trm(mhn_trm)
  tmp_samples_as_df <- genot_probs_2_pD_ordered_sample(x = mhn_probs,
                                                       ngenes = n_genes,
                                                       gene_names = gene_names,
                                                       N = N,
                                                       out = "data.frame"
                                                       )

  ## return(list(csdfreqs = tmp_samples_as_df,
  ##             data = freqs2csd(tmp_samples_as_df, gene_names)))

  return(list(csdfreqs = tmp_samples_as_df,
              data = genotypeCounts_to_data(tmp_samples_as_df, e = 0)))


}

get_dag_data <- function(data, parent_set, N = 10000){
  if (any(is.null(data))) stop("Data should be defined")
  gene_names <- names(parent_set)
  n_genes <- length(parent_set)

  dag_trm <- HESBCN_model_2_output(data, parent_set)$HESBCN_trans_rate_mat
  dag_probs <- probs_from_trm(dag_trm)
  tmp_samples_as_df <- genot_probs_2_pD_ordered_sample(x = dag_probs,
                                                       ngenes = n_genes,
                                                       gene_names = gene_names,
                                                       N = N,
                                                       out = "data.frame"
                                                       )
  ## ## FIXME: aqui. To rm later
  ## tmp_genotypes_sampled <- sample_to_pD_order(
  ##               sample(names(dag_probs), size = N,
  ##                      prob = dag_probs, replace = TRUE),
  ##               ngenes = n_genes, gene_names = gene_names)
  ## tmp_samples <- setNames(tmp_genotypes_sampled,
  ##           generate_pD_sorted_genotypes(n_genes, gene_names))
  ## tmp_samples <- tmp_samples[tmp_samples > 0]
  ## tmp_samples_as_df <- data.frame(Genotype=names(tmp_samples), Counts=tmp_samples)

  ## Note: no, the above call to genot_probs_2_pD ... does not remove
  ## lines with 0 counts. But we do not want to do that as that is a separate
  ## operation that has to with, maybe, generating the csd.

  ## return(list(csd_freqs = tmp_samples_as_df,
  ##             data = freqs2csd(tmp_samples_as_df, gene_names)))
  return(list(csd_freqs = tmp_samples_as_df,
              data = genotypeCounts_to_data(tmp_samples_as_df, e = 0)))

  
}

create_tabular_data <- function(data) {
    available_methods <- c("OT", "OncoBN", "CBN", "MHN", "HESBCN", "MCCBN")
    attr_to_make_tabular <- c("trans_mat", "trans_rate_mat",
                              "obs_genotype_transitions",
                              "predicted_genotype_freqs",
                              "sampled_genotype_freqs")
    tabular_data <- list()
    for (attr in attr_to_make_tabular){
        if (attr %in% c("sampled_genotype_freqs",
                        "predicted_genotype_freqs")) {
        all_counts <- data.frame(Genotype = c(NULL))

        for (method in available_methods){
          tmp_data <- data[[paste0(method, "_", attr)]]
          if (!any(is.null(tmp_data)) && !any(is.na(tmp_data))) {
            if (any(is.null(all_counts$Genotype))){
              all_counts <- data.frame(Genotype = names(tmp_data)) #They include all genotypes
              rownames(all_counts) <- all_counts$Genotype
            }
            all_counts[, method] <- round(tmp_data, 3)
          }
        }
        ## order_by_counts <- order(rowSums(all_counts[-1]),
        ##                          decreasing = TRUE)
        ## tabular_data[[attr]] <- all_counts[order_by_counts, ]
        all_counts <-
            data.frame(Index = standard_rank_genots_1(all_counts$Genotype),
                       all_counts)
        tabular_data[[attr]] <- all_counts[order(all_counts$Index), ]

      } else if (attr %in% c("trans_rate_mat",
                             "obs_genotype_transitions", "trans_mat")) {
          df <- data.frame(From = character(),
                 To = character(), 
                 OT = numeric(), 
                 OncoBN = numeric(), 
                 CBN = numeric(), 
                 HESBCN = numeric(), 
                 MCCBN = numeric(), 
                 MHN = numeric(), 
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
                  genes <- strsplit(transition, "->")[[1]]
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

          df <- data.frame(Index = standard_rank_genots_2(df$To, df$From),
                           df)
          ## Sorting 
          ## order_by_counts <- order(rowSums(df[,3:ncol(df)]),
          ##                          decreasing = TRUE)

          ## tabular_data[[attr]] <- df[order_by_counts, ]
          tabular_data[[attr]] <- df[order(df$Index), ]
      }
      # else if(type %in% c("lambdas")){
      #     lambda_field <- c("Lambdas", "OT_edgeWeight", "rerun_lambda", "Lambdas", "lambda", "Thetas")
      #     names(lambda_field) <- c("Source", "OT", "CBN", "HESBCN", "MCCBN")

      #     gene_names <- sort(unique(data$OT_model$To))
      #     all_counts <- data.frame(Gene = gene_names)
      #     for(name in names(data)){
      #         if(grepl("_model", name)){
      #             method_name <- strsplit(name, "_")[[1]][[1]]
      #             if(!is.null(data[[name]]) & !is.na(data[[name]])){
      #                 tmp_data <- data[[name]][[lambda_field[method_name]]]
      #                 names(tmp_data) <- data[[name]]$To
      #                 all_counts[[method_name]] <- round(tmp_data[all_counts$Gene], 2)
      #             }
      #         }
      #     }

      #     return(all_counts)
      # }
    }
    return(tabular_data)
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
      if(!all(unique(unlist(data$data)) %in% c(0, 1))){
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
    
#     model_data2plot <- process_data(cpm_output, mod)
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
