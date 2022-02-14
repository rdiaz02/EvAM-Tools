## Copyright 2021, 2022 Pablo Herrera Nieto

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Use plot matrix to plot the sampled genotypes
#' 
#' @param data data.frame object with cross sectional datas
plot_sampled_genots <- function(data) {
    d1 <- as.data.frame(sampledGenotypes(data))

    ## Reorder by number mutations and names, but WT always first
    which_wt <- which(d1[, 1] == "WT")
    if(length(which_wt)) {
        dwt <- d1[which_wt, ]
        dnwt <- d1[-which_wt, ]
    } else {
        dwt <- NULL
        dnwt <- d1
    }
    
    oo <- order(str_count(dnwt[, 1], ","), dnwt[, 1])
    dnwt <- dnwt[oo, ]
    d1 <- rbind(dwt, dnwt)
    
    m1 <- matrix(d1[, 2], ncol = 1)
    colnames(m1) <- "Freq."
    rownames(m1) <- d1[, 1]
    op <- par(mar = c(3, 5, 5, 3), las = 1)
    plot(m1, cex = 1.5, digits = 0, key = NULL,
         axis.col = list(side = 3), xlab = "", ylab = "",
         col = NULL,
         main  = "")
    par(op)
}


#' Order paths of a weigthed directd graph
#' 
#' Computes the relevance of the paths starting on WT based on edge weigth
#' 
#' @param graph igraph object with genotype transition. The graph is expect to be directed and weighted
#' @return List with all paths sorted in descreasing order of importance
rank_paths <- function(graph) {
    g <- graph
    all_paths <- list()
    ## Starting from WT --> get the most likely children
    ## recusively until reaching a leaf
    leaves <- V(g)[igraph::degree(g, mode = "out") == 0]
    all_paths <- igraph::all_simple_paths(g, from = "WT", to = leaves)
    cum_weigths <- vapply(all_paths, function(x){
        sum(igraph::E(g, path = x)$weight)
    }, numeric(1))

    return(all_paths[order(cum_weigths, decreasing = TRUE)])
}

#' Return the labels of most relevant paths starting from WT 
#' 
#' @param graph igraph object with genotype transitions
#' @param paths_from_graph List of path from WT to ending genotype
#' @param top_paths Int > 0. Include labels from vertex present in n top_paths
#' @param type String. Default genotype. Valid options are "genotypes" or "acquisition"
#' "genotype" option returns the genotype of the vertex.
#' "acquisition" option returns the genotype acquire along the path.
compute_vertex_labels <- function(graph, paths_from_graph, top_paths = NULL,
                                  type = "genotype") {
    if(is.null(top_paths)) top_paths <- length(paths_from_graph)
    else if(is.numeric(top_paths)){
        if(top_paths > length(paths_from_graph)) top_paths <- length(paths_from_graph)
        else if(top_paths > 0) top_paths <- top_paths
    } else top_paths <- length(paths_from_graph)

    paths_from_graph <- paths_from_graph[1:top_paths]
    edge_labels <- NULL
    
    # if(vertex_labels){
    nodes_in_top_paths <- unique(unlist(sapply(paths_from_graph,
        function(x) x$name
    )))
    vertex_labels <- vapply(V(graph)$name,
        function(x){
            if (x %in% nodes_in_top_paths){
                if (type == "genotype") return(x)
                else if (type == "acquisition") return(sprintf("+%s", tail(strsplit(x, ",")[[1]], n = 1)))
            } 
            else return("")
        },
        character(1))
    # } else vertex_labels <- NULL

    # if(edge_labels){
        ## TODO in each vertex say +NEW_GENE
    # }


    return(list(
        vertex_labels = vertex_labels,
        edge_labels = edge_labels
    ))
}

#' Computes the best layout for a genotype transtion matrix 
#' 
#' Genotypes are distributed in the x axis according the number 
#' of mutated genes they have
#' Distribution in the Y axis in done alphabetically.
#' 
#' @param graph igraph object with genotype transitions
#' 
#' @return dataframe the x and y position for each vertex in graph
cpm_layout <- function(graph){
    num_mutations <- igraph::distances(graph, v="WT", weights=NA)
    igraph::V(graph)$num_mutations <- num_mutations 
    lyt <- matrix(0, ncol = 2, nrow = length(V(graph)))
    lyt[, 2] <- num_mutations

    for (i in 0:max(num_mutations)) {
        level_idx <- which(num_mutations == i)
        gnt_names <- sort(igraph::V(graph)$name[level_idx], index.return = TRUE)
        spacing <- 6 / (length(level_idx) + 1)
        level_elements <- 1:length(level_idx) * spacing
        level_elements <-  rev(level_elements - max(level_elements))
        correction_rate <- max(level_elements) - min(level_elements) 
        level_elements <- level_elements + correction_rate/2
        lyt[, 1][level_idx[gnt_names$ix]] <- level_elements
    }

    ## Avoiding layout in one line
    if (all(lyt[, 1] == 0)) lyt[,1] <- rep(c(0,1,-1), ceiling(nrow(lyt)/3))[1:nrow(lyt)]
    return(lyt)
}

#' Plot hypercubic genotype transitions
#' 
#' @param trans_mat transitions matrix to plot: can contain row counts, probabilities... 
#' Entries will be normalized
#' @param observations Original cross sectional data used to compute the model. Optional.
#' @param freqs DataFrame with column $Genotype and $Freqs with their frequencies. Optional.
#' @param top_paths Int>0. Default NULL. Most transited paths to label
#' @param freq2label Int>0. Laberl genotypes with a frequency larger taht freqs2label
#' @param max_edge Int>0. Maximun width of edge. If NULL it will be infered from data.
#' @param min_edge Int>0. Minimum width of edge. If NULL it will be infered from data.
#' @param fixed_vertex_size Boolen. If TRUE, all nodes with all the same size and frequencies or observed data will not be used.
#' @examples 
#' \dontrun{
#' dB_c2 <- matrix(
#' c(
#'      rep(c(1, 0, 0, 0, 0), 300) #A
#'    , rep(c(0, 0, 1, 0, 0), 300) #C
#'    , rep(c(1, 1, 0, 0, 0), 200) #AB
#'    , rep(c(0, 0, 1, 1, 0), 200) #CD
#'    , rep(c(1, 1, 1, 0, 0), 100) #ABC
#'    , rep(c(1, 0, 1, 1, 0), 100) #ACD
#'    , rep(c(1, 1, 0, 0, 1), 100) #ABE
#'    , rep(c(0, 0, 1, 1, 1), 100) #CDE
#'    , rep(c(1, 1, 1, 0, 1), 100) #ABCE
#'    , rep(c(1, 0, 1, 1, 1), 100) #ACDE
#'    , rep(c(1, 1, 1, 1, 0), 50) # ABCD
#'    , rep(c(0, 0, 0, 0, 0), 10) # WT
#'  ), ncol = 5, byrow = TRUE
#' )
#' colnames(dB_c2) <- LETTERS[1:5]
#' out <- evam(dB_c2)
#' png("fluxes.png")
#' par(mfrow = c(1, 1))
#' plot_genot_fg(out$MHN_trans_mat, dB_c2)
#' dev.off()
#' }
plot_genot_fg <- function(trans_mat
    , observations = NULL
    , freqs = NULL
    , top_paths = NULL
    , freq2label = NULL
    , max_edge = NULL
    , min_edge = NULL
    , fixed_vertex_size = FALSE
    , label_type = "genotype"
    ) {
    if(is.null(trans_mat)){
        par(mar = rep(3, 4))
        plot(0, type = 'n', axes = FALSE, ann = FALSE)
        return()
    }

    unique_genes_names <- sort(unique(unlist(str_split(rownames(trans_mat)[-1], ", "))))
    rownames(trans_mat) <- colnames(trans_mat) <- str_replace_all(rownames(trans_mat), ", ", ",")
    graph <- igraph::graph_from_adjacency_matrix(trans_mat, weighted = TRUE)

    num_genes <- length(unique_genes_names)
    graph <- igraph::decompose(graph)[[1]] ## We do not want disconnected nodes
        
    if (!is.null(observations)){
        observations <- as.data.frame(sampledGenotypes(observations))
        observations$Abs_Freq <- observations$Freq / sum(observations$Freq)
        observations$Genotype <- str_replace_all(observations$Genotype, ", ", ",")
    }

    if (!is.null(freqs)){
        freqs$Abs_Freq <- freqs$Counts / sum(freqs$Counts)
        freqs$Genotype <- str_replace_all(freqs$Genotype, ", ", ",")
        rownames(freqs) <- freqs$Genotype
    }

    ## Layout
    lyt <- cpm_layout(graph)

    ## Labels
    sorted_paths <- rank_paths(graph)
    if(is.null(freq2label)){
        labels <- compute_vertex_labels(graph, sorted_paths, top_paths = top_paths, type=label_type)
    } else {
        labels <- vapply(igraph::V(graph)$name,
            function(x){
                if (freqs[x, ]$Abs_Freq >= freq2label) return(x)
                else return("")
            },
            character(1))
        labels <- list(vertex_labels = labels)
    }

    ## Vertex colors based on the presence/absence in the original data
    observed_color <- "#ff7b00"
    not_observed_color <- "#0892d0" 
    if(is.null(observations)){
        not_observed_color <- "#ff7b00"
    } 
    colors <- sapply(igraph::V(graph)$name, 
        function(gen){
            if (sum(match(observations$Genotype, gen, nomatch = 0)) == 1){
                return(observed_color)
            } 
            return(not_observed_color)
        })

    igraph::V(graph)$color <- colors
    igraph::V(graph)$frame.color <- colors

    ## Sizes based of frequency
    min_size <- 2
    max_size <- 40

    if(fixed_vertex_size){
        node_sizes <- vapply(igraph::V(graph)$name, 
        function(gen) min_size, numeric(1.0))
    } else if(!is.null(freqs)){
        node_sizes <- vapply(igraph::V(graph)$name, 
            function(gen){
                if (sum(match(freqs$Genotype, gen, nomatch = 0)) == 1)
                    return(freqs$Abs_Freq[which(freqs$Genotype == gen)])
                else 
                    return(min_size)
            }, numeric(1.0))
    } else if(!is.null(observations)){
        node_sizes <- vapply(igraph::V(graph)$name, 
            function(gen){
                if (sum(match(observations$Genotype, gen, nomatch = 0)) == 1) # Observed
                    return(observations$Abs_Freq[which(observations$Genotype == gen)])
                else # Not Observed
                    return(0)
            }, numeric(1.0))
    } else {
        node_sizes <- vapply(igraph::V(graph)$name, 
        function(gen) min_size, numeric(1.0))
    }

    node_sizes[node_sizes <= 0.01] <- 0.01
    if(length(unique(node_sizes))>1){ #Not unique values
        node_sizes <- (node_sizes - min(node_sizes))/(max(node_sizes) - min(node_sizes)) * (max_size - min_size) + min_size
    } else {
        node_sizes <- rep(15, length(node_sizes))
    }
    igraph::V(graph)$size <- node_sizes

    igraph::V(graph)$label.family <- "Helvetica"

    opx <- par(mar = c(3, 0, 3, 0))
    
    ## Weights proportional to the time they are transited
    w <- igraph::E(graph)$weight
    max_edge <- ifelse(is.null(max_edge), max(w), max_edge)
    min_edge <- ifelse(is.null(min_edge), min(w), min_edge)

    min_width <- 3
    max_width <- 13

    w2 <- NULL
    if(length(unique(w))>1){ #Not unique values
        w2 <- (w - min_edge)/(max_edge - min_edge) * (max_width - min_width) + min_width
    } else {
        w2 <- rep(min_width, length(node_sizes))
    }

    ## if(any(is.na(w2))) browser()
    transparent_w2 <-  w2/max(w2) * 0.9 + 0.1
    ## Actual plotting
    plot(graph
        , layout = lyt[, 2:1]
        , vertex.label.color = "black"
        , vertex.label.family = "Helvetica"
        , vertex.label.font = 2
        , vertex.label.cex = 1.2
        , vertex.label = labels$vertex_labels
        , font.best = 2
        , vertex.frame.width = 0
        # , edge.color = rgb(0.5, 0.5, 0.5, transparent_w2)
        , edge.color = rgb(0.5, 0.5, 0.5, 1) 
        ## Some error in my latop because, I cannot add alpha channel
        , edge.arrow.size = 0
        , xlab = "Number of features acquired"
        , edge.width = w2
    )

    margin <- -1.15
    lines(c(-1.2, 1.2), c(margin, margin), lwd = 2)
    node_depth <- vapply(igraph::V(graph)$name
                       , function(x)
                           distances <- igraph::distances(graph,
                                                          weights = NA,
                                                          algorithm = "unweighted",
                                                          to = x)["WT",],
                         3)
    max_node_depth <- max(node_depth[is.finite(node_depth)])
    axis(1
        , at = seq(-1, 1, length.out = max_node_depth + 1)
        , labels = 0:(max_node_depth)
        , lwd = 2
        , cex = 2
        , pos = margin)
    if(!(is.null(observations))){
        legend(-1, 1.6, c("Observed", "Not observed") 
            , box.lwd = 0, lty=c(NA, NA), lwd = c(NA, NA)
            , pch = c(21, 21), bty = "n"
            , col = c(observed_color, not_observed_color)
            , pt.bg = c(observed_color, not_observed_color)
            , pt.cex = c(2, 2), horiz = TRUE
            , x.intersp = c(0, 0)
            )
    }
    title(xlab = "Number of features acquired", line = 2)
}

#' Process data of CPMs to make it easier to plot
#' 
#' Extracts the outpus concerning a single CPM
#' 
#' @param data Complete CPM output
#' by calling sample_CPMs with the exact CPM output from data parameter
#' @param mod String for the CPM to process.
#' @param plot_type String for the plot_type to process.
#' @param sample_data Data form sampling, this must be generated 
#' @returns List with processed output of the CPM
process_data <- function(data, mod, plot_type, sample_data = NULL) {


    if (is.null(data)) {
        stop("data is NULL")
    }

    accepted_plot_types <- c("trans_mat",
                             "obs_genotype_transitions", "trans_rate_mat")
    
    if (!(plot_type %in% accepted_plot_types))
        stop("Incorrect plot_type. plot_type must be one of ",
             paste(accepted_plot_types, sep = ", ", collapse = ", "))
    

    all_data <- c(data, sample_data)
    
    dag_tree <- NULL
    tryCatch (expr = {
        dag_model <- get(paste(mod, "_model", sep = ""), data)
        dag_tree <- graph_from_data_frame(dag_model[, c(1, 2)])
    }, error = function(e){})

    if(!is.null(dag_tree)){
        model_info <- dag_tree
    }else if(!is.null(all_data[[paste0(mod, "_theta")]])){
        model_info <- all_data[[paste0(mod, "_theta")]]
    }

    return(list(
        model_info = model_info
        , data2plot = all_data[[paste0(mod, "_", plot_type)]]
        , parent_set = all_data[[paste0(mod, "_parent_set")]]
        , genotype_freqs = all_data[[paste0(mod, "_genotype_freqs")]]
        ))
}

dag_layout <- function(graph){ ## Avoiding lines
    lyt <- igraph::layout.reingold.tilford(graph)
    if(all(lyt[,1] == 0)) lyt[,1] <- rep(c(0,0.5,0,-0.5), ceiling(nrow(lyt)/3))[1:nrow(lyt)]
    return(lyt)
}

plot_model <- function(model_info, parent_set, mod){
    if (typeof(model_info) == "list") { ## Potting DAGs
        ## DAG relationships colors 
        standard_relationship <- "cornflowerblue"
        colors_relationships <- c(standard_relationship, standard_relationship,
                                "#E2D810", ##"#FBDE44FF",
                                "coral2")
        names(colors_relationships) <- c("Single", "AND", "OR", "XOR")
        g <- model_info
        if (!is.null(parent_set)) {
            for (i in igraph::E(g)) {
                igraph::E(g)[i]$color <-
                    colors_relationships[
                        parent_set[[igraph::head_of(g, igraph::E(g)[i])$name]]]
            }
        } else igraph::E(g)$color <- standard_relationship
        plot(g
            , layout = dag_layout
            , vertex.size = 50 
            , vertex.label.color = "black"
            , vertex.label.family = "Helvetica"
            , font.best = 2
            , vertex.frame.width = 0.5
            , vertex.color = "white"
            , vertex.frame.color = "black" 
            , vertex.label.cex = 1
            , edge.arrow.size = 0
            , edge.width = 5
            , main = mod)
        if(!is.null(parent_set)){
            legend("topleft", legend = names(colors_relationships),
                    col = colors_relationships, lty = 1, lwd = 5, bty = "n")
        }
    } else if(is.matrix(model_info)) { ##Plotting matrix
        op <- par(mar=c(3, 3, 5, 3), las = 1)
        plot(model_info, cex = 1.5, digits = 2, key = NULL
            , axis.col = list(side = 3)
            , xlab = "Effect of this (effector)"
            , ylab = " on this (affected)"
            , main = mod
            , mgp = c(2, 1, 0))
        par(op)
    } else {
        par(mar = rep(3, 4))
        plot(0, type = 'n', axes = FALSE, ann = FALSE)
    }
}

plot_CPMs <- function(cpm_output, samples = NULL, orientation = "horizontal", 
                        models = c("OT", "CBN", "OncoBN", "MCCBN", "MHN", "HESBCN"),
                        plot_type = "trans_mat", label_type="genotype",
                        fixed_vertex_size = FALSE,
                        top_paths = NULL) {

    if (!(plot_type %in% c("trans_mat", "trans_rate_mat", "obs_genotype_transitions"))){
        stop(sprintf("Plot type %s is not supported", plot_type))
    }

    if ((plot_type == "obs_genotype_transitions") & is.null(samples)){
        stop("obs_genotype_transitions needs you to pass the output ",
             "of a call to sample_CPMs")
    }
    
    ## List of available models
    available_models <- unique(models[
        vapply(models, function(mod) {
            attr_name <- ifelse(mod == "MHN", "theta", "model")
            return(any(!is.na(cpm_output[[sprintf("%s_%s", mod, attr_name)]])))
        }, logical(1))
    ])

    ## Shape of the plot
    l_models <- length(available_models)
    ## n_rows <- ifelse(old_plot_type == "matrix", 3, 2)
    n_rows <- 2
    
    if(orientation == "vertical") {
        op1 <- par(mfrow = c(l_models, n_rows))
    } else {
        op1 <- par(mfcol = c(n_rows, l_models))
    }
    par(mar = c(0.5, 1, 0.5, 0.5), mai = c(0.25, 0.25, 0.25, 0.25))

    ## Plotting CPMs
    ## For each CPM there are two plots
    ## 1. DAG or Matrix (for MHN)
    ## 2. Forward graph (HyperTRAPs style) of sampled transitions, transition matrix
    for(mod in available_models) {
        ## Processing data
        model_data2plot <- process_data(cpm_output, mod, plot_type, samples)
        
        ## Plotting model (DAG or MHN matrix)
        plot_model(model_data2plot$model_info, model_data2plot$parent_set, mod)

        ## Plotting forward graph
        if ((mod %in% c("OT")) &&
            (plot_type %in% c("trans_rate_mat", "obs_genotype_transitions"))) {
            par(mar = rep(3, 4))
            plot_sampled_genots(cpm_output$analyzed_data)
        } else { ## TODO remove all this
            # if (plot_type == "trans_mat") {## transition probabilities; was "probabilities"
            #     plot_genot_fg(model_data2plot$data2plot, 
            #                 cpm_output$analyzed_data,
            #                 top_paths = top_paths,
            #                 label_type = label_type,
            #                 fixed_vertex_size = fixed_vertex_size)
            # } else if (plot_type == "obs_genotype_transitions") { ## observed genotype trans. ; was "transitions"
            plot_genot_fg(model_data2plot$data2plot,
                        observations = cpm_output$analyzed_data,
                        model_data2plot$genotype_freqs,
                        top_paths = top_paths,
                        label_type = label_type,
                        fixed_vertex_size = fixed_vertex_size)
            # } else if (plot_type == "trans_rate_mat"){ ## transition rate matrix; was "trm"
            #     plot_genot_fg(model_data2plot$data2plot,
            #                 cpm_output$analyzed_data,
            #                 top_paths = top_paths,
            #                 label_type = label_type,
            #                 fixed_vertex_size = fixed_vertex_size)
            # }
        }
    }
    par(op1)
}

plot_genotypes_freqs <- function(data){
    if(nrow(data) == 0) return()
    par(las = 2, cex.main=1.6, cex.lab=1.5, cex.axis=1.2)
    barplot(data[, 2]
        , names = data$Genotype
        , ylab="Counts", main="Genotype Frequencies"
        , horiz = FALSE
        , panel.first=grid())
    grid(nx = NA, ny = NULL, col='gray', lwd = 2)
    ## TODO sort genotypes
}