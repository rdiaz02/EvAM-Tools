
## use plot matrix to plot the sampled genotypes
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
         main  = "")
    par(op)
}


## output from all_methods_2_trans_mat, data for methods -> figures
##   OT: tree of restrictions, transition probs between genotypes
##        table of genotype freqs. (for lack of better idea)            
##   CBN: DAG of restrictions, transition probs between genotypes,
##                             transition probs genots TD
##   MCCBN: same as CBN 
##   MHN: theta matrix, transition probs between genotypes,
##                             transition probs genots TD
cpm_layout <- function(graph){
    # V(graph)$num_mutations <<- 
    num_mutations <- sapply(V(graph)$name, function(x){
        ifelse(x == "WT", 0, nchar(x))
    })
    V(graph)$num_mutations <- num_mutations 
    
    lyt <- matrix(0, ncol = 2, nrow = length(V(graph)))
    lyt[, 2] <-  num_mutations

    for (i in 0:max(num_mutations)) {
        level_idx <- which(num_mutations == i)
        gnt_names <- sort(V(graph)$name[level_idx], index.return = TRUE)
        spacing <- 6 / (length(level_idx) + 1)
        level_elements <- 1:length(level_idx) * spacing
        level_elements <-  rev(level_elements - max(level_elements))
        correction_rate <- max(level_elements) - min(level_elements) 
        level_elements <- level_elements + correction_rate/2
        lyt[, 1][level_idx[gnt_names$ix]] <- level_elements
        # browser()
    }
    return(lyt)
}

#' Plot transitions between genotypes
#' 
#' @param trans_mat transitions matrix to plot: can contain row counts, probabilities... 
#' Entries will be normalized
#' @param observations Original cross sectional data used to compute the model. Optional.
#' @param freqs DataFrame with column $Genotype and $Freqs with their frequencies. Optional.
#' @examples
#' png("fluxes.png")
#' par(mfrow = c(1, 3))
#' plot_genot_fg(out$MHN_trans_mat, db2, sorted_observations)
#' title("Transition Matrix", line = -3)
#' plot_genot_fg(edge_transitions, db2, sorted_observations)
#' title("Fluxes", line = -3)
#' dev.off()
plot_genot_fg <- function(trans_mat
    , observations = NULL
    , freqs = NULL
    , simplify = TRUE
    , top_edge_per_node = 5){
    # trans_mat <- as.matrix(trans_mat)
    rownames(trans_mat) <- str_replace_all(rownames(trans_mat), ", ", "")
    colnames(trans_mat) <- str_replace_all(colnames(trans_mat), ", ", "")

    graph <- graph_from_adjacency_matrix(trans_mat, weighted = TRUE)
    # print(sprintf("->Number of nodes %s Number of edges %s", length(V(graph)), length(E(graph))))

    unique_genes_names <- sort(unique(str_split(paste(rownames(trans_mat)[-1], collapse=""), "")[[1]]))

    num_genes <- length(unique_genes_names)
    graph <- graph_from_adjacency_matrix(trans_mat, weighted = TRUE)
    # browser()
    ## TODO take into account that a genoytpe can be observed
    if(simplify){
        graph <- decompose(graph)[[1]]
        min_values <- sort(trans_mat[trans_mat > 0], decreasing = TRUE)
        thr <- num_genes * top_edge_per_node
        ifelse(length(min_values > thr)
            , min_value <- min_values[thr]
            , min_value <- min_values[-1]
        )
        trans_mat[trans_mat < min_value] = 0

        graph <- graph_from_adjacency_matrix(trans_mat, weighted = TRUE)
        
        subgraphs <- decompose(graph)
        graph <- subgraphs[[1]]
        for(i in subgraphs[-1]){
            if(length(V(i)) > 5){
                EL  = get.edgelist(graph)
                EL1 = get.edgelist(i)
                ELU = rbind(EL, EL1)
                ELU = ELU[!duplicated(ELU),]
                w <- c(get.edge.attribute(graph, "weight")
                    , get.edge.attribute(i, "weight"))
                graph <- graph_from_edgelist(ELU)
                graph <- set_edge_attr(graph, "weight", value = w)
            } 
        }
    } 
        
    if (!is.null(observations)){
        observations <- as.data.frame(sampledGenotypes(observations))
        observations$Abs_Freq <- observations$Freq / sum(observations$Freq)
        observations$Genotype <- str_replace_all(observations$Genotype, ", ", "")
    }

    if (!is.null(freqs)){
        freqs$Abs_Freq <- freqs$Counts / sum(freqs$Counts)
        freqs$Genotype <- str_replace_all(freqs$Genotype, ", ", "")
    }

    lyt <- cpm_layout(graph)
    # browser()

    observed_color <- "#ff7b00"
    not_observed_color <- "#0892d0" 
    if(is.null(observations)){
        not_observed_color <- "#ff7b00"
    } 
    colors <- sapply(V(graph)$name, 
        function(gen){
            if (sum(match(observations$Genotype, gen, nomatch = 0)) == 1){
                return(observed_color)
            } 
            return(not_observed_color)
        })
    V(graph)$color <- colors
    V(graph)$frame.color <- colors

    sizes <- vapply(V(graph)$name, 
        function(gen){
            if (sum(match(freqs$Genotype, gen, nomatch = 0)) == 1){
                return(freqs$Abs_Freq[which(freqs$Genotype == gen)] * 150)
            } else if ((sum(match(observations$Genotype, gen, nomatch = 0)) == 1)){
                return(observations$Abs_Freq[which(observations$Genotype == gen)] * 150)
            } else {
                return(7.5)
            }
        }, numeric(1.0))

    if(all(sizes == 7.5)) sizes <- rep(15, length(sizes))
    V(graph)$size <- sizes

    V(graph)$label.family <- "Helvetica"

    opx <- par(mar = c(1, 1, 1, 1))

    w <- E(graph)$weight
    w <- w / max(w) * 10
    if(all(w == 10)) w <- rep(1, length(w))
    plot(graph
        , layout = lyt[, 2:1]
        , vertex.label.color = "black"
        , vertex.label.family = "Helvetica"
        , font.best = 2
        , vertex.label.cex = 1
        , vertex.frame.width = 0
        , edge.color = rgb(0.5, 0.5, 0.5, 1)
        # , edge.color = rgb(0.5, 0.5, 0.5, E(graph)$weight/max(E(graph)$weight))
        , edge.arrow.size = 0.5
        , xlab = "Number of features acquired"
        , edge.width = w
    )

    margin <- -1.15
    lines(c(-1.2, 1.2), c(margin, margin), lwd = 2)
    node_depth <- sapply(V(graph)$name
        , function(x) distances <- distances(graph, algorithm = "unweighted", to = x)["WT",]
        )
    max_node_depth <- max(node_depth[is.finite(node_depth)])
    axis(1
        , at = seq(-1, 1, length.out = max_node_depth + 1)
        , labels = 0:(max_node_depth)
        , lwd = 2
        , cex = 2
        , pos = margin)
    if(!(is.null(observations))){
        legend("bottom", c("Observed", "Not observed") 
            , box.lwd = 0, lty=c(NA, NA), lwd = c(NA, NA)
            , pch = c(21, 21)
            , col = c(observed_color, not_observed_color)
            , pt.bg = c(observed_color, not_observed_color)
            , pt.cex = c(2, 2), horiz = TRUE
            , x.intersp = c(0, 0)
            )
    }
    par(opx)
    # title(xlab = "Number of features acquired", line = -3)
}

#' Plot results from CPMs
#' 
#' By default it create a top row with the DAG of the CPM 
#' or de transtionRateMatrix for MHN
#' The bottom row has a custom plot for transition between genotypes

#' @param x output from the cpm
#' @param data Original cross sectional data used to compute the model. Optional.
#' @param models Output of the CPMs to plot. Current support is for OT, CBN, DBN, MCCBN and MHN Optional.
#' @param orientation String. If it not "vertical" will be displayed with an horizontal layout. Optional.
#' @param plot_type String. You can choose between 4 options. Optional.
#' genotypes: DAG with genotypes transitions
#' matrix: respresents the transtion matrix as a heatmap
#' transitions: shows the transitions count between genotypes in HyperTraPS style
#'              running simulations is needed before this
#' trans_mat: HyperTraps-like representation of the transition matrix
#' 
#' @examples
#' out <- all_methods_2_trans_mat(dB_c1, do_MCCBN = TRUE)
#' png("trans_at.png", width = 1000, height = 600, units = "px")
#' plot_DAG_fg(out, dB_c1, plot_type = "trans_mat")
#' dev.off()
#' out2 <- run_all_simulations(out, 100, n_genes = 5)
#' png("graph.png", width = 1000, height = 600, units = "px")
#' plot_DAG_fg(out, dB_c1, plot_type = "transitions")
#' dev.off()
plot_DAG_fg <- function(x, data, orientation = "horizontal", 
                        models = c("OT", "CBN", "DBN", "MCCBN", "MHN", "HESBCN"),
                        plot_type = "trans_mat",
                        prune_edges = TRUE) {
    
    if (!(plot_type %in% c("matrix", "transitions", "trans_mat", "genotypes"))){
        stop(sprintf("Plot type %s is not supported", plot_type))
    }
    # plot_fg <- function(fg) {
    #     ## Ideas from: https://stackoverflow.com/a/48368540
    #     lyt <- layout.reingold.tilford(fg)
    #     opx <- par(mar=c(2, 0.5, 2, 0.5))
    #     plot(fg,
    #          ## If I plot them sideways, labels in self-transitions
    #          ## overlap. FIXME. This sucks, I want them sideways!
    #          ## layout = -lyt[, 2:1],
    #          layout = lyt,
    #          edge.label = round(E(fg)$weight, 2),
    #          vertex.color = "SkyBlue2",
    #          edge.label.color = "black")
    #     par(opx)
    # }

    process_data <- function(mod) {
        dag_tree <- NULL
        
        dag_tree <- NULL
        tryCatch (expr = {
            dag_model <- get(paste(mod, "_model", sep = ""), x)
            dag_tree <- graph_from_data_frame(dag_model[, c(1, 2)])
        }, error = function(e){})

        dag_trans_mat <- get(paste(mod, "_trans_mat", sep = ""), x)
        fg <- graph_from_adjacency_matrix(dag_trans_mat, weighted = TRUE)

        if(prune_edges) {
            dag_trans_mat[dag_trans_mat < 0.01] <- 0
        }

        if(plot_type == "matrix") {
            dag_trans_mat <- as.matrix(dag_trans_mat)
            dag_trans_mat <- dag_trans_mat[rowSums(dag_trans_mat) > 0, colSums(dag_trans_mat) > 0]
        }

        td_trans_mat <- NULL
        td_fg <- NULL
        tryCatch(expr = {
            td_trans_mat <- get(paste(mod, "_td_trans_mat", sep = ""), x)
            td_fg <- graph_from_adjacency_matrix(td_trans_mat, weighted = TRUE)
            if(prune_edges) {
                td_trans_mat[td_trans_mat < 0.01] <- 0
            }
            if (plot_type == "matrix"){
                td_trans_mat <- as.matrix(td_trans_mat)
                td_trans_mat <- td_trans_mat[rowSums(td_trans_mat) > 0, colSums(td_trans_mat) > 0]
            }
        }, error = function(e){ })

        theta <- NULL
        tryCatch(expr ={
            theta <- get(paste(mod, "_theta", sep=""), x)
        }, error = function(e) { })

        return(list(dag_tree = dag_tree
            , dag_trans_mat = dag_trans_mat
            , fg = fg
            , td_trans_mat = td_trans_mat
            , td_fg = td_fg
            , theta = theta
            , parent_set = x[[sprintf("%s_parent_set", mod)]]
            , transitions = x[[sprintf("%s_genotype_transitions", mod)]]
            ))
    }

    ## List of available models
    available_models <- models[
        vapply(models, function(mod) {
            attr_name <- ifelse(mod == "MHN", "theta", "model")
            return(any(!is.na(x[[sprintf("%s_%s", mod, attr_name)]])))
        }, logical(1))
    ]
    print(available_models)

    ## DAG relationships colors 
    standard_relationship <- "gray73"
    colors_relationships <- c(standard_relationship, "coral2", "cornflowerblue", "darkolivegreen3")
    names(colors_relationships) <- c("Single", "AND", "OR", "XOR")
    
    ## Shape of the plot
    l_models <- length(available_models)
    n_rows <- ifelse(plot_type == "matrix", 3, 2)

    if(orientation == "vertical") {
        op1 <- par(mfrow = c(l_models, n_rows))
    } else {
        op1 <- par(mfcol = c(n_rows, l_models))
    }
    par(mar = c(0.5, 1, 0.5, 0.5), mai = c(0.25, 0.25, 0.25, 0.25))

    ## Plotting models
    for(mod in available_models) {
        ## Processing data
        model_data2plot <- process_data(mod)
        ## Plotting data
        if(!is.null(model_data2plot$dag_tree)) {
            if(!is.null(model_data2plot$parent_set)){
                for(i in names(model_data2plot$parent_set)){
                    E(model_data2plot$dag_tree)[.to(i)]$color <- colors_relationships[model_data2plot$parent_set[[i]]]
                }
            } else E(model_data2plot$dag_tree)$color <- standard_relationship
            plot(model_data2plot$dag_tree
                , layout = layout.reingold.tilford
                , vertex.size = 30 
                , vertex.label.color = "black"
                , vertex.label.family = "Helvetica"
                , font.best = 2
                , vertex.frame.width = 0.5
                , vertex.color = "white"
                , vertex.frame.color = "black" 
                , vertex.label.cex = 1
                , edge.arrow.size = 0.45
                , edge.width = 5
                , main = mod)
            if(!is.null(model_data2plot$parent_set)){
                legend("topleft", legend = names(colors_relationships),
                    col = colors_relationships, lty = 1, lwd = 2)
            }
        }else if(!is.null(model_data2plot$theta)) {
            op <- par(mar=c(3, 3, 5, 3), las = 1)
            plot(model_data2plot$theta, cex = 1.5, digits = 2, key = NULL
                , axis.col = list(side = 3)
                , xlab = "Effect of this (effector)"
                , ylab = " on this (affected)"
                , main = mod
                , mgp = c(2, 1, 0))
            par(op)
        }

        if (plot_type == "matrix"){
            plot(model_data2plot$dag_trans_mat
                , digits = 1, xlab = "", ylab = ""
                , axis.col = list(side = 1, las = 2)
                , axis.row = list(side = 2, las = 1) 
                , main = paste(mod, ": trans matrix", sep = " ")
                , cex.axis = 0.7
                , mgp = c(2, 1, 0)
                , key = NULL)
            
            if (!is.null(model_data2plot$td_trans_mat)){
                plot(model_data2plot$td_trans_mat, 
                    digits = 1, cex.axis = 0.7,
                    main = paste(mod, ": trans td matrix", sep = " "), 
                    xlab = "", ylab = "",
                    axis.col = list(side = 1, las = 2), axis.row = list(side = 2, las = 1), 
                    mgp = c(2, 1, 0), key = NULL)
            }
        } else if (plot_type == "genotypes") {
            plot_genot_fg(as_adjacency_matrix(model_data2plot$fg), simplify = FALSE)
        } else if (plot_type == "transitions") {
            plot_genot_fg(model_data2plot$transitions, data)
        }else if (plot_type == "trans_mat"){
            plot_genot_fg(model_data2plot$dag_trans_mat, data, simplify = FALSE)
        }

        if ((mod %in% c("OT")) & (plot_type == "matrix")) {
            par(mar = rep(3, 4))
            plot_sampled_genots(data)
        }
    }
    par(op1)
}