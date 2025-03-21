run_MHN <- function(x, opts) {
    RhpcBLASctl::omp_set_num_threads(opts$omp_threads)
    time_out <- system.time({
        out <- do_MHN2(x, lambda = opts$lambda)
        out <- c(out,
            predicted_genotype_freqs =
                list(probs_from_trm(out$transitionRateMatrix))
        )
    })["elapsed"]

    return(list(time_out = time_out, out = out))
}

run_HESBCN <- function(x, opts) {
    time_out <- system.time({
        out <- do_HESBCN(x,
            MCMC_iter = opts$MCMC_iter,
            seed = opts$seed,
            silent = opts$silent,
            reg = opts$reg
        )
        out <- c(out, cpm2tm(out))
        out <- c(out,
            td_trans_mat =
                trans_rate_to_trans_mat(out[["weighted_fgraph"]],
                    method = "uniformization"
                )
        )
        out <- c(out,
            predicted_genotype_freqs =
                list(probs_from_trm(out$weighted_fgraph))
        )
    })["elapsed"]

    return(list(time_out = time_out, out = out))
}

run_CBN <- function(x, opts) {
    time_out <- system.time({
        out <- try(cbn_proc(x,
            addname = "tmpo",
            init.poset = opts$init_poset,
            nboot = 0,
            parall = TRUE,
            omp_threads = opts$omp_threads
        ))
        out <- c(out, cpm2tm(out))
        out <- c(out,
            td_trans_mat =
                trans_rate_to_trans_mat(out[["weighted_fgraph"]],
                    method = "uniformization"
                )
        )
        out <- c(out,
            predicted_genotype_freqs =
                list(probs_from_trm(out$weighted_fgraph))
        )
    })["elapsed"]

    return(list(time_out = time_out, out = out))
}

run_MCCBN <- function(x, opts) {
    if (opts$model == "OT-CBN") {
        time_out <-
            system.time(out <- try(do_MCCBN_OT_CBN(x)))["elapsed"]
    } else if (opts$model == "H-CBN2") {
        mccbn_hcbn2_opts_2 <- opts
        mccbn_hcbn2_opts_2$model <- NULL
        time_out <-
            system.time(
                out <- try(
                    do_MCCBN_HCBN2(x,
                        mccbn_hcbn2_opts = mccbn_hcbn2_opts_2
                    )
                )
            )["elapsed"]
    }
    time_out2 <- system.time({
        out <- c(out, cpm2tm(out))
        out <- c(out,
            td_trans_mat =
                trans_rate_to_trans_mat(out[["weighted_fgraph"]],
                    method = "uniformization"
                )
        )
        out <- c(out,
            predicted_genotype_freqs =
                list(probs_from_trm(out$weighted_fgraph))
        )
    })["elapsed"]
    time_out <- time_out + time_out2

    return(list(time_out = time_out, out = out))
}

run_OT <- function(x, opts) {
    time_out <- system.time({
        out <- try(
            suppressMessages(
                ot_proc(x,
                    nboot = 0,
                    distribution.oncotree = TRUE,
                    with_errors_dist_ot = opts$with_errors_dist_ot
                )
            )
        )
        out <- c(out, cpm2tm(out))
    })["elapsed"]

    return(list(time_out = time_out, out = out))
}

run_OncoBN <- function(x, opts) {
    time_out <- system.time({
        out <- do_OncoBN(x,
            model = opts$model,
            algorithm = opts$algorithm,
            k = opts$k,
            epsilon = opts$epsilon,
            silent = opts$silent
        )
        out <- c(out, cpm2tm(out))
    })["elapsed"]

    return(list(time_out = time_out, out = out))
}

decode_state <- function(state, num_features, feature_labels) {
    binary_rep <- intToBits(state)[1:num_features]  # Convert to binary
    active_features <- feature_labels[as.logical(rev(as.numeric(binary_rep)))]  # Extract active features
    if (length(active_features) == 0) return("WT")  # If no active features, return "WT"
    return(paste(active_features, collapse = ", "))  # Return active features
}


run_HyperTraps <- function(x, opts) {
    time_out <- system.time({
        opts <- c(list(obs = x), opts)
        out <- invisible(do.call(
            hypertrapsct::HyperTraPS,
            opts
        ))
    })["elapsed"]

    num_features <- ncol(x)
    feature_labels <- colnames(x)

    states <- unique(c(out$dynamics$trans$From, out$dynamics$trans$To))
    decoded_states <- sapply(states, decode_state, num_features = num_features, feature_labels = feature_labels)

    state_mapping <- data.frame(
        State = states,
        Features = decoded_states,
        stringsAsFactors = FALSE
    )

    # Create a transition matrix (initialized as zeros)
    trans_matrix <- Matrix(0, nrow = length(states), ncol = length(states), sparse = TRUE)

    for (i in 1:nrow(out$dynamics$trans)) {
        from_state <- out$dynamics$trans$From[i]
        to_state <- out$dynamics$trans$To[i]
        probability <- out$dynamics$trans$Probability[i]
  
        from_decoded <- decode_state(from_state, num_features, feature_labels)
        to_decoded <- decode_state(to_state, num_features, feature_labels)
  
        trans_matrix[from_state + 1, to_state + 1] <- probability
    }

    rownames(trans_matrix) <- decoded_states
    colnames(trans_matrix) <- decoded_states
    
    out$predicted_genotype_freqs = probs_from_trm(trans_matrix)
    out$td_trans_mat <- trans_matrix

    return(list(time_out = time_out, out = out))
}


run_BML <- function(x, opts) {
    time_out <- system.time({
        opts <- c(list(dataset = x), opts)
        out <- invisible(do.call(
            BML::bml,
            opts
        ))
    })["elapsed"]

    out$DAG$labels <- sub("Normal", "WT", out$DAG$labels)
    out$DAG$labels <- sapply(out$DAG$labels, function(label) {
        # If label is not "WT", sort the components and join them with a comma and space
        if (label != "WT") {
        label_parts <- strsplit(label, ",")[[1]]
        sorted_label <- paste(sort(trimws(label_parts)), collapse = ", ")
        return(sorted_label)
        }
        # If the label is "WT", return it as is
        return(label)
    }) 

    unname(out$DAG$labels)

    out$trans_mat <- Matrix::Matrix(BML::adjacency_matrix(out))
    out$trans_rate_mat <- Matrix::Matrix(BML::adjacency_matrix(out))
 
    for (i in seq_len(nrow(out$trans_mat))) {
        for (j in seq_len(ncol(out$trans_mat))) {
            if (out$trans_mat[i, j] == 1) {
                # Scale the transition rate by the probability for node i
                out$trans_rate_mat[i, j] <- out$DAG$probs[i]
            } else {
                # If no edge exists, transition rate remains 0
                out$trans_rate_mat[i, j] <- 0
            }
        }
    }

    genotype_names = colnames(x)
    combinations = c("WT", unlist(lapply(1:length(genotype_names), function(n) combn(genotype_names, n, FUN = function(x) paste(x, collapse = ", "))), use.names = FALSE))
    freqs = rep(0, length(combinations))
    out$predicted_genotype_freqs <- as.data.frame(t(freqs))
    colnames(out$predicted_genotype_freqs) <- combinations

   for (i in seq_along(out$DAG$labels)) {
        label <- out$DAG$labels[i]
        prob <- 0
  
        out$predicted_genotype_freqs[label] <- prob
    } 

    return(list(time_out = time_out, out = out))
}

run_method <- function(method, x, opts) {
    if (method == "MHN") {
        result <- run_MHN(x, opts$mhn_opts)
    } else if (method == "HESBCN") {
        result <- run_HESBCN(x, opts$hesbcn_opts)
    } else if (method == "CBN") {
        result <- run_CBN(x, opts$cbn_opts)
    } else if (method == "MCCBN") {
        result <- run_MCCBN(x, opts$mccbn_opts)
    } else if (method == "OT") {
        result <- run_OT(x, opts$ot_opts)
    } else if (method == "OncoBN") {
        result <- run_OncoBN(x, opts$oncobn_opts)
    } else if (method == "HyperTraps") {
        result <- run_HyperTraps(x, opts$hyper_traps_opts)
    } else if (method == "BML") {
        result <- run_BML(x, opts$bml_opts)
    }
    

    time_out <- result$time_out
    out <- result$out

    message(paste0("time ", method, ": ", round(time_out, 3)))
    return(c(out, elapsed_time = time_out[["elapsed"]]))
}
