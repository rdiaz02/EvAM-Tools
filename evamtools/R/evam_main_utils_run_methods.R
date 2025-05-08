## Copyright 2022, 2025 Ramon Diaz-Uriarte, Javier Pérez de Lema Díez

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
    binary_rep <- intToBits(state)[1:num_features]
    active_features <- feature_labels[as.logical(rev(as.numeric(binary_rep)))]
    if (length(active_features) == 0) return("WT")
    return(paste(active_features, collapse = ", "))
}

## Like decode state, but using the same logic as
## hypertrapsct (e.g., see function hypertrapsct:::prob.by.time )
## They should give identical output. I am not using this function per se,
## but I am testing decode_state gives the same as this.
decode_state_ht <- function(state, num_features, feature_labels) {
    binary_string <- hypertrapsct:::DecToBin(state, num_features)
    binary_vector <- as.numeric(strsplit(binary_string, "")[[1]])
    if (sum(binary_vector) == 0) return("WT")
    return(paste(feature_labels[which(binary_vector == 1)], collapse = ", "))
}


run_HyperTraPS <- function(x, opts) {
  if (opts$seed == -1) opts$seed <- round(runif(1, 1, 1e9))
  time_out <- system.time({
        opts_rm <- which(names(opts) %in% c("nsampl",
                                            "cores") )
        opts_call <- opts[-opts_rm]
        opts_call <- c(list(obs = x), opts_call)
        out <- invisible(do.call(hypertrapsct::HyperTraPS, opts_call))
  })["elapsed"]

    num_features <- ncol(x)
    feature_labels <- colnames(x)

    states <- unique(c(out$dynamics$trans$From, out$dynamics$trans$To))
    decoded_states <- vapply(states, decode_state,
                             character(1),
                           num_features = num_features,
                           feature_labels = feature_labels)

    state_mapping <- data.frame(
        State = states,
        Features = decoded_states,
        stringsAsFactors = FALSE
    )

  ## Create a transition matrix (initialized as zeros)
  trans_mat <- Matrix(0,
                         nrow = length(states),
                         ncol = length(states),
                         sparse = TRUE)

  rownames(trans_mat) <- decoded_states
  colnames(trans_mat) <- decoded_states

    for (i in 1:nrow(out$dynamics$trans)) {
        from_state <- out$dynamics$trans$From[i]
        to_state <- out$dynamics$trans$To[i]
        probability <- out$dynamics$trans$Probability[i]
  
        from_decoded <- decode_state(from_state, num_features, feature_labels)
        to_decoded <- decode_state(to_state, num_features, feature_labels)
  
      trans_mat[from_decoded, to_decoded] <- probability
    }
    
  ## out$predicted_genotype_freqs = probs_from_trm(trans_mat)


    ## This could be inside the first loop, but this is simpler
    ## for now.
    ## Before calling it, a sanity check
    ## stopifnot(identical(colnames(x), out$featurenames))
    ## time2 <- system.time({
    ##     predicted_genotype_freqs <- probs_from_HT(out,
    ##                                               gene_names = colnames(x),
    ##                                                   nsampl = opts$nsampl,
    ##                                                   cores = opts$cores)
    ## })["elapsed"]

    ## message("HyperTraPS times: time_out = ",  round(time_out, 3),
    ##         ". time2 = ", round(time2, 3))

    return(list(time_out = time_out, ## + time2,
                out = c(primary_output = list(out),
                        trans_mat = list(trans_mat)
                        ## , predicted_genotype_freqs = list(predicted_genotype_freqs)
                        )))
}


run_BML <- function(x, opts) {
  ## FIXME: rm when settled
  ## FIXME: why do this? Just to obtain EdgeProbabilities?
  ## if (opts$rep == 0) {
  ##   opts$rep <- 1
  ## }
    
    time_out <- system.time({
        opts <- c(list(dataset = x), opts)
        out <- invisible(do.call(
            BML::bml,
            opts
        ))
    })["elapsed"]

  ## This is about right, but it is much better to provide
  ## their native output instead of try to get right the
  ## transition matrix. We could recover this code if we wanted.
  ## But NOT the transition rate matrix, as
  ## that makes no sense. To rm code later but not rm the useful things
  ## I make this dead code.

  if (FALSE) {
    out$adjacency_mat <- Matrix::Matrix(BML::adjacency_matrix(out))
    out$trans_mat <- Matrix::Matrix(BML::adjacency_matrix(out))
    ## out$td_trans_mat <- Matrix::Matrix(BML::adjacency_matrix(out))
    ## Misra does not give a transition rate matrix
    ## out$trans_rate_mat <- Matrix::Matrix(BML::adjacency_matrix(out))

    ## There are, I think, cleaner ways of getting this probs
    ## which are, actually P(g)/m_k
    ## We could hack BML and not divide, but watch out for
    ## possible partial sums.

        ## This breaks the primary output, I think. Watch out
        out$DAG$labels <- sub("Normal", "WT", out$DAG$labels)
        out$DAG$labels <- sapply(out$DAG$labels, function(label) {
            ## If label is not "WT", sort the components and join them with a comma and space
            if (label != "WT") {
                label_parts <- strsplit(label, ",")[[1]]
                sorted_label <- paste(sort(trimws(label_parts)), collapse = ", ")
                return(sorted_label)
            }
                                        # If the label is "WT", return it as is
            return(label)
        })
        out$DAG$labels <- unname(out$DAG$labels)

    for (i in seq_len(nrow(out$adjacency_mat))) {
      for (j in seq_len(ncol(out$adjacency_mat))) {
        if (out$adjacency_mat[i, j] == 0) {
          out$trans_mat[i, j] <- 0
          ## out$trans_rate_mat[i, j] <- 0
        } else {
          name = colnames(out$adjacency_mat)[j]

          idx <- which(out$DAG$labels == name)
          ## The next is most likely right
          ## the probs come from the BML code that
          ## already divides by max, max2 and max3.
          ## But the sum for all descendants can sometimes be > 1.
          ## Makes no sense but this happens with the original
          ## software too. I think those P(g)/m_k are for
          ## coloring the figures; one is not supposed to use
          ## the P(g) directly from these figures.
          out$trans_mat[i, j] = out$DAG$probs[idx]
          ## FIXME: rm later
          ## This makes no sense: EdgeProbabilites are only for pairs
          ## and singletons
          ## if (name %in% rownames(out$bootstrap$EdgeProbabilities)) {
          ##     out$trans_mat[i, j] <- mean(out$bootstrap$EdgeProbabilities[name, ])
          ##     # out$trans_mat[i, j] <- 1
          ## }
        }
      }
    }
  } ## end dead code block

  ## This is done so that all predicted genotype frequencies are 0.
  ## But I think this ain't needed. FIXME: rm when settled not needed.
  ## genotype_names = colnames(x)
  ## combinations = c("WT", unlist(lapply(1:length(genotype_names), function(n) combn(genotype_names, n, FUN = function(x) paste(x, collapse = ", "))), use.names = FALSE))
  ## freqs = rep(0, length(combinations))
  ## out$predicted_genotype_freqs <- as.data.frame(t(freqs))
  ## colnames(out$predicted_genotype_freqs) <- combinations

  ##  for (i in seq_along(out$DAG$labels)) {
  ##       label <- out$DAG$labels[i]
  ##       prob <- 0

  ##       out$predicted_genotype_freqs[label] <- prob
  ##   }

    return(list(time_out = time_out, out = c(primary_output = list(out))))
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
    } else if (method == "HyperTraPS") {
      result <- run_HyperTraPS(x, opts$hyper_traps_opts)
    } else if (method == "BML") {
        result <- run_BML(x, opts$bml_opts)
    }
    

    time_out <- result$time_out
    out <- result$out

    message(paste0("time ", method, ": ", round(time_out, 3)))
    return(c(out, elapsed_time = time_out[["elapsed"]]))
}
