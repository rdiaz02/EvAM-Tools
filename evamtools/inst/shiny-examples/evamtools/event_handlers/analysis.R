## Copyright 2022, 2025 Pablo Herrera-Nieto, Ramon Diaz-Uriarte,
## Javier Pérez de Lema Díez

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


source("utils/dataModal.R")

large_gene_number_check <- function(data, max_genes, warn = FALSE) {
    large_gene_number <- ncol(data) >= max_genes

    if (warn && large_gene_number) {
        showModal(
            dataModal(
                paste(
                    "Beware! You are analyzing data ",
                    "with 7 or more genes. ",
                    "This can take longer than usual ",
                    "and plots may be crowded. "
                ),
                ## "We recommend using top_paths options in ",
                ## "the Results' tab.",
                type = "Warning: "
            )
        )
    }

    return(large_gene_number)
}

method_validation <- function(methods, stop_if_empty = FALSE) {
    if (stop_if_empty) {
        if (is.null(methods) || (length(methods) == 1 && is.na(methods))) { 
            stop(
                "You must use at least one method ",
                "(check 'CPMs to use' under 'Advanced options ",
                "and CPMs to use')."
            )
        }
    }

    ## methods <- .ev_SHINY_dflt$cpms2run
    if (!is.null(methods)) {
        valid_methods <- unique(methods)
    }

    return(valid_methods)
}

get_mhn_args <- function(input) {
    mhn_opts = list()

    if(!is.na(input$MHN_lambda)) {
        mhn_opts$lambda <- input$MHN_lambda
    }

    return(mhn_opts)
}

get_ot_args <- function(input) {
    ot_opts <- list()
    
    if (input$OT_with_error == "TRUE") {
        ot_opts$with_errors_dist_ot <- TRUE
    } else {
        ot_opts$with_errors_dist_ot <- FALSE
    }

    return(ot_opts)
}

get_cbn_args <- function(input) {
    cbn_opts <- list(init_poset = input$CBN_init_poset,
                                 omp_threads = input$CBN_omp_threads) 

    return(cbn_opts)
}

get_hesbcn_args <- function(input) {
    hesbcn_opts <- list(
        MCMC_iter = input$HESBCN_MCMC_iter,
        reg = input$HESBCN_reg
    )

    if (!is.na(input$HESBCN_seed)) {
        hesbcn_opts$seed <- input$HESBCN_seed
    }

    return(hesbcn_opts)
}

get_oncobn_args <- function(input) {
    oncobn_opts <- list(
        model = input$OncoBN_model,
        algorithm = input$OncoBN_algorithm,
        k = input$OncoBN_k
    )

    if (!is.na(input$OncoBN_epsilon)) {
        oncobn_opts$epsilon <- input$OncoBN_epsilon
    }

    return(oncobn_opts)
}

get_mccbn_args <- function(input) {
    mccbn_opts <- list(
        model = input$MCCBN_model,
        L = input$MCCBN_L,
        sampling = input$MCCBN_sampling,
        max.iter = input$MCCBN_max_iter,
        update.step.size = input$MCCBN_update_step_size,
        tol = input$MCCBN_tol,
        max.lambda.val = input$MCCBN_max_lambda_val,
        T0 = input$MCCBN_T0,
        adap.rate = input$MCCBN_adapt_rate,
        max.iter.asa = input$MCCBN_max_iter_asa,
        neighborhood.dist = input$MCCBN_neighborhood_dist
    )

    if (!is.na(input$MCCBN_seed)) {
        mccbn_opts$seed <- input$MCCBN_seed
    }

    if (!is.na(input$MCCBN_acceptance_rate)) {
        mccbn_opts$acceptance.rate <- input$MCCBN_acceptance_rate
    }

    if (!is.na(input$MCCBN_step_size)) {
        mccbn_opts$step.size <- input$MCCBN_acceptance_rate
    }
    if (input$MCCBN_adaptive == "TRUE") {
        mccbn_opts$adaptive <- TRUE
    } else {
        mccbn_opts$adaptive <- FALSE
    }
 
    return(mccbn_opts)
}

get_hyper_traps_args <- function(input) {
    hyper_traps_opts <- list(
        length = input$HyperTraPS_length,
        kernel = input$HyperTraPS_kernel,
        walkers = input$HyperTraPS_walkers,
        samplegap = input$HyperTraPS_samplegap,
        losses = as.numeric(input$HyperTraPS_losses),
        apm_type = as.numeric(input$HyperTraPS_apm),
        sa = as.numeric(input$HyperTraPS_sa),
        sgd = as.numeric(input$HyperTraPS_sgd),
        pli = as.numeric(input$HyperTraPS_pli)
    )

    if (!is.na(input$HyperTraPS_seed)) {
        hyper_traps_opts$seed <- input$HyperTraPS_seed
    }

    return(hyper_traps_opts)
}

get_bml_args <- function(input) {
    bml_opts <- list(
        ntree = input$BML_ntree,
        threshold = input$BML_threshold,
        rep = input$BML_rep
    )

    return(bml_opts)
}

parse_opts <- function(input) {
    opts <- list()

    opts$mhn_opts <- get_mhn_args(input)
    opts$ot_opts <- get_ot_args(input)
    opts$cbn_opts <- get_cbn_args(input)
    opts$hesbcn_opts <- get_hesbcn_args(input)
    opts$oncobn_opts <- get_oncobn_args(input)
    opts$mccbn_opts <- get_mccbn_args(input)
    opts$hyper_traps_opts <- get_hyper_traps_args(input)
    opts$bml_opts <- get_bml_args(input)

    return(opts)
}

validate_data <- function(data) {
    if (ncol(data) < 1) {
        stop(
            "Your data contains less than one column (genes, events). ",
            "Do you have a single genotype? Maybe only WT?"
        )
    }
}

run_analysis <- function(data, input, disp_freqs_ret, EVAM_MAX_ELAPSED) {
    large_gene_number_check(
        data = data$data,
        max_genes = 7,
        warn = TRUE
    )

    methods <- method_validation(
        methods = input$cpm_methods,
        stop_if_empty = TRUE
    )

    shinyjs::disable("analysis")

    progress <- shiny::Progress$new()
    ## Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Running evamtools", value = 0)

    opts <- parse_opts(input)

    progress$inc(1 / 5, detail = "Setting up data")

    data2run <- evamtools:::genotypeCounts_to_data(disp_freqs_ret, e = 0)

    validate_data(data2run)

    progress$inc(2 / 5, detail = "Running CPMs")
    
    cpm_output <- R.utils::withTimeout(
        {
            evam(data2run,
                methods = methods,
                paths_max = input$return_paths_max,
                mhn_opts = opts$mhn_opts,
                ot_opts = opts$ot_opts,
                cbn_opts = opts$cbn_opts,
                hesbcn_opts = opts$hesbcn_opts,
                oncobn_opts = opts$oncobn_opts,
                mccbn_opts = opts$mccbn_opts,
                hyper_traps_opts = opts$hyper_traps_opts,
                bml_opts = opts$bml_opts,
                only_used_methods = FALSE
            )
        },
        elapsed = EVAM_MAX_ELAPSED,
        timeout = EVAM_MAX_ELAPSED,
        cpu = Inf,
        onTimeout = "silent"
    )

    if (is.null(cpm_output)) {
        stop(
            "Error running evam. ",
            "Most likely you exceeded maximum ",
            "allowed time (EVAM_MAX_ELAPSED)."
        )
    }


    sampled_from_CPMs <- NULL

    do_sampling <- input$do_sampling == "TRUE"
    if (do_sampling) {
        n_samples <- input$sample_size
        if ((is.null(n_samples)) ||
            (!is.numeric(n_samples)) ||
            (n_samples < 100)) {
            n_samples <- .ev_SHINY_dflt$cpm_samples
        }
        progress$inc(3 / 5, detail = paste("Running ", n_samples, " samples"))
        ## if (input$do_genotype_transitions) {
        ##     ## disabled when removal_note_sogt_1
        ##     sout <- c("sampled")
        ## }

        sampled_from_CPMs <-
            sample_evam(cpm_output,
                N = n_samples, methods = methods,
                output = "sampled_genotype_counts",
                ## ## disabled when removal_note_sogt_1
                ## if (input$do_genotype_transitions) {
                ##     c("sampled_genotype_counts",
                ##       "obs_genotype_transitions")
                ##                                } else {
                ##           "sampled_genotype_counts"
                ##       },
                obs_noise = input$sample_noise
            )
    }

    progress$inc(4 / 5, detail = "Post processing data")

    orig_data <- list(
        data = data2run, name = data$name,
        type = input$input2build, gene_names = data$gene_names,
        thetas = data$thetas, lambdas = data$lambdas,
        dag = data$dag, DAG_parent_set = data$DAG_parent_set
    )

    tabular_data <- evamtools:::create_tabular_data(c(cpm_output, sampled_from_CPMs))
    
    all_evam_output <- list("cpm_output" = c(cpm_output, sampled_from_CPMs)
                                      , "orig_data" = orig_data
                                      , "tabular_data" = tabular_data
                                      , "do_sampling" = do_sampling
                                        ) 
    
    progress$inc(5/5, detail = "You can see your result by going to the Results tab")

    shinyjs::enable("analysis")

    return(all_evam_output)
}