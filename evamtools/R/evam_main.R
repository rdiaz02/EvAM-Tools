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



common_preprocess <- function(x, max_cols) {
    ## ########      Preprocessing: common to all methods
    x <- df_2_mat_integer(x)
    xoriginal <- x

    x <- add_pseudosamples(x)
    ## remove.constant makes no difference IFF we add pseudosamples, as
    ## there can be no constant column when we add pseudosamples
    x <- pre_process(x, remove.constant = FALSE,
                     min.freq = 0, max.cols = max_cols)

    # Return both the processed and original x
    return(list(processed = x, original = xoriginal))
}

evam <- function(x,
                 methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN",
                             "MCCBN", "HyperTraPS", "BML"),
                 max_cols = 15,
                 cores = detectCores(),
                 paths_max = FALSE,
                 mhn_opts = list(lambda = 1/nrow(x),
                                 omp_threads = ifelse(cores > 1, 1, detectCores())),
                 ot_opts = list(with_errors_dist_ot = TRUE),
                 cbn_opts = list(
                   omp_threads = 1,
                   init_poset = "OT"
                 ),
                 hesbcn_opts = list(
                   MCMC_iter = 100000,
                   seed = NULL,
                   reg = c("bic", "aic", "loglik"),
                   silent = TRUE
                 ),
                 oncobn_opts = list(
                   model = "DBN",
                   algorithm = "DP",
                   k = 3,
                   epsilon = min(colMeans(x)/2),
                   silent = TRUE
                 ),
                 mccbn_opts = list(
                   model = "OT-CBN",
                   tmp_dir = NULL,
                   addname = NULL,
                   silent = TRUE,
                   L = 100,
                   sampling = c("forward", "add-remove",
                                "backward", "bernoulli", "pool"),
                   max.iter = 100L,
                   update.step.size = 20L,
                   tol = 0.001,
                   max.lambda.val = 1e+06,
                   T0 = 50,
                   adap.rate = 0.3,
                   acceptance.rate = NULL,
                   step.size = NULL,
                   max.iter.asa = 10000L,
                   neighborhood.dist = 1L,
                   adaptive = TRUE,
                   thrds = 1L,
                   verbose = FALSE,
                   seed = NULL),
                 hyper_traps_opts = list(
                   initialstates = NULL,
                   priors = NULL,
                   starttimes = NULL,
                   endtimes = NULL,
                   length = 3,
                   kernel = 5,
                   samplegap = -1,
                   losses = 0,
                   apm_type = 0,
                   sa = 0,
                   sgd = 0,
                   sgd_scale = 0.01,
                   seed = 1,
                   outputinput = 0,
                   regularise = 0,
                   penalty = 0,
                   lasso = 0,
                   model = 2,
                   pli = 0,
                   walkers = 200,
                   full_analysis = 1,
                   limited_output = 0,
                   output_transitions = 1,
                   samples_per_row = 10,
                   featurenames = NULL # Replace NULL with the actual character vector if available
                 ),
                 bml_opts = list(
                   ntree = 100, ## what they used in the paper
                   ## 0.3 is what they said, in their README, they used in the ms.
                   threshold = 0.3,
                   rep = 10 ## In the paper, 1000. That is way too many for webapp
                 ),
                 only_used_methods = TRUE
                 ) {

    preprocessed_x <- common_preprocess(x, max_cols)
    x <- preprocessed_x$processed
    xoriginal <- preprocessed_x$original

    ## Dealing with the multiple spellings of HyperTraPS
    methods <- gsub("HyperTraps", "HyperTraPS", methods, fixed = TRUE)
    methods <- gsub("Hypertraps", "HyperTraPS", methods, fixed = TRUE)
    methods <- gsub("hypertraps", "HyperTraPS", methods, fixed = TRUE)


    ## Get default parametrs
    default_opts <- create_default_opts(x)

    opts <- list(
        mhn_opts = fill_args_default(mhn_opts, default_opts$mhn_opts),
        ot_opts = fill_args_default(ot_opts, default_opts$ot_opts),
        cbn_opts = fill_args_default(cbn_opts, default_opts$cbn_opts),
        hesbcn_opts = fill_args_default(hesbcn_opts, default_opts$hesbcn_opts),
        oncobn_opts = fill_args_default(oncobn_opts, default_opts$oncobn_opts),
        mccbn_opts = fill_args_default(mccbn_opts, default_opts$mccbn_opts),
        hyper_traps_opts = fill_args_default(hyper_traps_opts, default_opts$hyper_traps_opts),
        bml_opts = fill_args_default(bml_opts, default_opts$bml_opts)
    ) 



    check_cbn_opts_init_poset(opts$cbn_opts$init_poset)

    if ("MCCBN" %in% methods) {
        MCCBN_INSTALLED <- requireNamespace("mccbn", quietly = TRUE)
        
        if (!MCCBN_INSTALLED) {
            warning("MCCBN method requested, but mccbn packaged not installed. ",
                    "Removing MCCBN from list of requested methods.")
            methods <- setdiff(methods, "MCCBN")
        } else {
            stopifnot(opts$mccbn_opts$model %in% c("OT-CBN", "H-CBN2"))
        }
    }

    if ("OncoBN" %in% methods) {
        stopifnot(opts$oncobn_opts$model %in% c("DBN", "CBN"))
    }


    if ("HyperTraPS" %in% methods) {
        if (is.null(opts$hyper_traps_opts$featurenames)) {
            opts$hyper_traps_opts$featurenames <- colnames(x)
        }
    }

    # Sanity check of gene names
    sanity_check_colnames(colnames(x), TRUE)

    ## Sanity checks of methods
    methods <- sanity_check_methods(methods)


    sanity_check_data_dimensions(x)

    all_out <- mclapply(methods, run_method, x=x, opts=opts, mc.cores = cores)
    names(all_out) <- methods

    get_output <- function(method, component) {
        if (!exists(method, all_out)) return(NA)
        if (!exists(component, all_out[[method]])) return(NA)
        return(all_out[[method]][[component]])
    }

    get_all_method_output <- function(method) {
        if (!exists(method, all_out)) return(NA)
        return(all_out[[method]])
    }

    ## To avoid repeating code
    get_paths_max <- function(method) {
        if (paths_max) {
            trans_mat_name <- ifelse(method == "MHN",
                                     "transitionMatrixCompExp",
                                     "trans_mat_genots")
            trans_mat <- get_output(method, trans_mat_name)
            if ((length(trans_mat) == 1) && is.na(trans_mat)) return(NA)
            return(paths_probs_2_df(trans_mat_2_paths_probs(trans_mat),
                                    order = "prob"))
        } else {
            return(NA)
        }
    }

    tmpr <- list(
        OT_model = get_output("OT", "edges"),
        OT_f_graph = get_output("OT", "weighted_fgraph"),
        OT_trans_mat = get_output("OT", "trans_mat_genots"),
        OT_predicted_genotype_freqs = get_output("OT",
                                                 "predicted_genotype_freqs"),
        OT_eps = get_output("OT", "eps"),
        OT_fit = get_output("OT", "ot_fit"),
        OT_paths_max = get_paths_max("OT"),
        OT_elapsed_time = get_output("OT", "elapsed_time"),
       
        CBN_model = get_output("CBN", "edges"),
        CBN_trans_rate_mat = get_output("CBN", "weighted_fgraph"),
        CBN_trans_mat = get_output("CBN", "trans_mat_genots"),
        CBN_td_trans_mat = get_output("CBN", "td_trans_mat"),
        CBN_predicted_genotype_freqs = get_output("CBN",
                                                  "predicted_genotype_freqs"),
        CBN_paths_max = get_paths_max("CBN"),
        CBN_elapsed_time = get_output("CBN", "elapsed_time"),
        
        MCCBN_model = get_output("MCCBN", "edges"),
        MCCBN_trans_rate_mat = get_output("MCCBN", "weighted_fgraph"),
        MCCBN_trans_mat = get_output("MCCBN", "trans_mat_genots"),
        MCCBN_td_trans_mat = get_output("MCCBN", "td_trans_mat"),
        MCCBN_predicted_genotype_freqs = get_output("MCCBN",
                                                    "predicted_genotype_freqs"),
        MCCBN_paths_max = get_paths_max("MCCBN"),
        MCCBN_elapsed_time = get_output("MCCBN", "elapsed_time"),
        
        MHN_theta = get_output("MHN", "theta"),
        MHN_trans_rate_mat = get_output("MHN", "transitionRateMatrix"),
        MHN_trans_mat = get_output("MHN", "transitionMatrixCompExp"),
        MHN_td_trans_mat = get_output("MHN", "transitionMatrixTimeDiscretized"),
        MHN_exp_theta = exp(get_output("MHN", "theta")),
        MHN_predicted_genotype_freqs = get_output("MHN",
                                                  "predicted_genotype_freqs"),
        MHN_paths_max = get_paths_max("MHN"),
        MHN_elapsed_time = get_output("MHN", "elapsed_time"),
        
        OncoBN_model = get_output("OncoBN", "edges"),
        OncoBN_likelihood = get_output("OncoBN", "likelihood"),
        OncoBN_f_graph = get_output("OncoBN", "weighted_fgraph"), 
        OncoBN_trans_mat = get_output("OncoBN", "trans_mat_genots"),
        OncoBN_predicted_genotype_freqs = get_output("OncoBN",
                                                     "predicted_genotype_freqs"),
        OncoBN_fitted_model = get_output("OncoBN", "model"),
        OncoBN_epsilon = get_output("OncoBN", "epsilon"),
        OncoBN_parent_set = get_output("OncoBN", "parent_set"),
        OncoBN_fit = get_output("OncoBN", "fit"),
        OncoBN_paths_max = get_paths_max("OncoBN"),
        OncoBN_elapsed_time = get_output("OncoBN", "elapsed_time"),
        
        HESBCN_model = get_output("HESBCN", "edges"),
        HESBCN_parent_set = get_output("HESBCN", "parent_set"),
        HESBCN_trans_rate_mat = get_output("HESBCN", "weighted_fgraph"),
        HESBCN_trans_mat = get_output("HESBCN", "trans_mat_genots"),
        HESBCN_td_trans_mat = get_output("HESBCN", "td_trans_mat"),
        HESBCN_predicted_genotype_freqs = get_output("HESBCN",
                                                     "predicted_genotype_freqs"),
        HESBCN_command = get_output("HESBCN", "command"),
        HESBCN_paths_max = get_paths_max("HESBCN"),
        HESBCN_elapsed_time = get_output("HESBCN", "elapsed_time"),

        ## FIXME: For HyperTraPS and BML we return
        ## pieces, and then all output. So some things
        ## are returned in two places.
        HyperTraPS_model = get_output("HyperTraPS", "model"),
        HyperTraPS_edges = get_output("HyperTraPS", "edges"),
        HyperTraPS_posterior_samples = get_output("HyperTraPS", "posterior.samples"),
        HyperTraPS_dynamics = get_output("HyperTraPS", "dynamics"),
        HyperTraPS_best = get_output("HyperTraPS", "best"),
        HyperTraPS_lik_traces = get_output("HyperTraPS", "lik.traces"),
        HyperTraPS_bubbles = get_output("HyperTraPS", "bubbles"),
        HyperTraPS_routes = get_output("HyperTraPS", "routes"),
        HyperTraPS_times = get_output("HyperTraPS", "times"),
        HyperTraPS_timediffs = get_output("HyperTraPS", "timediffs"),
        HyperTraPS_timehists = get_output("HyperTraPS", "timehists"),
        HyperTraPS_trans_mat = get_output("HyperTraPS", "td_trans_mat"),
        HyperTraPS_predicted_genotype_freqs = get_output("HyperTraPS", "predicted_genotype_freqs"),
        HyperTraPS_elapsed_time = get_output("HyperTraPS", "elapsed_time"),
        HyperTraPS_post = get_all_method_output("HyperTraPS"),

    BML_model = get_output("BML", "model"),
    BML_trans_mat = NA,
    BML_predicted_genotype_freqs = NA,
        BML_elapsed_time = get_output("BML", "elapsed_time"),
    BML_bootstrap = ifelse(exists("BML", all_out), opts$bml_opts$rep, NA),
        BML_output = get_all_method_output("BML"),

    methods = methods,
        original_data = xoriginal,
        analyzed_data = x,
        genotype_id_ordered =
            stats::setNames(1:(2^ncol(x)),
                            genotypes_standard_order(colnames(x))),
        all_options = list(
            mhn_opts = opts$mhn_opts,
            ot_opts = opts$ot_opts,
            cbn_opts = opts$cbn_opts,
            hesbcn_opts = opts$hesbcn_opts,
            oncobn_opts = opts$oncobn_opts,
            mccbn_opts = opts$mccbn_opts,
            hyper_traps_opts = opts$hyper_traps_opts
            )
    )
    
    if (only_used_methods) {
        tmpr_rm <- lapply(tmpr, function(x) (length(x) == 1) && (is.na(x)))
        tmpr <- tmpr[which(!unlist(tmpr_rm))]
    }

    return(tmpr)
}