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
                 methods = c("CBN", "OT", "HESBCN", "MHN", "OncoBN", "MCCBN", "HyperTraps", "BML"),
                 max_cols = 15,
                 cores = detectCores(),
                 paths_max = FALSE,
                 mhn_opts = list(),
                 ot_opts = list(), 
                 cbn_opts = list(),
                 hesbcn_opts = list(),
                 oncobn_opts = list(),
                 mccbn_opts = list(),
                 hyper_traps_opts = list(),
                 bml_opts = list(),
                 only_used_methods = FALSE
) {
    preprocessed_x <- common_preprocess(x, max_cols)
    x <- preprocessed_x$processed
    xoriginal <- preprocessed_x$original

    # Get default parametrs
    default_opts <- evamtools:::create_default_opts(x)

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



    evamtools:::check_cbn_opts_init_poset(opts$cbn_opts$init_poset)

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

    if ("HyperTraps" %in% methods) {
        if (is.null(opts$hyper_traps_opts$featurenames)) {
            opts$hyper_traps_opts$featurenames <- colnames(x)
        }
    }

    # Sanity check of gene names
    evamtools:::sanity_check_colnames(colnames(x), TRUE)

    ## Sanity checks of methods 
    methods <- evamtools:::sanity_check_methods(methods)


    evamtools:::sanity_check_data_dimensions(x)

    all_out <- mclapply(methods, evamtools:::run_method, x=x, opts=opts, mc.cores = cores)
    names(all_out) <- methods

    get_output <- function(method, component) {
        if (!exists(method, all_out)) return(NA)
        if (!exists(component, all_out[[method]])) return(NA)
        return(all_out[[method]][[component]])
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
        MCCBN_td_trans_mat = get_output("CBN", "td_trans_mat"),
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

        HyperTraps_model = get_output("HyperTraps", "model"),
        HyperTraps_edges = get_output("HyperTraps", "edges"),
        HyperTraps_posterior_samples = get_output("HyperTraps", "posterior.samples"),
        HyperTraps_dynamics = get_output("HyperTraps", "dynamics"),
        HyperTraps_best = get_output("HyperTraps", "best"),
        HyperTraps_lik_traces = get_output("HyperTraps", "lik.traces"),
        HyperTraps_bubbles = get_output("HyperTraps", "bubbles"),
        HyperTraps_routes = get_output("HyperTraps", "routes"),
        HyperTraps_times = get_output("HyperTraps", "times"),
        HyperTraps_timediffs = get_output("HyperTraps", "timediffs"),
        HyperTraps_timehists = get_output("HyperTraps", "timehists"),
        HyperTraps_trans_mat = get_output("HyperTraps", "td_trans_mat"),
        HyperTraps_predicted_genotype_freqs = get_output("HyperTraps", "predicted_genotype_freqs"),
        HyperTraps_elapsed_time = get_output("HyperTraps", "elapsed_time"),
        HyperTraps_post = all_out$HyperTraps,

        BML_trans_mat = get_output("BML", "trans_mat"),
        BML_predicted_genotype_freqs = get_output("BML", "predicted_genotype_freqs"),
        BML_trans_rate_mat = get_output("BML", "trans_rate_mat"),
        BML_elapsed_time = get_output("BML", "elapsed_time"),
        BML_output = all_out$BML,


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