# default-method-params.R

# Default options for each method

create_default_opts <- function(x) {
    list(
        mhn_opts = list(
            lambda = 1 / nrow(x),
            omp_threads = ifelse(detectCores() > 1, 1, detectCores())
        ),

        ot_opts = list(
            with_errors_dist_ot = TRUE
        ),

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
            epsilon = min(colMeans(x) / 2),
            silent = TRUE
        ),

        mccbn_opts = list(
            model = "OT-CBN",
            tmp_dir = NULL,
            addname = NULL,
            silent = TRUE,
            L = 100,
            sampling = c("forward", "add-remove", "backward", "bernoulli", "pool"),
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
            seed = NULL
        ),

        hyper_traps_opts = list(
            # obs , is x in evam
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
        )
    )
}
