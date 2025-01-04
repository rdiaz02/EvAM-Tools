# Check if any column names contain a comma
check_comma <- function(x) {
    gn_comma <- stringi::stri_count_fixed(x, ",")
    return(!any(gn_comma)) # TRUE if no commas, FALSE if any commas
}

# Check if any column names contain a backslash
check_backslash <- function(x) {
    gn_backslash <- stringi::stri_count_fixed(x, "\\")
    return(!any(gn_backslash)) # TRUE if no backslashes, FALSE if any backslashes
}

# Check if any column names contain a space
check_space <- function(x) {
    gn_space <- stringi::stri_count_regex(x, "[\\s]")
    return(!any(gn_space)) # TRUE if no spaces, FALSE if any spaces
}

# Check if any column name is "WT"
check_WT <- function(x) {
    return(!any(x == "WT")) # TRUE if no "WT", FALSE if "WT" is found
}

check_valid_methods <- function(methods) {
    accepted_methods <- c("OT", "OncoBN", "CBN", "MCCBN", "MHN", "HESBCN", "HyperTraps")
    not_valid_methods <- which(!(methods %in% accepted_methods))
    
    return(not_valid_methods)
}
check_cbn_opts_init_poset <- function(init_poset) {
    if (!(init_poset %in% c("OT", "linear"))) {
        stop("CBN's init_poset must be one of OT or linear. ",
                " Custom not allowed in call from evam.")
    }
}


sanity_check_colnames <- function(x, stop_on_error = TRUE) {
    valid <- TRUE # Assume the input is valid initially
    error_messages <- character() # Initialize a character vector to store error messages

    # Check for commas
    if (!check_comma(x)) {
        error_messages <- c(error_messages, "At least one of your gene names has a comma. That is not allowed")
        valid <- FALSE
    }

    # Check for backslashes
    if (!check_backslash(x)) {
        error_messages <- c(error_messages, "At least one of your gene names has a backslash. That is not allowed")
        valid <- FALSE
    }

    # Check for spaces
    if (!check_space(x)) {
        error_messages <- c(error_messages, "At least one of your gene names has a space. That is not allowed")
        valid <- FALSE
    }

    # Check for "WT"
    if (!check_WT(x)) {
        error_messages <- c(error_messages, "One of your genes is called WT. That is not allowed")
        valid <- FALSE
    }

    # If there were errors, return the messages
    if (!valid) {
        # Concatenate all error messages into one string
        error_message <- paste(error_messages, collapse = "\n")

        if (stop_on_error) {
            stop(error_message) # Stop execution if stop_on_error is TRUE
        }
    }

    # If all checks pass, return TRUE
    return(TRUE)
}

sanity_check_methods <- function(methods) {
    methods <- unique(methods)

    not_valid_methods <- check_valid_methods(methods)

    if (length(not_valid_methods)) {
        warning("Method(s) ",
                paste(methods[not_valid_methods], sep = ", ", collapse = ", "),
                " not among the available methods.",
                " Ignoring the invalid method.")
        methods <- methods[-not_valid_methods]
    }

    if ("MCCBN" %in% methods) {
        MCCBN_INSTALLED <- requireNamespace("mccbn", quietly = TRUE)
        if (!MCCBN_INSTALLED) {
            warning("MCCBN method requested, but mccbn packaged not installed. ",
                    "Removing MCCBN from list of requested methods.")
            methods <- setdiff(methods, "MCCBN")
        }
    }

    if (length(methods) == 0) {
        stop("No valid methods given.")
    }

    return(methods)
}

sanity_check_data_dimensions <- function(x) {
    if (ncol(x) < 2) {
        stop("Fewer than 2 columns in the data. ",
             "There must be at least two genes ",
             "and two different genotypes to run evam ",
             "(and remember that any genes that are ",
             "completely aliased, i.e., indistinguishable, ",
             "because they have identical patterns ",
             "---identical columns in the data matrix--- ",
             "are regarded as a single gene).")
    }
}
