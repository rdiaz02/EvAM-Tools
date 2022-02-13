
## list, list -> list

## If any entry en defaults is not in inputl, add to inputl
## If inputl has an name component not in defaults, give warning.

## I use it to add defaults (possibly not explicitly passed) to the
##      list inputl, if the inputl has no value


fill_args_default <- function(inputl, defaults) {

    not_valid <- which(!(names(inputl) %in% names(defaults)))

    if(length(not_valid)) {
        name_input <- deparse(substitute(inputl))
        not_valid_opt <- names(inputl)[not_valid]
        warning("Option(s) ", paste(not_valid_opt, sep = ", ", collapse = ", "),
                " in argument ", name_input, " is(are) invalid.",
                " They will be ignored.")
        inputl <- inputl[-not_valid]
    }
        
    not_passed <- vapply(names(defaults),
                         function(x) !exists(x, where = inputl),
                         TRUE)
    inputl <- c(inputl, defaults[not_passed])
    return(inputl)
}

## This can introduce bugs
## ## wrap the above, expecting the defaults to be called d_name
## ## For this to work, the default arguments must start with "d_"
## wrap_fill_args_default <- function(name) {
##     name0 <- paste0("d_", deparse(substitute(name)))
##     tmp <- fill_args_default(name, get(name0, envir = parent.frame(n = 1)))
##     return(tmp)
## }
