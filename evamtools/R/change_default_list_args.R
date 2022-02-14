
## list, list -> list
## if any entry en defaults is not in inputl, add to inputl
## I us it to add defaults (possibly not explicitly passed) to the
##      list inputl, if the defaults have no value
fill_args_default <- function(inputl, defaults) {

    not_passed <- vapply(names(mccbn_hcbn2_opts),
                         function(x) !exists(x, where = l3),
                         TRUE)

    inputl <- c(input, defaults[not_passed])
    return(inputl)
}
